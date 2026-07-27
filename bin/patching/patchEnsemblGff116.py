# patchEnsemblGff116.py
# 
# Patches GFF3 genome annotation files from Ensembl for version 116 (and possibly others - we'll see).
# 
# 1. Patching the files (adding "projection_parent_gene=").
# The problem is the files from Ensembl do not contain project_parent_gene attributes for 
# protein coding genes (but it does for non-coding genes). The projection parent is the key 
# piece used by multiple downstream processes for linking the Ensembl gene models with MGI genes. 
# This script patches the files by adding the missing attribute, inferring the Ensembl ID from
# other information available in the gene's col 9 attributes:
# - NCBI Gene ID: The feature's description attribute may contain a substring like this:
#   "[Source:NCBI gene (formerly Entrezgene);Acc:226304]"
# - MGI ID: The description may contain a substring like this:
#   "[Source:MGI Symbol;Acc:MGI:3801960]"
# - Mouse symbol: The feature may have a Name attribute containing a mouse symbol.
# The script uses this information to try to derive the projection_parent_gene using MGI data.
# If successful, it adds the missing projection parent attribute to the output.
#
# 2. Limiting association counts.
# Another (separate) issue with v116 data is that there are cases of multiple strain genes (sometimes many) 
# within a single genome that refer to the same projection parent gene. The sequence group discussed this and 
# decided to cap (at 3) the max number of strain genes (in one genome) that can be associated with 
# one projection parent. If there are more than the limit, then none of the strain gene is associated with 
# the projection parent. 
# 
# The script imposes the limit by removing project_parent_genes attributes as needed. Specifically, it 
# first counts the number of strain genes per projection parent ID and gets the list of IDs over the limit.
# Then it removes all occurrences projection_parent_gene attributes that specify one of these IDs.
#
# 3. Old MGP ids (adding "mgp=").
# Prior to release 114, a strain's genome features were assigned IDs like "MGP_AJ_123456", ie "MGP", an 
# abbreviation of strain name, and a number. Starting in release 114, Ensembl changed the ID scheme
# to look like any other Ensembl mouse gene ID, like "ENSMUSG000001234". Unfortunately, there's no mapping
# between the new and the old IDs.
#
# The patching script uses data archived from MGI (see WTS2-1665) to recover the MGPs for many genes. 
# The archived data is a table, one row per strain gene, with columns that include the
# MGP id and the Ensembl id of the associated MGI gene. The Ensembl id is used to match against the 
# projection_parent_genes (either supplied or inferred in part 1). Matched MGP ids are added as a new 
# attribute ("mgp") in column 9.
#

#
import os
import sys
import re
from db import sql
import gff3lite
from urllib.request import urlopen
import argparse


PPG = 'projection_parent_gene'
MGP = 'mgp'
STRAIN = ''
PATCH_ARCHIVED_MGP_IDS=os.environ['PATCH_ARCHIVED_MGP_IDS']

class Patcher :
    def __init__ (self):
        self.LIMIT = 3
        self.ifd = sys.stdin
        self.ofd = sys.stdout
        self.PPG2count = {}
        # [Source:NCBI gene (formerly Entrezgene);Acc:226304]
        self.RE1 = re.compile(r'\[Source:NCBI.*Acc:([0-9]+)\]')
        # [Source:MGI Symbol;Acc:MGI:3801960] 
        self.RE2 = re.compile(r'\[Source:MGI.*(MGI:[0-9]+)\]')


    def getArgs (self) :
        parser = argparse.ArgumentParser(description='Patch Ensemble gff3 files by adding ' +
            'projection_parent_gene attributes where possible, and apply max association limit.')
        parser.add_argument('-i', '--input', metavar='FILE',
            help='Input gff3 file. If not specified, reads from stdin.')
        parser.add_argument('-o', '--output', metavar='FILE',
            help='Output gff3 file. If not specified, writes to stdout.')
        parser.add_argument('-L', '--limit', metavar='INT', type=int, default=self.LIMIT,
            help='Association count limit.')

        args = parser.parse_args()
        if args.input:
            self.ifd = open(args.input, 'r')
        if args.output:
            self.ofd = open(args.output, 'w')
        self.LIMIT = args.limit
        return args


    def log (self, s) :
        sys.stderr.write(s + '\n')

    def incPPGCount(self, ppgs) :
        for ppg in ppgs:
            ppg = re.sub(r'[.][0-9]*$', '', ppg)
            ppgCount = self.PPG2count.setdefault(ppg,0)
            self.PPG2count[ppg] = ppgCount + 1

    def getMgiGeneModelIds (self):
        symbol2mgi = {}
        entrez2mgi = {}
        mgi2ensembl = {}
        q1 = '''
            SELECT m.symbol, a.accid
            FROM mrk_marker m, acc_accession a
            WHERE m._organism_key = 1
            AND m._marker_status_key = 1
            AND m._marker_key = a._object_key
            AND a._mgitype_key = 2
            AND a._logicaldb_key = 1
            AND a.preferred = 1
            '''
        for r in sql(q1):
            symbol2mgi[r['symbol']] = r['accid']

        q2 = '''
            SELECT ma.accid as mgiid, me.accid, me._logicaldb_key
            FROM acc_accession ma, acc_accession me
            WHERE ma._mgitype_key = 2
            AND ma._logicaldb_key = 1
            AND ma.preferred = 1
            AND ma._object_key = me._object_key
            AND me._mgitype_key = 2
            AND me._logicaldb_key in (59,60)
            AND me.preferred = 1
            '''
        for r in sql(q2):
            if r['_logicaldb_key'] == 59:
                entrez2mgi[r['accid']] = r['mgiid'] 
            else:
                mgi2ensembl.setdefault(r['mgiid'],[]).append(r['accid'])

        return symbol2mgi, entrez2mgi, mgi2ensembl

    def getMGPids (self, fname) :
        ensembl2mgp = {}
        with open(fname, 'r') as fd:
            for line in fd:
                if '|' in line and 'strain' not in line:
                    fields = line.strip().split('|')
                    strain = fields[2].strip().replace('/','').lower()
                    mgpid = fields[3].strip()
                    ensemblid = fields[4].strip()
                    ensembl2mgp.setdefault(strain,{}).setdefault(ensemblid, []).append(mgpid)
        self.log(str(list(ensembl2mgp['aj'].items())[:25]))
        return ensembl2mgp


    def process(self, entrez2mgi, mgi2ensembl, symbol2mgi, ensembl2mgp) :
        global STRAIN
        results = []
        for line in self.ifd:
            # comment lines get printed as is
            if line.startswith('#') :
                if line.startswith('#!genome-build '):
                    #!genome-build  A_J_v3
                    gb = line.strip().split()[1]
                    i = gb.rindex('_')
                    STRAIN = gb[:i].replace('_','').lower()
                    self.log('STRAIN = ' + STRAIN)
                results.append(line)
                continue
            # Any feature line that is not a top level feature, or is top level but has no projection parent,
            # gets printed as is.
            f = gff3lite.parseLine(line)
            attrs = f[8] # get a handle on the column 9 attributes
            if PPG in attrs:
                # strip the version number
                attrs[PPG] = [attrs[PPG].split('.')[0]]
                self.incPPGCount(attrs[PPG])
            elif 'Parent' in attrs or 'ID' not in attrs:
                # If non-top-level feature, or biological region (or other exotic) feature, just print it
                pass
            # if here, top level feature with no PPG. Try to infer one
            else:
                if 'description' in attrs:
                    description = attrs['description']
                    m = self.RE1.search(description)
                    m2 = self.RE2.search(description)
                    if m:
                        # Has an Entrez Gene ID in Source description
                        # Lookup to get MGI number, then lookup MGI number to get Ensembl ID
                        entrezId = m.group(1)
                        if entrezId in entrez2mgi:
                            mgiId = entrez2mgi[entrezId]
                            if mgiId in mgi2ensembl:
                                attrs[PPG] = sorted(mgi2ensembl[mgiId])
                                self.incPPGCount(attrs[PPG])
                    elif m2:
                        # Has an MGI Gene ID in Source description
                        # Lookup MGI id to get Ensembl
                        mgiId = m2.group(1)
                        if mgiId in mgi2ensembl:
                            attrs[PPG] = sorted(mgi2ensembl[mgiId])
                            self.incPPGCount(attrs[PPG])
                # Top level feaure with no usable IDs in the attributes
                if 'Name' in attrs and PPG not in attrs:
                    # Lookup the symbol.
                    mgiId = symbol2mgi.get(attrs['Name'], None)
                    ensemblId = mgi2ensembl.get(mgiId, None)
                    if ensemblId:
                        attrs[PPG] = sorted(ensemblId)
                        self.incPPGCount(attrs[PPG])
            # If feature has a PPG, try to find old MGP ids and attach them as a new "mgp" attribute
            if PPG in attrs:
                mgps = set()
                for ppg in attrs[PPG]:
                    for mgp in ensembl2mgp.get(STRAIN,{}).get(ppg, []):
                        mgps.add(mgp)
                if len(mgps) > 0:
                    attrs[MGP] = sorted(mgps)
            #
            results.append(gff3lite.formatLine(f))
        #
        return results


    def doPatching (self) :
        ensembl2mgp = self.getMGPids(PATCH_ARCHIVED_MGP_IDS)
        symbol2mgi, entrez2mgi, mgi2ensembl = self.getMgiGeneModelIds()
        lines = self.process(entrez2mgi, mgi2ensembl, symbol2mgi, ensembl2mgp)
        return lines

    def applyLimit (self, lines) :
        for (ppg, n) in self.PPG2count.items() :
            if n > self.LIMIT:
                self.log('%s %d' % (ppg, n))

        for (i,line) in enumerate(lines):
            if PPG in line:
                f = gff3lite.parseLine(line)
                attrs = f[8]
                ppg = attrs[PPG]
                if self.PPG2count.get(ppg,0) > self.LIMIT:
                    attrs.pop(PPG)
                    if MGP in attrs:
                        attrs.pop(MGP)
                    lines[i] = gff3lite.formatLine(f)
        return lines

    def writeOutput (self, lines) :
        for l in lines:
            self.ofd.write(l)

    def main (self) :
        args = self.getArgs()
        results = self.doPatching()
        results = self.applyLimit(results)
        self.writeOutput(results)

if __name__ == '__main__':
    Patcher().main()

