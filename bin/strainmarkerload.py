#!/usr/local/bin/python
#
#  strainmarkerload.py
###########################################################################
#
#  Purpose:
#
#      Validate/QC input and create Strain Marker bcp file
#
#  Usage:
#
#      strainmarkerload.py
#
#  Inputs:
#
#	1. 18 strain specific GFF3 files
#	2. For B6 MGI database
#	3. Configuration - see strainmarkerload.config
#
#  Outputs:
#
#       1. MRK_StrainMarker.bcp
#	2. ACC_Accession.bcp
#	3. ACC_AccessionReference.bcp
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#      2:  bcp fails

#  Assumes:
#
#	
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Validate the arguments to the script.
#      2) Perform initialization steps.
#      3) Open the input/output files.
#      4) Run the QC  checks
#      5) Run the load if QC checks pass
#      6) Close the input/output files.
#      7) Delete existing strain marker records
#      8) BCP in new strain marker records
#
#  Notes:  None
#
###########################################################################
#
#  Modification History:
#
#  Date        SE   Change Description
#  ----------  ---  -------------------------------------------------------
#
#  04/09/2018  sc  Initial development
#
###########################################################################

import sys
import os
import string
import Set
import re

import db
import mgi_utils
import loadlib
import accessionlib

#
#  CONSTANTS
#
TAB = '\t'
CRT = '\n'
DATE = mgi_utils.date("%m/%d/%Y")
USAGE='strainmarkerload.py'

# sequence description templates

# Chr<chr>:<start>-<end>, <strand> strand. MGI derived this sequence for the C57BL/6J strain version of Gene: <marker symbol>, Gene type: <feature type>, from outermost boundary coordinates of combined annotations to mouse reference assembly GRCm38 provided by: <comma delimited provider:ID>
 
b6DescriptTemplate = "Chr%s:%s-%s, %s strand. MGI derived this sequence for the C57BL/6J strain version of Gene: %s, Gene type: %s, from outermost boundary coordinates of combined annotations to mouse reference assembly GRCm38 provided by:  %s"

# Chr<chr>:<start>-<end>, <strand> strand. MGI derived this sequence for the C57BL/6J strain version of Gene: <marker symbol>, Gene type: <feature type>, from outermost boundary coordinates of combined BLAT alignments to the mouse reference assembly GRCm38 for sequences: <comma delimited list of genbank ids>.

b6BlatDescriptTemplate =   "Chr%s:%s-%s, %s strand. MGI derived this sequence for the C57BL/6J strain version of Gene: %s, Gene type: %s, from outermost boundary coordinates of combined BLAT alignments to the mouse reference assembly GRCm38 for sequences: %s."

mgpDescriptTemplate = '"Chr%s:%s-%s, %s strand. xxx %s strain version of Gene: %s, Gene type: %s."'

loaddate = loadlib.loaddate
#
#  GLOBALS
#

DEBUG = 0

# true if we want to run QC and not load relationships
QC_ONLY = os.environ['QC_ONLY']

# if true only load B6 strain markers
loadOnlyB6  = os.environ['B6_ONLY']

# min number of records expected
MIN_RECORDS = int(os.environ['MIN_RECORDS'])

# accession ID logicalDB keys
mgpLDBKey = 209 	# MGP
smLDBKey = 212		# strain marker
ensLDBKey = 60		# ENSEMBL
ncbiLDBKey = 59		# NCBI
mirLDBKey = 83		# miRBase
gbLDBKey = 9		# GenBank

# MRK_StrainMarker MGI Type
mgiTypeKey = 44

# if true do not delete/reload
isFatal = 0

# input
infileDir = os.environ['INPUTDIR']
inputFileString =  os.environ['INPUT_MGP_DIR_LIST']
b6InputFile = os.environ['INPUT_MGI_GFF_FILE']
fpB6Inputfile = ''

# output files
curLog = os.getenv('LOG_CUR')

# output bcp files
outputDir = os.environ['OUTPUTDIR']

smBcpFile =   os.environ['SM_BCP_FILE']
strainMarkerFile = '%s/%s' % (outputDir, smBcpFile)
fpStrainMarkerFile = ''

accBcpFile = os.environ['ACC_BCP_FILE']
accFile = '%s/%s' % (outputDir, accBcpFile)
fpAccFile = ''

accRefBcpFile = os.environ['ACC_REF_BCP_FILE']
accRefFile = '%s/%s' % (outputDir, accRefBcpFile)
fpAccRefFile = ''

# output for downstream loads
gmMgpFile = os.environ['GM_MGP_INPUT_FILE']
fpGmMgpFile = ''
biotypeMgpFile = os.environ['GM_MGP_BIOTYPE_FILE']
fpBiotypeMgpFile = ''

gmB6File = os.environ['GM_B6_INPUT_FILE']
fpGmB6File = ''
biotypeB6File = os.environ['GM_B6_BIOTYPE_FILE']
fpBiotypeB6File = ''

# QC reporting data structures
qcDict = {} 
messageMap = {}

# Lookups
strainTranslationLookup = {} # {badName: _Strain_key, ...}
markerLookup = {}            # {MGI ID: Marker ...}
chrLookup = {}		     # {Mouse, Laboratory chr: chrKey, ...}
biotypeLookup = {}           # {raw biotype:feature type, ...}
mcvTermLookup = []	     # lookup feature type (B6)

# the MGP Strain Marker reference key for J:259852
mgpRefsKey = 282407

# the MGI B6 Strain Marker reference key for J:260092
b6RefsKey =  282660

# the MGI C57BL/6J strain key
b6StrainKey = 38048

# Data structure for B6 lines parsed from input by MGI ID
b6ToLoadDict = {}

# strainmarkerload user key
userKey = 1600

# database primary keys, will be set to the next available from the db
nextSMKey = None	# MRK_StrainMarker._StrainMarker_key
nextAccKey = None	# ACC_Accession._Accession_key

# for bcp
bcpin = '%s/bin/bcpin.csh' % os.environ['PG_DBUTILS']
server = os.environ['MGD_DBSERVER']
database = os.environ['MGD_DBNAME']
strainmarker_table = 'MRK_StrainMarker'
acc_table = 'ACC_Accession'
accref_table = 'ACC_AccessionReference'

# count of total strain markers that will be loaded
totalLoadedCt = 0
b6LoadedCt = 0 # b6 only count
mgpFileCt = 0 # total in mgp file
mgpLoadCt = 0 # gotal mgp loaded
mgpSkipCt = 0 # total mgp skipped
mgpNoMarkerCt = 0 # total mgp w/no marker
ctByStrain = {} # {strain: ct, ...}

class Marker:
    # Is: data object for a marker
    # Has: a set of marker attributes
    # Does: provides direct access to its attributes
    #
    def __init__ (self):
        # Purpose: constructor
        self.markerID = None
	self.symbol = None
        self.markerPreferred = None
        self.markerStatus = None

    def toString(self):
	return '%s, %s, %s, %s' % (self.markerID, self.symbol, self.markerPreferred, self.markerStatus)
# end class Marker ---------------------------------

class StrainMarker:
    # Is: data object for a strain marker
    # Has: a set of strain marker attributes
    # Does: provides direct access to its attributes
    #
    def __init__ (self):
        # Purpose: constructor
	self.markerID = ''
	self.markerKey = ''
	self.strainKey = ''
	self.mgpIDs = set()
	self.chr = ''
	self.start = ''
	self.end = ''
	self.strand = ''
	self.description = ''
	self.biotype = ''
        
    def toString(self):
        return '%s, %s, %s, %s %s %s %s %s %s' % (self.markerID, self.markerKey, self.strainKey, self.mgpIDs, self.chr, self.start, self.end, self.strand, self.description)
# end class StrainMarker ----------------------------

def checkArgs ():
    # Purpose: Validate the arguments to the script.
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: exits if unexpected args found on the command line
    # Throws: Nothing

    if len(sys.argv) != 0:
        print USAGE
        return 1

    return 0

# end checkArgs() -------------------------------

def init():
    # Purpose: check args, create lookups, open files, create db connection, 
    #   gets max keys from the db
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened,
    #  creates files in the file system, creates connection to a database

    global nextSMKey, nextAccKey, strainTranslationLookup, markerLookup, chrLookup
    global biotypeLookup, mcvTermLookup, messageMap
   
    checkArgs()

    #
    # Open input and output files
    #
    openFiles()

    #
    # create database connection
    #
    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    #
    # get next MRK_StrainMarker key
    #
    results = db.sql('''select max(_StrainMarker_key) + 1 as nextSMKey
	    from MRK_StrainMarker''', 'auto')
    if results[0]['nextSMKey'] is None:
	nextSMKey = 1000
    else:
	nextSMKey = results[0]['nextSMKey']

    #
    # get next ACC_Accession key
    #
    results = db.sql('''select max(_Accession_key) + 1 as nextAccKey
	    from ACC_Accession''', 'auto')
    nextAccKey = results[0]['nextAccKey']

    # load qcDict with keys
    qcDict['chr_m'] = []     # chr is missing, report/skip 
    qcDict['chr_u'] = []     # chromosome unresolved, report/skip
    qcDict['start'] = []     # startCoordinate is missing, report/skip
    qcDict['end'] = []       # endCoordinate is missing, report/skip
    qcDict['start/end'] = [] # start > end, report/skip
    qcDict['strand'] = []    # strand is missing, report/skip
    qcDict['biotype_m'] = [] # biotype is missing, report/skip
    qcDict['biotype_u'] = {} # biotype unresolved {biotype:count}, report/skip
    qcDict['strain_u'] = []  # strain unresolved, fatal
    qcDict['mgp'] = []       # mgp is missing, report/skip
    qcDict['mgi_u'] = []     # marker ID unresolved, report create strain marker with null marker
    qcDict['mgi_s'] = []     # marker ID secondary, report, load
    qcDict['mgi_no'] = []    # marker ID not official, report, skip
    qcDict['mgi_mgp'] = []   # list of {mgiID: ([set of mpIDs]), ...}, one for each strain file
			     # used to a) determine multiple MGP IDs (within a given strain) per marker, report/skip
			     # b) write out to bcp files

    messageMap['chr_m'] = 'Chromosome missing from input, record(s) skipped'
    messageMap['chr_u'] = 'Chromosome from input unresolved, record(s) skipped'
    messageMap['start'] = 'Start Coordinate missing from input, record(s) skipped' 
    messageMap['end'] = 'End Coordinate missing from input, record(s) skipped'
    messageMap['start/end'] = 'Start Coordinate > End Coordinate, record(s) skipped'
    messageMap['strand'] = 'Strand missing from input, record(s) skipped'
    messageMap['biotype_m'] = 'Biotype missing from input, record(s) skipped'
    messageMap['biotype_u'] = 'Biotype from input unresolved, record(s) skipped'
    messageMap['strain_u'] = 'Strain from input unresolved, load fails'
    messageMap['mgp'] = 'MGP ID missing from input, record(s) skipped'
    messageMap['mgi_u'] = 'Marker from input unresolved, strain marker created with null marker'
    messageMap['mgi_s'] = 'Marker from input is secondary ID, strain marker created'
    messageMap['mgi_no'] = 'Marker from input is not official, records(s) skipped'
    messageMap['mgi_mgp'] = 'Markers from input with > 1 Strain specific MGP ID, record(s) skipped'  # when reporting multis. 

    #
    # create lookups
    #

    # load lookup of strain translations
    results = db.sql('''select badName, _Object_key as strainKey
	from MGI_Translation
	where _TranslationType_key = 1021''', 'auto')
    for r in results:
	strainTranslationLookup[r['badName']] = r['strainKey']

    # load lookup of all marker MGI IDs
    results = db.sql('''select m._Marker_key, m.symbol, s.status as markerStatus,
	    a.accid as mgiID, a.preferred
	from ACC_Accession a, MRK_Marker m, MRK_Status s
	where a. _MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a._Object_key = m._Marker_key
	and m._Organism_key = 1
	and m._Marker_Status_key = s._Marker_Status_key''', 'auto')
    for r in results:
	m = Marker()
	m.markerKey = r['_Marker_key']
	m.markerID = r['mgiID']
	m.symbol = r['symbol']
	m.markerStatus = r['markerStatus']
	m.markerPreferred = r['preferred']

	markerLookup[m.markerID] = m
    
    # load lookup of mouse, laboratory chromosomes
    results = db.sql('''select chromosome, _Chromosome_key
	from MRK_Chromosome
	where _Organism_key = 1''', 'auto')
    for r in results:
	chrLookup[r['chromosome']] = r['_Chromosome_key']

    # load lookup of raw MGP biotype to feature type
    results = db.sql('''select t1.term as rawBiotype, t2.term as mcvTerm
        from VOC_Term t1, VOC_Term t2, 
	    MRK_BiotypeMapping m
        where t1._Vocab_key = 136 --Biotype MGP
        and t1._Term_key = m._BiotypeTerm_key
	and m._mcvTerm_key = t2._Term_key''', 'auto')
    for r in results:
	biotypeLookup[r['rawBiotype'].lower()] = r['mcvTerm']
    #print 'biotypeLookup: %s' % biotypeLookup

    # load lookup of feature type vocabulary
    results = db.sql('''select term from VOC_Term
	where _Vocab_key = 79''')
    for r in results:
	mcvTermLookup.append(r['term'].lower())

    return 0

# end init() -------------------------------

def openFiles ():
    # Purpose: Open input/output files.
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened, 
    #  creates files in the file system

    global fpStrainMarkerFile, fpAccFile, fpAccRefFile, fpLogCur
    global fpGmMgpFile, fpBiotypeMgpFile
    global fpB6InputFile, fpGmB6File, fpBiotypeB6File

    try:
        fpStrainMarkerFile  = open(strainMarkerFile, 'w')
    except:
        print 'ERROR: Cannot open Strain Marker BCP file: %s' % strainMarkerFile
        sys.exit(1)
    try:
        fpAccFile  = open(accFile, 'w')
    except:
        print 'ERROR: Cannot open ACC_Accession bcp  file: %s' % accFile
        sys.exit(1)
    try:
        fpAccRefFile  = open(accRefFile, 'w')
    except:
        print 'ERROR: Cannot open ACC_AccessionReference bcp file: %s' % accRefFile
        sys.exit(1)

    try:
	fpLogCur = open(curLog, 'a')
    except:
	print 'ERROR: Cannot open Curator Log file: %s' % curLog
        sys.exit(1)
    try:
        fpGmMgpFile = open(gmMgpFile, 'w')
    except:
        print 'ERROR: Cannot open MGP Gene Model file: %s' % gmMgpFile
        sys.exit(1)
    try:
        fpBiotypeMgpFile = open(biotypeMgpFile, 'w')
    except:
        print 'ERROR: Cannot open MGP Gene Model Biotype file: %s' % biotypeMgpFile
        sys.exit(1)
    try:
        fpB6InputFile  = open(b6InputFile, 'r')
    except:
        print 'ERROR: Cannot open MGI.gff3 B6 file: %s' % b6InputFile
        sys.exit(1)
    try:
        fpGmB6File = open(gmB6File, 'w')
    except:
        print 'ERROR: Cannot open MGI B6 Gene Model file: %s' % gmB6File
        sys.exit(1)
    try:
        fpBiotypeB6File = open(biotypeB6File, 'w')
    except:
        print 'ERROR: Cannot open MGI B6 Gene Model Biotype file: %s' % biotypeB6File
        sys.exit(1)

    return 0

# end openFiles() -------------------------------

def closeFiles ():
    # Purpose: Close all file descriptors
    # Returns: 1 if error, else 0
    # Assumes: all file descriptors were initialized
    # Effects: Nothing
    # Throws: Nothing
    
    try:
	fpStrainMarkerFile.close()
	fpAccFile.close()
        fpAccRefFile.close()
	fpGmMgpFile.close()
	fpBiotypeMgpFile.close()
        fpB6InputFile.close()
        fpGmB6File.close()
        fpBiotypeB6File.close()
    except:
	return 1
    return 0

# end closeFiles() -------------------------------

def parseMGPFiles( ): 
    # Purpose: parses input files looping, opening and closing
    # Returns: 1 if error, else 0
    # Assumes: output file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing
   
    global qcDict, mgpFileCt, mgpLoadCt, mgpSkipCt, ctByStrain, mgpNoMarkerCt

    inputFileList = string.split(inputFileString, ' ')
    
    for file in inputFileList:
 	recordCt = 0  # current number of records in this file
        file = file.strip() # in case extra spaces btwn filenames in config
	inputFile = '%s/%s.txt' % (infileDir, file)
   	print 'inputFile: %s' % inputFile
        fpIn = open(inputFile, 'r')

	# strain is the first line in the file
        strain = string.strip(fpIn.readline())

	# resolve strain with translation lookup
	if strain not in strainTranslationLookup:
	    qcDict['strain_u'].append(strain)
	    continue 
	# get the strain key
        strainKey = strainTranslationLookup[strain] 
	print 'strain: %s strainKey: %s' % (strain, strainKey)
	# build this as we parse each file - adding mgpIDs to strainMarkerObject if mgiID
	# found > 1 in file
	strainMarkerDict = {}   # {mgiID:strainMarkerObject, ...}

	# after file parsed copy the strainMarkerObjects from strainMarkerDict to this mapping by strain
	strainMarkerInput = {} 		# {strain: list of strainMarkerObjects, ...}
	strainMarkerInput[strain] = []	# initialize 

	# iterate thru lines in this strain file
	for line in fpIn.readlines():
	    recordCt +=1
	    mgpFileCt += 1
	    isSkip = 0  # if errors, set to true (1) to skip this record
	    #print 'line: %s' % line
	    tokens = line.split('\t')
	    chr = tokens[0]
	    start = tokens[3]
	    end = tokens[4]
	    strand = tokens[6]
	    # example col9:
	    # ID=gene:MGP_CAROLIEiJ_G0013919;Name=Xkr4;biotype=protein_coding;description=X-linked Kx blood group related 4 [Source:MGI Symbol%3BAcc:MGI:3528744];gene_id=MGP_CAROLIEiJ_G0013919;logic_name=mgp_import;projection_parent_gene=ENSMUSG00000051951.5;version=1
	    col9 = tokens[8]
            #print 'col9: %s' % col9
	    tokens2 = col9.split(';')
	    #print 'tokens2: %s' % tokens2
	    mgpID = ''
	    mgiID = ''
	    symbol = ''
	    markerKey = ''
	    biotype = ''
	    for t in tokens2: # col9 tokens
		#print 'token in col9: %s' % t
		# 'ID=gene:MGP_CAROLIEiJ_G0013919'
		if t.find('ID=gene:') != -1:
		    mgpID = t.split(':')[1]
		# 'description=X-linked Kx blood group related 4 [Source:MGI Symbol%3BAcc:MGI:3528744]'
		elif t.find('Acc:MGI:') != -1:
		    mgiRE = re.compile('(MGI:[0-9]*)')
		    match = mgiRE.search(t)
		    mgiID = match.group(1)
		elif t.find('biotype=') != -1:
		    biotype = t.split('=')[1]
		    #print 'biotype=%s' % biotype
	    #print 'mgpID: %s mgiID: %s biotype: %s' % (mgpID, mgiID, biotype)
	    #print 'start: "%s" end: "%s"' % (start, end)
	    if chr == '':
		qcDict['chr_m'].append(line)
		isSkip = 1
	    if chr != '' and chr not in chrLookup:
		qcDict['chr_u'].append(line)
	    if start == '':
		qcDict['start'].append(line)
		isSkip = 1
  	    if end == '':
		qcDict['end'].append(line)
		isSkip = 1
	    if start != '' and end != '' and (int(start) > int(end)):
		qcDict['start/end'].append(line)
		isSkip = 1
	    if strand == '':
		qcDict['strand'].append(line)
		isSkip = 1
	    if biotype == '':
		qcDict['biotype_m'].append(line)
		isSkip = 1
	    if mgpID == '':
		qcDict['mgp'].append(line)
		isSkip = 1

            # check that biotype is in the database
	    # don't report missing biotype as unresolved
	    biotypeLower = biotype.lower().strip()
            if biotypeLower != '' and biotypeLower not in biotypeLookup:
		if  biotype not in qcDict['biotype_u']:     
		    qcDict['biotype_u'][biotype] = 1
		else:
		    qcDict['biotype_u'][biotype] += 1
		isSkip = 1
	    else: # resolve the rawbiotype to feature type term
		biotype = biotypeLookup[biotypeLower]
	    # Resolve MGI ID 
	    if mgiID != '':
		if mgiID not in markerLookup:
		    qcDict['mgi_u'].append('%s : %s' % (mgiID, line))
		    #print '%s NOT IN MGI in strain file: %s' % (mgiID, strain)
		else:
		    marker = markerLookup[mgiID]
		    markerKey = marker.markerKey 
		    symbol = marker.symbol

		    # not a preferred ID (secondary ID) - report and load
		    if marker.markerPreferred != 1:
			qcDict['mgi_s'].append('%s : %s' % \
			    (mgiID, line))
			#print '%s SECONDARY ID in strain file: %s' % (mgiID, strain)

		    # not official marker status - report and skip
		    if marker.markerStatus != 'official':
			markerKey = '' # Don't load unofficial markers
			qcDict['mgi_no'].append('%s : %s' % (mgiID, line))
			#print '%s MARKER NOT OFFICIAL in strain file: %s' % (mgiID, strain)
			isSkip = 1
	    if isSkip:
		mgpSkipCt += 1
	    if not isSkip:
		mgpLoadCt += 1
		if markerKey == '': # count them
		    mgpNoMarkerCt += 1
		    # create a temporary ID ofr markerless strain marker object - each must have its own uniq
		    # set of coordinate attributes
		    mgiID = 'TEMP:%s' % mgpNoMarkerCt
		description = mgpDescriptTemplate % (chr, start, end, strain, strain, symbol, biotype)
		strainMarkerObject = StrainMarker() # default to a new one
	        if mgiID not in strainMarkerDict:
		    strainMarkerObject.markerID = mgiID
		    strainMarkerObject.markerKey = markerKey
		    strainMarkerObject.strainKey = strainKey
		    strainMarkerObject.mgpIDs.add(mgpID)
		    strainMarkerObject.chr = chr
		    strainMarkerObject.start = start
		    strainMarkerObject.end = end
		    strainMarkerObject.strand = strand
		    strainMarkerObject.biotype = biotype
		    strainMarkerObject.description = description
		else:   # later we will eliminate any strainmarkers w/multi
			# MGP Ids
		    strainMarkerObject = strainMarkerDict[mgiID]
                    strainMarkerObject.mgpIDs.add(mgpID)
		#print strainMarkerObject.toString()
		strainMarkerDict[mgiID] =  strainMarkerObject
	ctByStrain[strain] =  recordCt
	# ------------ end for line in file --------------------

	# copy strainMarker objects from strainMarkerDict(by MGI ID) to strainMarkerInput (By strain)
	for mgiID in strainMarkerDict:
	    strainMarkerInput[strain].append(strainMarkerDict[mgiID])

	# add this strain to the qcDict
        qcDict['mgi_mgp'].append(strainMarkerInput)
    return 0

# end parseMGPFiles() -------------------------------------

def writeMGPOutput():
    # Purpose: writes to Accession, AccessionReference & StrainMarker 
    #	BCP file and Gene Model and GM Assoc files if there are no errors
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global nextSMKey, nextAccKey, mgpSkipCt, totalLoadedCt

    strainMarkerInputList = qcDict['mgi_mgp']
    for strainMarkerInputDict in strainMarkerInputList:
        for strain in strainMarkerInputDict:
            strainMarkerObjectsList = strainMarkerInputDict[strain]
            for strainMarkerObject in strainMarkerObjectsList:
                mgpList = list(strainMarkerObject.mgpIDs) # get the list of mgp IDs
                if len(mgpList) > 1:
		    mgpSkipCt += len(mgpList)
		    continue  # these were reported
	        # otherwise write out to bcp file
	  	mgpID = mgpList[0]
		mgiID = strainMarkerObject.markerID 
		markerKey = strainMarkerObject.markerKey
		if mgiID.find('TEMP:') == 0: # temp ID for no marker strain marker
		    markerKey = ''
		    
		strainKey = strainMarkerObject.strainKey
		chr = strainMarkerObject.chr
		start = strainMarkerObject.start
		end = strainMarkerObject.end
		strand = strainMarkerObject.strand
		description = strainMarkerObject.description
		biotype = strainMarkerObject.biotype

		totalLoadedCt += 1
		fpStrainMarkerFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextSMKey, TAB, strainKey, TAB, markerKey, TAB, mgpRefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

		prefixPart, numericPart = accessionlib.split_accnum(mgpID)

		fpAccFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s0%s1%s%s%s%s%s%s%s%s%s' \
		% (nextAccKey, TAB, mgpID, TAB, prefixPart, TAB, numericPart, TAB, mgpLDBKey, TAB, nextSMKey, TAB, mgiTypeKey, TAB, TAB, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))
		fpAccRefFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' \
		% (nextAccKey, TAB, mgpRefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))
		fpGmMgpFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' % (mgpID, TAB, chr, TAB, start, TAB, end, TAB, strand, TAB, description, CRT))
		fpBiotypeMgpFile.write('%s%s%s%s' % (mgpID, TAB, biotype, CRT))
		nextSMKey += 1
		nextAccKey += 1

    return 0

# end writeMGPOutput() ---------------------------------------------------

def parseB6File( ):
    # Purpose: parses  MGI.gff3 file
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing
        # example with dbXref in col9
        # 1       MGI     gene    4807560 4848410 .       +       .       ID=MGI:1344588;Name=Lypla1;mgi_type=protein coding gene;so_term_name=protein_coding_gene;gene_id=MGI:1344588;curie=MGI:1344588;strain_gene_id=MGI_C57BL6J_1344588;Dbxref=ENSEMBL:ENSMUSG00000025903,NCBI_Gene:18777,NCBI_Gene:105243856;description=lysophospholipase 1

        # example  with no dbXref in col9
        # 1       MGI     gene    5156597 5158506 .       +       .       ID=MGI:2443922;Name=B230334L07Rik;curie=MGI:2443922;so_term_name=gene;mcv_type=unclassified gene;description=RIKEN cDNA B230334L07 gene

        # example BlatAlignment for above example
        # 1       BlatAlignment   match   5156597 5158506 .       +       .       ID=MGI:2443922.m2;Name=AK046028.1;Parent=MGI:2443922;qSize=1910;pctLen=100.0;matches=1909;qName=AK046028.1;mgi_id=MGI:2443922;matchLen=1910;pctIdentity=99.9476439791;qEnd=1910;qStart=0

    global b6ToLoadDict 

    # iterate thru lines in the B6 file
    for line in fpB6InputFile.readlines():
	if line.find('#') == 0:
	    continue
	
	tokens = line.split('\t')
	if tokens[1] == 'BlatAlignment':
	    feature = tokens[1]
	else:
	    feature = tokens[2]
	# these are the only features we are interested in
	if feature not in ['gene', 'pseudogene', 'BlatAlignment']:
	    #print 'badLine: %s' % line
	    continue

	#print 'goodLine: %s' % line
	col9 = tokens[8]
	#print 'col9: %s' % col9
	tokens2 = col9.split(';')

	mgiID = ''
	for t in tokens2: # col9 tokens
	    #print 'token in col9: %s' % t
	    if t.find('curie=MGI:') != -1: # top level feature
		mgiID = t.split('=')[1]
		#print 'found curie=MGI: %s' % mgiID
	    elif t.find('mgi_id=') != -1:  # blat hit
		mgiID = t.split('=')[1]
		#print 'found mgi_id: %s' % mgiID
	if mgiID != '':
	    if mgiID not in b6ToLoadDict:
		b6ToLoadDict[mgiID] = []
	    b6ToLoadDict[mgiID].append(line)
    return 0

# end parseB6File() ---------------------------------------------------

def parseB6Feature(line, type):
    #print 'In parseB6Feature type: %s line: %s' % (type, line)
    chr = ''
    start = ''
    end = ''
    strand = ''
    tokens = line.split('\t')
    chr = tokens[0]
    start = tokens[3]
    end = tokens[4]
    strand = tokens[6]
    col9 = tokens[8]
    #print 'col9: %s' % col9
    tokens2 = col9.split(';')
    #print 'tokens2: %s' % tokens2
    # This is the strain/marker ID in col9 e.g. strain_gene_id=MGI_C57BL6J_1344588
    smID = ''
    mgiID = ''
    markerKey = ''
    symbol = ''
    biotype = ''
    gmIdString = ''
    qName = ''
    description = ''
    for t in tokens2: # col9 tokens
	#print 'token in col9: %s' % t
	# 'ID=MGI:1344588'
	if t.find('curie=MGI:') != -1:
	    mgiID = t.split('=')[1]
	    #print 'found curie=MGI: %s' % mgiID
	elif t.find('ID=') != -1:
	    smID = t.split('=')[1]
	    #print 'found strain_gene_id: %s' % smID
	elif t.find('mgi_type=') != -1:
	    biotype = t.split('=')[1]
	    #print 'found mgi_type=%s' % biotype
	# Dbxref=miRBase:MI0005004,ENSEMBL:ENSMUSG00000076010,NCBI_Gene:751557
	elif t.find('Dbxref=') != -1:
	    gmIdString = t.split('=')[1]
	elif t.find('qName=') != -1:
	    qName = t.split('=')[1]
	    qName = qName.split('.')[0]
    #print 'smID: %s mgiID: %s biotype: %s gmIdString: %s' % (smID, mgiID, biotype, gmIdString)
    #print 'chr: "%s" start: "%s" end: "%s" strand: "%s"' % (chr, start, end, strand)
    # IMPLEMENTATION NOTE: We expect no errors given that this data is from 
    # Joel's gff3 file. I added print statements to confirm there were no errors
    # there were not, we don't expect any, so won't take the time to log, just 
    # print
    if chr == '':
	#qcDict['chr_m'].append(line)
	if type == 'f':
	    print 'ERROR: no chr'
    if  chr not in chrLookup:
	#qcDict['chr_u'].append(line)
	if type == 'f':
	    print 'ERROR: chr not resolved'
    if start == '':
	#qcDict['start'].append(line)
	if type == 'f':
	    print 'ERROR: no start coord'
    if end == '':
	#qcDict['end'].append(line)
	if type == 'f':
	    print 'ERROR: no end coord'
    if start != '' and end != '' and (int(start) > int(end)):
	#qcDict['start/end'].append(line)
	if type == 'f':
	    print 'ERROR: start > end'
    if strand == '':
	#qcDict['strand'].append(line)
	if type == 'f':
	    print 'ERROR: missing strand'
    if biotype == '':
	#qcDict['biotype_m'].append(line)
	if type == 'f':
	    print 'ERROR: missing biotype'
    if smID == '':
	#qcDict['mgp'].append(line)
	if type == 'f':
	    print 'ERROR: missing strain/marker id'
    if mgiID == '':
	if type == 'f' or type == 'bf':
	    print 'ERROR: missing mgi ID'
    else:
	symbol = markerLookup[mgiID].symbol
    if qName == '' and type == 'b':
	print 'ERROR: missing qName: %s' % qName
    # check that feature is in the database
    if biotype.lower().strip() not in mcvTermLookup:
	if type == 'f':
	    print 'ERROR: biotype not in MGI'
    # Nothing is reported from the adhoc QC checks above
    
    # calculate the description
    if type == 'bf': # blat feature
	description = b6BlatDescriptTemplate % (chr, start, end, strand, symbol, biotype, gmIdString ) # need to get prefixed gmIds in gmIdString
    elif type == 'f': # feature
	description = b6DescriptTemplate % (chr, start, end, strand, symbol, biotype, qName)
    else:
	pass # type == 'b' gets default empty string 
	
    return [chr, start, end, strand, smID, mgiID, biotype, gmIdString, qName, description]

# end parseB6Feature(line) ---------------------------------------------------

def writeB6Output():
    # Purpose: parses the output line dictionary
    # writes to Accession, AccessionReference & StrainMarker BCP file 
    # if there are no errors
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global  nextSMKey, nextAccKey, totalLoadedCt, b6LoadedCt

    #print 'in writeB6Output len of b6ToLoadDict: %s' % len(b6ToLoadDict)

    description = ''
    for mgiID in b6ToLoadDict:
	lineList = b6ToLoadDict[mgiID]
	qNameSet = set()
	#print mgiID
	#print lineList
	# Resolve MGI ID
	if mgiID not in markerLookup:
	    #qcDict['mgi_u'].append('%s : %s' % (mgiID, line))
	    print '%s in MGI GFF File, but NOT IN MGI' % (mgiID)
	    continue
	else:
	    marker = markerLookup[mgiID]
	    markerKey = marker.markerKey
	    symbol = marker.symbol
	if len(lineList) == 1:  # This is non-BlatAlignment gene/pseudogene
	    #print 'NOT blat alignment'
	    line = lineList[0]
	    chr, start, end, strand, smID, mgiID, biotype, gmIdString, qName, description  = parseB6Feature(line, 'f')
	    fpStrainMarkerFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextSMKey, TAB, b6StrainKey, TAB, markerKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    prefixPart, numericPart = accessionlib.split_accnum(smID)

	    fpAccFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s0%s1%s%s%s%s%s%s%s%s%s' \
		% (nextAccKey, TAB, smID, TAB, prefixPart, TAB, numericPart, TAB, smLDBKey, TAB, nextSMKey, TAB, mgiTypeKey, TAB, TAB, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    fpAccRefFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' \
		% (nextAccKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))
	    
	    fpGmB6File.write('%s%s%s%s%s%s%s%s%s%s%s%s' % (smID, TAB, chr, TAB, start, TAB, end, TAB, strand, TAB, description, CRT))
            fpBiotypeB6File.write('%s%s%s%s' % (smID, TAB, biotype, CRT))
  
	    nextAccKey += 1
	
	    # get gmIDs from the input file for the sequence description
	    # gmIDs example:
	    # Dbxref=miRBase:MI0005004,ENSEMBL:ENSMUSG00000076010,NCBI_Gene:751557
	    if gmIdString == '': 
		continue
	    gmIdList = gmIdString.split(',')
	    #
	    # 6/12 GF-184, remove all associated IDs from strain gene
	    #

	    #for id in gmIdList:
	    #	provider, ID = id.split(':')
	    #	if provider == 'miRBase':
	    #	    ldbKey = mirLDBKey
	    #	elif provider == 'ENSEMBL':
	    #	    ldbKey = ensLDBKey
	    #	elif provider == 'NCBI_Gene':
	    #	    ldbKey = ncbiLDBKey
		
		#prefixPart, numericPart = accessionlib.split_accnum(ID)
		# 6/12 GF-184, remove all associated IDs from strain ge
		# fpAccFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s0%s1%s%s%s%s%s%s%s%s%s' \
		#    % (nextAccKey, TAB, ID, TAB, prefixPart, TAB, numericPart, TAB, ldbKey, TAB, nextSMKey, TAB, mgiTypeKey, TAB, TAB, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))
		#fpAccRefFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' \
		#    % (nextAccKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

		#nextAccKey += 1
	    nextSMKey += 1
	else: # This is BlatAlignment set
	    #print 'BlatAlignment:'
	    #print mgiID
	    # The first line is feature line, the following are 
	    # BlatAlignments
	    featureLine = lineList[0]
	    #print 'featureLine: %s' % featureLine

	    chr, start, end, strand, smID, mgiID, biotype, gmIdString, qName, description  = parseB6Feature(featureLine, 'bf')	    

	    # remainder of the list are blat hits - save them to a set
	    lineList = lineList[1:] # remove the blat feature
	    for line in lineList:
		#print 'blatLine: %s' % line
		j1, j2, j3, j4, j5, j6, j7, j8, qName, j9 = parseB6Feature(line, 'b')
		qNameSet.add(qName.strip())
	    #print 'qNames: %s ' % qNameSet

	    #
	    # Create the strain marker and its accession ID
	    #
	    fpStrainMarkerFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextSMKey, TAB, b6StrainKey, TAB, markerKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    prefixPart, numericPart = accessionlib.split_accnum(smID)
	    fpAccFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s0%s1%s%s%s%s%s%s%s%s%s' \
		% (nextAccKey, TAB, smID, TAB, prefixPart, TAB, numericPart, TAB, smLDBKey, TAB, nextSMKey, TAB, mgiTypeKey, TAB, TAB, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    fpAccRefFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' \
		% (nextAccKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    nextAccKey += 1

	    fpGmB6File.write('%s%s%s%s%s%s%s%s%s%s%s%s' % (smID, TAB, chr, TAB, start, TAB, end, TAB, strand, TAB, description, CRT))
            fpBiotypeB6File.write('%s%s%s%s' % (smID, TAB, biotype, CRT))
	    
	    # 6/12 GF-184, remove all associated IDs from strain gene	
	    # for blat hits in the input file associate the GenBank IDs that was
	    # blatted for coordinates
	    #ldbKey = gbLDBKey
	    #for ID in qNameSet:
	    #	prefixPart, numericPart = accessionlib.split_accnum(ID)
	       
	    #	fpAccFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s0%s1%s%s%s%s%s%s%s%s%s' \
	    #	    % (nextAccKey, TAB, ID, TAB, prefixPart, TAB, numericPart, TAB, ldbKey, TAB, nextSMKey, TAB, mgiTypeKey, TAB, TAB, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    #	fpAccRefFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' \
	    #	    % (nextAccKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	    #	nextAccKey += 1
	    nextSMKey += 1
	totalLoadedCt += 1
        b6LoadedCt += 1
    return 0

# end writeB6Output() ---------------------------------------------------

def writeCuratorLog():
    # Purpose: writes QC errors to the curator log
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global isFatal, mgpLoadCt, mgpNoMarkerCt
    #
    # special handling for MGI IDs mapped to multiple IDs for a given strain. 
    # Each dict in mgiMgpList represents the ids of a strain 
    #
    
    # unresolved strain is fatal
    if len(qcDict['strain_u']):
	isFatal = 1
	
    count = 0
    mgpSkipped = 0
    multiList = []
    strainMarkerInputList = qcDict['mgi_mgp']
    for strainMarkerInputDict in strainMarkerInputList:
	for strain in strainMarkerInputDict:
	    #print 'strain: %s' % strain
	    strainMarkerObjectsList = strainMarkerInputDict[strain]
	    for strainMarkerObject in strainMarkerObjectsList:
		mgpList = list(strainMarkerObject.mgpIDs) # get the list of mgp IDs
		#print 'strainMarkerObject.markerID: %s mgpList: %s' % (strainMarkerObject.markerID, mgpList)
		if strainMarkerObject.markerID != '' and len(mgpList) > 1:  # if > 1 ID, report
		    #print 'mgpList: %s' % mgpList
		    count += 1

		    # subtract from mgp no marker count if this strainmarker has
		    # no marker
		    if strainMarkerObject.markerKey == '': 
			mgpNoMarkerCt = mgpNoMarkerCt - 1
		    mgpSkipped += len(mgpList)
		    multiList.append('%s: %s' % (strainMarkerObject.markerID, string.join(mgpList, ', ')))

    # if we have multiple mpg IDs/marker within a strain, report
    if len(multiList):
	mgpLoadCt = mgpLoadCt - mgpSkipped # subtract multiples, won't load them
	msg = messageMap['mgi_mgp']
	fpLogCur.write('%s%s'% (msg,CRT))
	fpLogCur.write('-' * 80 + CRT)
	fpLogCur.write('%s%s' % (string.join(multiList, CRT), CRT))
	fpLogCur.write('Total: %s%s%s' % (count, CRT, CRT))

    #
    # special handling for unresolved biotypes
    # if we have unresolved biotypes print out the biotype with its count of MGP IDs
    #
    biotypeDict = qcDict['biotype_u']
    if len(biotypeDict):
	msg = messageMap['biotype_u']
        fpLogCur.write('%s%s'% (msg,CRT))
        fpLogCur.write('-' * 80 + CRT)
        for b in biotypeDict: 
    	    fpLogCur.write('%s: %s%s' % (b, biotypeDict[b], CRT))
	fpLogCur.write('Total: %s%s%s' % (len(biotypeDict), CRT, CRT))

    #
    #  process remaining QC in order specified by richard
    #
    order = ['strain_u', 'biotype_m', 'mgp', 'chr_u', 'chr_m', 'start', 'end', 'strand', 'start/end', 'mgi_no', 'mgi_u', 'mgi_s']
    for key in order:
        qcList = qcDict[key]
	if qcList == []:
	    continue
	msg = messageMap[key]
	fpLogCur.write('%s%s'% (msg,CRT))
	fpLogCur.write('-' * 80 + CRT)
	fpLogCur.write(string.join(qcList, CRT))
	fpLogCur.write('%sTotal: %s%s%s' % (CRT, len(qcList), CRT, CRT))
    
    if isFatal:
	return 1

    return 0

# end writeCuratorLog() -------------------------------

def doDeletes(refsKeys):
    # Purpose: deletes all MGI_Relationships created by this load
    # Returns: 1 if error, else 0
    # Assumes:  database connection exists
    # Effects: queries a database, writes number deleted to curation log
    # Throws: Nothing

    db.sql('''select _StrainMarker_key
	into temporary table toDelete
	from MRK_StrainMarker
	where _Refs_key in (%s)''' % refsKeys, None)

    db.sql('''create index idx1 on toDelete(_StrainMarker_key)''', 'auto')

    results = db.sql('''select count(*) as deleteCt
	from toDelete''', 'auto')

    deleteCt = 0
    for r in results:
	deleteCt = r['deleteCt']

    fpLogCur.write('\nDeleting %s Strain Markers\n\n' % deleteCt)
    db.sql('''delete from MRK_StrainMarker sm
	using toDelete d
	where d._StrainMarker_key = sm._StrainMarker_key''', None)
    db.commit()
    db.useOneConnection(0)

    return 0

# end doDeletes() -------------------------------------

def doBcp():
    # Purpose: executes bcp
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % (bcpin, server, database, strainmarker_table, outputDir, smBcpFile)
    print bcpCmd
    rc = os.system(bcpCmd)
    
    if rc:
	return rc
    if loadOnlyB6 == 'false':
	bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % (bcpin, server, database, acc_table, outputDir, accBcpFile)
	print bcpCmd
	rc = os.system(bcpCmd)
	
	if rc:
	    return rc

	bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % (bcpin, server, database, accref_table, outputDir, accRefBcpFile)
	print bcpCmd
	rc = os.system(bcpCmd)
	
	if rc:
	    return rc

    fpLogCur.write('\nLoaded %s Strain Markers\n\n' % totalLoadedCt)
    fpLogCur.write('Total MGP in input: %s\n\n' % mgpFileCt)
    fpLogCur.write('Total MGP  skipped: %s\n\n' % mgpSkipCt)
    fpLogCur.write('Total MGP Strain Markers Loaded: %s\n\n' % mgpLoadCt)
    fpLogCur.write('Total MGP Strain Markers Loaded with no Marker: %s\n\n' % mgpNoMarkerCt)
    fpLogCur.write('Total B6 Strain Markers Loaded: %s\n\n' % b6LoadedCt)

    fpLogCur.close()
    return 0


# end doBcp() -----------------------------------------

#####################
#
# Main
#
#####################

print '%s' % mgi_utils.date()
print 'running init'
if init() != 0:
    print 'Initialization failed'
    closeFiles()
    sys.exit(1)

# parse MGP input files, write MGP QC and write MGP BCP, unless we are only reloading B6
if loadOnlyB6 == 'false':
    print '%s' % mgi_utils.date()
    print 'running parseMGPFiles'
    if parseMGPFiles() != 0:
	print 'Parsing MGP Files failed'
	closeFiles()
	sys.exit(1)

    # write QC
    print '%s' % mgi_utils.date()
    print 'running writeCuratorLog'
    if writeCuratorLog() != 0:
	print 'Fatal error writing Curator Log - see %s' % curLog
	closeFiles()
	sys.exit(1)

    # write MGP to the bcp files
    print '%s' % mgi_utils.date()
    print 'running writeMGPOutput'
    if writeMGPOutput() != 0:
	print 'Writing MGP Output Files failed'
	closeFiles()
	sys.exit(1)

print '%s' % mgi_utils.date()
print 'running parseB6File'
if parseB6File() != 0:
    print 'Parsing MGI GFF file failed'
    closeFiles()
    sys.exit(1)

print '%s' % mgi_utils.date()
print 'running writeB6Output'
if writeB6Output () != 0:
    print 'Writing B6 Output Files failed'
    closeFiles()
    sys.exit(1)

if QC_ONLY == 'false':

    print '%s' % mgi_utils.date()
    refsKeys = '%s, %s' % (b6RefsKey, mgpRefsKey) # default is B6 and MGP
    if loadOnlyB6 == 'true': # load only B6
	refsKeys = b6RefsKey

    print 'running doDeletes(%s)' % refsKeys

    if doDeletes(refsKeys) != 0:
	print 'Deleting Strain Markers failed'
	sys.exit(1)
	
# close all output files
print '%s' % mgi_utils.date()
print 'running closeFiles()'
if closeFiles() != 0:
    print 'Closing Files failed'
    sys.exit(1)

if  QC_ONLY == 'false':
    # execute bcp
    print '%s' % mgi_utils.date()
    print 'running doBcp()'
    if doBcp() != 0:
	print 'Do BCP failed'
	sys.exit(1)

sys.exit(0)
