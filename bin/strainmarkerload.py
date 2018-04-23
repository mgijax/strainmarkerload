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

# MGP logicalDB Name
logicalDBKey = 209

# MRK_StrainMarker MGI Type
mgiTypeKey = 44

# if true do not delete/reload
isFatal = 0

# input
infileDir = os.environ['INPUTDIR']
inputFileString =  os.environ['INPUT_DIR_LIST']

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

# QC reporting data structures
qcDict = {} 
messageMap = {}

# Lookups
strainTranslationLookup = {} # {badName: _Strain_key, ...}
markerLookup = {}            # {MGI ID: Marker ...}
chrLookup = {}		     # {Mouse, Laboratory chr: chrKey, ...}
biotypeLookup = []           # lookup raw biotypes from the database

# the MGP Strain Marker reference key for J:259852
mgpRefsKey = 282407

# the MGI B6 Strain Marker reference key for J:260092
b6RefsKey =  282660

# the MGI C57BL/6J strain key
b6StrainKey = 38048

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
b6LoadCt = 0
mgpFileCt = 0
mgpLoadCt = 0
mgpSkipCt = 0
mgpNoMarkerCt = 0
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
        
    def toString(self):
        return '%s, %s, %s, %s' % (self.mgiID, self.markerKey, self.strainKey, self.mgpIDs)
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
    global biotypeLookup, messageMap
   
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

    # load lookup of rawbiotypes
    results = db.sql('''select distinct t.term
	from MRK_BiotypeMapping b, VOC_Term t
	where b._BiotypeVocab_key in (103, 104)
	and b._BiotypeTerm_key = t._Term_key''', 'auto')
    for r in results:
	biotypeLookup.append(r['term'].lower())
    #print 'biotypeLookup: %s' % biotypeLookup

    return 0

# end init() -------------------------------

def openFiles ():
    # Purpose: Open input/output files.
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened, 
    #  creates files in the file system

    global fpStrainMarkerFile, fpAccFile, fpAccRefFile, fpLogCur

    try:
        fpStrainMarkerFile  = open(strainMarkerFile, 'w')
    except:
        print 'Cannot open Strain Marker BCP file: %s' % strainMarkerFile
        sys.exit(1)
    try:
        fpAccFile  = open(accFile, 'w')
    except:
        print 'Cannot open ACC_Accession bcp  file: %s' % accFile
        sys.exit(1)
    try:
        fpAccRefFile  = open(accRefFile, 'w')
    except:
        print 'Cannot open ACC_AccessionReference bcp file: %s' % accRefFile
        sys.exit(1)

    try:
	fpLogCur = open(curLog, 'a')
    except:
	print 'Cannot open Curator Log file: %s' % curLog
        sys.exit(1)

    return 0

# end openFiles() -------------------------------

def closeFiles ():
    # Purpose: Close all file descriptors
    # Returns: 1 if error, else 0
    # Assumes: all file descriptors were initialized
    # Effects: Nothing
    # Throws: Nothing
    global fpStrainMarkerFile, fpAccFile
    try:
	fpStrainMarkerFile.close()
	fpAccFile.close()
        fpAccRefFile.close()
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
	    markerKey = ''
	    biotype = ''
	    for t in tokens2:
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
            if biotype.lower().strip() != '' and biotype.lower().strip() not in biotypeLookup:
		if  biotype not in qcDict['biotype_u']:     
		    qcDict['biotype_u'][biotype] = 1
		else:
		    qcDict['biotype_u'][biotype] += 1
		isSkip = 1

	    # Resolve MGI ID 
	    if mgiID != '':
		if mgiID not in markerLookup:
		    qcDict['mgi_u'].append('%s : %s' % (mgiID, line))
		    #print '%s NOT IN MGI in strain file: %s' % (mgiID, strain)
		    isSkip = 1
		else:
		    marker = markerLookup[mgiID]
		    markerKey = marker.markerKey 
		    
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

		strainMarkerObject = StrainMarker() # default to a new one
	        if mgiID not in strainMarkerDict:
		    strainMarkerObject.markerID = mgiID
		    strainMarkerObject.markerKey = markerKey
		    strainMarkerObject.strainKey = strainKey
		    strainMarkerObject.mgpIDs.add(mgpID)
		else:
		    strainMarkerObject = strainMarkerDict[mgiID]
                    strainMarkerObject.mgpIDs.add(mgpID)
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
    # Purpose: writes to Accession, AccessionReference & StrainMarker BCP file if there are no errors
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
                if strainMarkerObject.markerID != '' and len(mgpList) > 1:
		    mgpSkipCt += len(mgpList)
		    continue  # these were reported
	        # otherwise write out to bcp file
		mgiID = strainMarkerObject.markerID 
		markerKey = strainMarkerObject.markerKey
		strainKey = strainMarkerObject.strainKey
		for mgpID in mgpList: # One except for the no marker
		    totalLoadedCt += 1
		    fpStrainMarkerFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextSMKey, TAB, strainKey, TAB, markerKey, TAB, mgpRefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

		    prefixPart, numericPart = accessionlib.split_accnum(mgpID)

		    fpAccFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s0%s1%s%s%s%s%s%s%s%s%s' \
		    % (nextAccKey, TAB, mgpID, TAB, prefixPart, TAB, numericPart, TAB, logicalDBKey, TAB, nextSMKey, TAB, mgiTypeKey, TAB, TAB, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))
		    fpAccRefFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' \
		    % (nextAccKey, TAB, mgpRefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))
		    nextSMKey += 1
		    nextAccKey += 1

    return 0

# end writeMGPOutput() ---------------------------------------------------

def writeB6Output():
    # Purpose: writes to Accession, AccessionReference & StrainMarker BCP file if there are no errors
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global nextSMKey, totalLoadedCt, b6LoadCt

    results = db.sql('''select distinct _Marker_Key
			from MRK_Location_Cache lc
			where lc.provider in ('NCBI', 'Ensembl')''', 'auto')
    print 'len(results): %s' % len(results)
    for r in results:
	markerKey = r['_Marker_key']

	fpStrainMarkerFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextSMKey, TAB, b6StrainKey, TAB, markerKey, TAB, b6RefsKey, TAB, userKey, TAB, userKey, TAB, loaddate, TAB, loaddate, CRT))

	nextSMKey += 1
    totalLoadedCt += len(results)
    b6LoadCt = len(results)
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
    print 'rc: %s' % rc
    if rc:
	return rc
    if loadOnlyB6 == 'false':
	bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % (bcpin, server, database, acc_table, outputDir, accBcpFile)
	print bcpCmd
	rc = os.system(bcpCmd)
	print 'rc: %s' % rc
	if rc:
	    return rc

	bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % (bcpin, server, database, accref_table, outputDir, accRefBcpFile)
	print bcpCmd
	rc = os.system(bcpCmd)
	print 'rc: %s' % rc
	if rc:
	    return rc

    fpLogCur.write('\nLoaded %s Strain Markers\n\n' % totalLoadedCt)
    fpLogCur.write('Total MGP in input: %s\n\n' % mgpFileCt)
    fpLogCur.write('Total MGP  skipped: %s\n\n' % mgpSkipCt)
    fpLogCur.write('Total MGP Strain Markers Loaded: %s\n\n' % mgpLoadCt)
    fpLogCur.write('Total MGP Strain Markers Loaded with no Marker: %s\n\n' % mgpNoMarkerCt)
    fpLogCur.write('Total B6 Strain Markers Loaded: %s\n\n' % b6LoadCt)

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
