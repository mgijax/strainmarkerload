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

#
#  CONSTANTS
#
TAB = '\t'
CRT = '\n'
DATE = mgi_utils.date("%m/%d/%Y")
USAGE='strainmarkerload.py'

#
#  GLOBALS
#

DEBUG = 0

# true if we want to run QC and not load relationships
QC_ONLY = os.environ['QC_ONLY']

# min number of records expected
MIN_RECORDS = int(os.environ['MIN_RECORDS'])

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

# QC reporting data structures
qcDict = {} 
messageMap = {}

# Lookups
strainTranslationLookup = {} # {badName: _Marker_key, ...}
markerLookup = {}            # {MGI ID: Marker ...}
biotypeLookup = []           # lookup raw biotypes from the database

# count of strain markers that will be loaded
loadedCt = 0

# the MGP Strain Marker reference key for J:259852
mgpRefsKey = 282407

# the MGI B6 Strain Marker reference key for J:260092
b6RefsKey =  282660

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
	self.markerID = None
	self.mgpIDs = set()
        
    def toString(self):
        return '%s, %s' % (self.mgiID, self.mgpIDs)
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

    global nextSMKey, nextAccKey, strainTranslationLookup, markerLookup
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
    nextAcc = results[0]['nextAccKey']

    # load qcDict with keys
    qcDict['chr'] = []       # chr is missing 
    qcDict['start'] = []     # startCoordinate is missing 
    qcDict['end'] = []       # endCoordinate is missing 
    qcDict['strand'] = []    # strand is missing 
    qcDict['biotype_m'] = [] # biotype is missing
    qcDict['biotype_u'] = {} # biotype unresolved {biotype:count}
    qcDict['strain_u'] = []  # strain unresolved
    qcDict['mgp'] = []       # mgp is missing
    qcDict['mgi_u'] = []     # marker ID unresolved
    qcDict['mgi_s'] = []     # marker ID secondary
    qcDict['mgi_no'] = []    # marker ID not official
    qcDict['mgi_mgp'] = []   # list of {mgiID: ([set of mpIDs]), ...}, one for each strain file
			     # used to a) determine multiple MGP IDs (within a given strain) per marker
			     # b) write out to bcp files

    messageMap['chr'] = 'Chromosome missing from input'
    messageMap['start'] = 'Start Coordinate missing from input' 
    messageMap['end'] = 'End Coordinate missing from input'
    messageMap['strand'] = 'Strain missing from input'
    messageMap['biotype_m'] = 'Biotype missing from input'
    messageMap['biotype_u'] = 'Biotype from input unresolved'
    messageMap['strain_u'] = 'Strain from input unresolved'
    messageMap['mgp'] = 'MGP ID missing from input'
    messageMap['mgi_u'] = 'Marker from input unresolved'
    messageMap['mgi_s'] = 'Marker from input is secondary ID'
    messageMap['mgi_no'] = 'Marker from input is not official'
    messageMap['mgi_mgp'] = 'Markers from input with > 1 Strain specific MGP ID'  # when reporting multis. 
    #
    # create lookups
    #

    # load lookup of strain translations
    results = db.sql('''select badName, _Object_key as strainKey
	from MGI_Translation
	where _TranslationType_key = 1021''', 'auto')
    for r in results:
	mgpStrain = r['badName']
	strainKey = r['strainKey']
	strainTranslationLookup[mgpStrain] = strainKey

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

    global fpStrainMarkerFile, fpAccFile
    global fpLogCur

    try:
        fpStrainMarkerFile  = open(strainMarkerFile, 'w')
    except:
        print 'Cannot open Strain Marker BCP file: %s' % strainMarkerFile
        sys.exit(1)
    try:
        fpAccFile  = open(accFile, 'w')
    except:
        print 'Cannot open Accession BCP file: %s' % accFile
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

    try:
	fpStrainMarkerFile.close()
	fpLogCur.close()
    except:
	return 1
    return 0

# end closeFiles() -------------------------------

def parseFiles( ): 
    # Purpose: parses input files looping, opening and closing
    # Returns: 1 if error, else 0
    # Assumes: output file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing
   
    global qcDict

    inputFileList = string.split(inputFileString, ' ')
    #first = 1
    for file in inputFileList:
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
        
        # add the last strainMarkerInput strain file to the qcDict
        #if first == 0:
        #   qcDict['mgi_mgp'].append(strainMarkerInput)
        #if first == 1:
        #    first = 0

	# build this as we parse each file - adding mgpIDs to strainMarkerObject if mgiID
	# found > 1 in file
	strainMarkerDict = {}   # {mgiID:strainMarkerObject, ...}

	# after file parsed copy the strainMarkerObjects from strainMarkerDict to this mapping by strain
	strainMarkerInput = {} 		# {strain: list of strainMarkerObjects, ...}
	strainMarkerInput[strain] = []	# initialize 

	# iterate thru lines in this strain file
	for line in fpIn.readlines():
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
	    if chr == '':
		qcDict['chr'].append(line)
		isSkip = 1
	    if start == '':
		qcDict['start'].append(line)
		isSkip = 1
  	    if end == '':
		qcDict['end'].append(line)
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
            #  check that biotype is in the database
            if biotype.lower().strip() not in biotypeLookup:
		if  biotype not in qcDict['biotype_u']:     
		    qcDict['biotype_u'][biotype] = 1
		else:
		    qcDict['biotype_u'][biotype] += 1
                #print '%s BIOTYPE NOT IN MGI: %s' %(mgpID, biotype)
# UNCOMMENT THIS AFTER DEV!!
		isSkip = 1

	    # resolve strain with translation lookup
	    if strain not in strainTranslationLookup:
		qcDict['strain_u'].append(strain)
		isSkip = 1
	    # Resolve MGI ID 
	    if mgiID != '':
		if mgiID not in markerLookup:
		    qcDict['mgi_u'].append('%s : %s' % (mgiID, line))
		    #print '%s NOT IN MGI in strain file: %s' % (mgiID, strain)
		    isSkip = 1
		else:
		    marker = markerLookup[mgiID]
		    markerKey = marker.markerKey # load official markers
		    
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
	    
	    if not isSkip:
		strainMarkerObject = StrainMarker() # default to a new one
	        if mgiID not in strainMarkerDict:
		    strainMarkerObject.markerID = mgiID
		    strainMarkerObject.mgpIDs.add(mgpID)
		else:
		    strainMarkerObject = strainMarkerDict[mgiID]
                    strainMarkerObject.mgpIDs.add(mgpID)
		strainMarkerDict[mgiID] =  strainMarkerObject

	# copy strainMarker objects from strainMarkerDict to strainMarkerInput
	for mgiID in strainMarkerDict:
	    strainMarkerInput[strain].append(strainMarkerDict[mgiID])

	# add to the qcDict
        qcDict['mgi_mgp'].append(strainMarkerInput)

    return 0

# end parseFiles() -------------------------------------
def writeBcp():
    # Purpose: writes to BCP files if there are no errors
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing
    strainMarkerInputList = qcDict['mgi_mgp']
    for strainMarkerInputDict in strainMarkerInputList:
        for strain in strainMarkerInputDict:
            strainMarkerObjectsList = strainMarkerInputDict[strain]
            for strainMarkerObject in strainMarkerObjectsList:
                mgpList = list(strainMarkerObject.mgpIDs) # get the list of mgp IDs
                if strainMarkerObject.markerID != '' and len(mgpList) > 1:
		    continue  # these were reported
	        # otherwise write out to bcp file
		mgiID = strainMarkerObject.markerID 
		    #for mgpID in mgpList:
			# write out to bcp files one per mgpID, there will be one mgpID for all but the strain markers
			# with no canonical gene
	#START HERE: 4/13/18 5:318 pm
  	    #print 'writeBcp mgiID: %s mpgID: %s' % (mgiID, mpID)
    return 0

# end writeBcp() ---------------------------------------------------

def writeCuratorLog():
    # Purpose: writes QC errors to the curator log
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    #
    # special handling for MGI IDs mapped to multiple IDs for a given strain. 
    # Each dict in mgiMgpList represents the ids of a strain 
    #
    count = 0
    multiList = []
    strainMarkerInputList = qcDict['mgi_mgp']
    for strainMarkerInputDict in strainMarkerInputList:
	for strain in strainMarkerInputDict:
	    print 'strain: %s' % strain
	    strainMarkerObjectsList = strainMarkerInputDict[strain]
	    for strainMarkerObject in strainMarkerObjectsList:
		mgpList = list(strainMarkerObject.mgpIDs) # get the list of mgp IDs
		#print 'strainMarkerObject.markerID: %s mgpList: %s' % (strainMarkerObject.markerID, mgpList)
		if strainMarkerObject.markerID != '' and len(mgpList) > 1:  # if > 1 ID, report
		    #print 'mgpList: %s' % mgpList
		    count += 1
		    multiList.append('%s: %s' % (strainMarkerObject.markerID, string.join(mgpList, ', ')))

    # if we have multiple mpg IDs/marker within a strain, report
    if len(multiList):
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
    #  process remaining QC, we have removed the special handling cases from qcDict
    #
    for key in qcDict:
	if  key in ('mgi_mgp', 'biotype_u'):
	    continue    # we already processed these above, mgi_mgp needed by writeBcp() 

        qcList = qcDict[key]
	if qcList == []:
	    continue
	msg = messageMap[key]
	fpLogCur.write('%s%s'% (msg,CRT))
	fpLogCur.write('-' * 80 + CRT)
	fpLogCur.write(string.join(qcList, CRT))
	fpLogCur.write('%sTotal: %s%s%s' % (CRT, len(qcList), CRT, CRT))

    return 0

# end writeCuratorLog() -------------------------------

def doDeletes():
    # Purpose: deletes all MGI_Relationships created by this load
    # Returns: 1 if error, else 0
    # Assumes:  database connection exists
    # Effects: queries a database, writes number deleted to curation log
    # Throws: Nothing

    results = db.sql('''select count(*) as deleteCt 
	from MGI_Relationship 
	where _CreatedBy_key = %s ''' % userKey, 'auto')
    deleteCt = 0
    for r in results:
	deleteCt = r['deleteCt']
    fpLogCur.write('\nDeleting %s Relationships\n\n' % deleteCt)
    db.sql('''delete from MGI_Relationship where _CreatedBy_key = %s ''' % userKey, None)
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

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % (bcpin, server, database, table, outputDir, bcpFile)
    rc = os.system(bcpCmd)
    return rc

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

# parse the input files
print '%s' % mgi_utils.date()
print 'running parseFiles'
if parseFiles() != 0:
    print 'Parsing Files failed'
    closeFiles()
    sys.exit(1)

# write QC
print '%s' % mgi_utils.date()
print 'running writeCuratorLog'
if writeCuratorLog() != 0:
    print 'Writing the Curator Log failed'
    closeFiles()
    sys.exit(1)

# write to the bcp files
print '%s' % mgi_utils.date()
print 'running writeBcp'
if writeBcp() != 0:
    print 'Writing to the BCP Files failed'
    closeFiles()
    sys.exit(1)

#if QC_ONLY == 'false':
#    # delete existing relationships
#    print '%s' % mgi_utils.date()
#    print 'running doDeletes()'
#    if doDeletes() != 0:
#	print 'Do Deletes failed'
#	sys.exit(1)

# close all output files
print '%s' % mgi_utils.date()
print 'running closeFiles()'
if closeFiles() != 0:
    print 'Closing Files failed'
    sys.exit(1)

#if  QC_ONLY == 'false':
#    # execute bcp
#    print '%s' % mgi_utils.date()
#    print 'running doBcp()'
#    if doBcp() != 0:
#	print 'Do BCP failed'
#	sys.exit(1)

sys.exit(0)
