#!/usr/local/bin/python
#
#  mp_emapaload.py
###########################################################################
#
#  Purpose:
#
#      Validate input and create feature relationships bcp file
#
#  Usage:
#
#      mp_emapaload.py 
#
#  Inputs:
#
#	1. MP OWL file, UBERON and EMAPA OBO files
#
#	2. Configuration - see mp_emapaload.config
#
#  Outputs:
#
#       1. MGI_Relationship.bcp
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
#      7) Delete existing relationships
#      8) BCP in new relationships:
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
#  11/03/2017  sc  Initial development
#
###########################################################################

import sys
reload(sys)
sys.setdefaultencoding('utf8')

import os
import string
import Set

import db
import mgi_utils

from xml.dom.minidom import parse
from xml import *
#
#  CONSTANTS
#
TAB = '\t'
CRT = '\n'
DATE = mgi_utils.date("%m/%d/%Y")
USAGE='mp_emapaload.py'

#
#  GLOBALS
#

# if true, write out to mp, uberon and emapa files
DEBUG = 0

# true if we want to run QC and not load relationships
QC_ONLY = os.environ['QC_ONLY']

# min number of records expected
MIN_RECORDS = int(os.environ['MIN_RECORDS'])

# input files
inMP = os.environ['INPUT_FILE_MP']
inUBERON = os.environ['INPUT_FILE_UBERON']
inEMAPA = os.environ['INPUT_FILE_EMAPA']

# input file descriptors
fpUin = ''
fpEin = ''

# output files
outMP = os.environ['OUTPUTFILE_MTOU']
outUberon =  os.environ['OUTPUTFILE_UTOE']
outEmapa =  os.environ['OUTPUTFILE_EMAPA']
curLog = os.getenv('LOG_CUR')

# output file descriptors
# we will use these for debugging
fpMtoU = ''
fpUtoE = ''
fpEmapa = ''

# output bcp files
bcpFile =   os.environ['RELATIONSHIP_BCP']
outputDir = os.environ['OUTPUTDIR']
relationshipFile = '%s/%s' % (outputDir, bcpFile)
fpRelationshipFile = ''

#
# For sanity/QC checks
#

# reporting data structures
mpNotInDatabase = []
emapaNotInDatabase = []
mpNoEmapa = []
obsAltUberonInMP = []
obsAltEmapaInUberon = []
oneMpMultiUberon = []
oneUberonMultiEmapa = []
someValuesFromLost = []  # list of emapaIds lost because we are using 
		     # Description field instead of someValuesFrom
# MP to UBERON from MP file {mp:MP relationship object, ...}
mpDict = {}

# UBERON to EMAPA from UBERON file {uberon:uberon relationship object, ...}
uberonDict = {} 

# EMAPA dict from file )emapaId:emapa relationship object with null id2
emapaDict = {}

# list of alt_id values from the uberon file
uberonAltIdList = []

# list of alt_id values from the emapa file
emapaAltIdList = []

# Lookups
mpLookup = {} # {mpID:[termKey, isObsolete], ...}
emapaLookup = {} # {emapaId:[termKey, isObsolete], ...}

# count of relationships that will be loaded
loadedCt = 0
distinctMpLoaded = set([])
distinctEmapaLoaded =  set([])
#
# For loading relationships
#

# The mp_emapa relationship category key 'mp_to_emapa'
catKey = 1007

# the mp_emapa relationship term key 'mp_to_emapa'
relTermKey = 37085930

# the mp_emapa qualifier key 'Not Specified'
qualKey = 11391898

# the mpo_emapa evidence key 'Not Specified'
evidKey = 17396909

# the mp_emapa reference key 'J:247341'
refsKey = 257191

# mp_emapaload user key
userKey = 1576

# database primary keys, will be set to the next available from the db
nextRelationshipKey = None	# MGI_Relationship._Relationship_key

# for bcp
bcpin = '%s/bin/bcpin.csh' % os.environ['PG_DBUTILS']
server = os.environ['MGD_DBSERVER']
database = os.environ['MGD_DBNAME']
table = 'MGI_Relationship'

class Relationship:
    # Is: data object for relationship between two vocab terms
    #     this is used by the QC checks
    # Has: a set of term attributes
    # Does: provides direct access to its attributes
    #
    def __init__ (self):
        # Purpose: constructor
        # Returns: nothing
        # Assumes: nothing
        # Effects: nothing
        # Throws: nothing
        self.id1 = None # from input file
	self.preferred = None # from database
        self.term = None # from input file
	self.termKey = None # null for uberon, null if not in database for mp/emapa
        self.isObsolete = None # from uberon file or db for emapa/mp
        self.id2 = [] # list of IDs, null for emapa
	self.id3 = [] # for uberon only - From owl file <rdf:Description rdf:about
# end class Relationships -----------------------------------------

def checkArgs ():
    # Purpose: Validate the arguments to the script.
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: exits if unexpected args found on the command line
    # Throws: Nothing

    if len(sys.argv) != 1:
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

    global nextRelationshipKey, mpLookup, emapaLookup
   
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
    # get next MGI_Relationship key
    #
    results = db.sql('''select max(_Relationship_key) + 1 as nextKey
	    from MGI_Relationship''', 'auto')
    if results[0]['nextKey'] is None:
	nextRelationshipKey = 1000
    else:
	nextRelationshipKey = results[0]['nextKey']

    #
    # create lookups
    #
    # lookup of MP terms
    results = db.sql('''select a.accid, a.preferred, t.term, t.isObsolete, t._Term_key
        from VOC_Term t, ACC_Accession a
        where t._Vocab_key = 5
        and t._Term_key = a._Object_key
        and a._MGIType_key = 13
        and a._LogicalDB_key = 34
        and a.preferred = 1''', 'auto')

    for r in results:
        #mpId = string.lower(r['accid'])
	mpId = r['accid']
        termKey = r['_Term_key']
	isObsolete = r['isObsolete']
	preferred = r['preferred']
        mpLookup[mpId] = [termKey, isObsolete, preferred]

    # load lookup of EMAPA terms
    results = db.sql('''select a.accid, a.preferred, t.term, t.isObsolete, t._Term_key
        from VOC_Term t, ACC_Accession a
        where t._Vocab_key = 90
        and t._Term_key = a._Object_key
        and a._MGIType_key = 13
        and a._LogicalDB_key = 169''', 'auto')

    for r in results:
        #emapaId = string.lower(r['accid'])
	emapaId = r['accid']
        termKey = r['_Term_key']
        isObsolete = r['isObsolete']
	preferred = r['preferred']
        emapaLookup[emapaId] = [termKey, isObsolete, preferred]

    return 0

# end init() -------------------------------

def openFiles ():
    # Purpose: Open input/output files.
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened, 
    #  creates files in the file system

    global fpUin, fpEin, fpMtoU, fpUtoE, fpEmapa, fpRelationshipFile
    global fpLogCur

    try:
        fpUin = open(inUBERON, 'r')
    except:
        print 'Cannot open UBERON input file: %s' % inUBERON
        sys.exit(1)
    try:
	fpEin = open(inEMAPA, 'r')
    except:
        print 'Cannot open EMAPA input file: %s' % inEMAPA
        sys.exit(1)
    try:
        fpMtoU = open(outMP, 'w')
    except:
        print 'Cannot open MP/Uberon output file: %s' % inUBERON
        sys.exit(1)
    try:
        fpUtoE = open(outUberon, 'w')
    except:
        print 'Cannot open Uberon/EMAPA output file: %s' % inUBERON
        sys.exit(1)
    try:
        fpEmapa = open(outEmapa, 'w')
    except:
        print 'Cannot open EMAPA output file: %s' % inUBERON
        sys.exit(1)

    try:
        fpRelationshipFile = open(relationshipFile, 'w')
    except:
        print 'Cannot open Feature relationships bcp file: %s' % relationshipFile
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
	fpUin.close()
	fpEin.close()
	fpMtoU.close() 
	fpUtoE.close()
	fpEmapa.close()

	fpRelationshipFile.close()
	fpLogCur.close()
    except:
	return 1
    return 0

# end closeFiles() -------------------------------

def parseFiles( ): 
    # Purpose: parses input files into data structures 
    #	(and intermediate files representing parsed input files when in 
    #    DEBUG mode)
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    try:
	parseMPFile()
	parseUberonFile()
	parseEmapaFile()
    except:
	return 1
    return 0

# end parseFiles() -------------------------------------

def parseMPFile():
    # Purpose: parses MP input file into data structures, doing QC on the input
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    global mpDict, mpNotInDatabase

    # --Parse the MP owl file
    # the parse function can take an open file object or a file name
    try:
	dom = parse(inMP) # dom is the complete file loaded into memory, you can view it via dom.toxml()
    except:
	# sanity #1
	fpLogCur.write('MP OWL file not in correct format, mp_emapaload failed\n')
	sys.exit('MP OWL file not in correct format, mp_emapaload failed')
    recordCt = 0
    for e in dom.getElementsByTagName('owl:Class'): # iterate over the NodeList, e is a single Node
 	recordCt += 1
	rel = Relationship()
	mpID = ''
	mpTerm = ''
	mp = e.getAttribute('rdf:about') # the line with the MP URL
	if string.find(mp, 'MP_') == -1:
	    continue
	mpID = string.split(mp, '/')[-1]
	mpID = mpID.replace('_', ':')
	termNodes = e.getElementsByTagName('rdfs:label') # get the node list that contains the MP term
	mpTerm = ''
	try:
	    mpTerm = termNodes[0].firstChild.data
 	except:
	    mpTerm = termNodes[1].firstChild.data
	rel.id1 = mpID
	rel.term = mpTerm
	if mpID in mpLookup:
	    termKey, isObsolete, preferred = mpLookup[mpID]
	    rel.termKey = termKey
	    rel.isObsolete = isObsolete
	    rel.preferred = preferred
	else:
	    # report #3
	    mpNotInDatabase.append('%s %s' % (mpID, mpTerm))
	    continue

	# uberon IDs are found in two places in an owl stanza
   	# a) If 1 only (rdf:Description rdf:about), use 1
	# b) If 2 only (owl:someValuesFrom rdf:resource), use 2
	# c) If both 1 and 2 (Description & resource) use 1
	uberonDescNodes =  e.getElementsByTagName('rdf:Description')

    	for n in uberonDescNodes:
	    url =  n.getAttribute('rdf:about')
	    if string.find(url, 'UBERON_') != -1:
                uberonDescID =  string.split(url, '/')[-1]
                uberonDescID =  uberonDescID.replace('_', ':')
		if uberonDescID not in rel.id3: # don't want duplicates
		    rel.id3.append(uberonDescID)
	uberonSomeNodes = e.getElementsByTagName('owl:someValuesFrom') # node list containing uberon id
	for n in uberonSomeNodes:
	    url = n.getAttribute('rdf:resource')
	    if string.find(url, 'UBERON_') != -1:
		uberonSomeID =  string.split(url, '/')[-1]
		uberonSomeID =  uberonSomeID.replace('_', ':')
		if uberonSomeID not in rel.id2: # don't want duplicates
		    rel.id2.append(uberonSomeID)
	mpDict[mpID] = rel
    
    if len(mpDict) < MIN_RECORDS:
	# sanity #2
	fpLogCur.write('MP File has less than the configured minimum records: %s, number of records: %s\n' % (MIN_RECORDS, len(mpDict)))
        sys.exit('MP File has less than the configured minimum records: %s, number of records: %s' % (MIN_RECORDS, len(mpDict)))

    # write out to file for DEBUG
    if DEBUG:
	print 'in DEBUG'
	fpMtoU.write('MP ID%spreferred%smpTerm%sisObsolete%sUberon Some IDs%sUberon Desc IDs%s' % (TAB, TAB, TAB, TAB, TAB, CRT))
        keys = mpDict.keys()
        keys.sort()
        for key in keys:
	    print 'key: %s' % key
	    rel = mpDict[key]
	    mpID = rel.id1
   	    print 'mpID: %s' % mpID	
	    preferred = rel.preferred
	    mpTerm = rel.term
	    isObsolete = rel.isObsolete
	    uberonDescIds = string.join(rel.id3, ', ')
	    print 'uberonDescIds: %s' % uberonDescIds
	    uberonSomeIds = string.join(rel.id2, ', ')
	    print 'uberonSomeIds: %s' % uberonSomeIds
	    print '%s%s%s%s%s%s%s%s%s%s%s%s' % (mpID, TAB, preferred, TAB, mpTerm, TAB, isObsolete, TAB, uberonSomeIds, TAB, uberonDescIds, CRT)
	    fpMtoU.write('%s%s%s%s%s%s%s%s%s%s%s%s' % (mpID, TAB, preferred, TAB, mpTerm, TAB, isObsolete, TAB, uberonSomeIds, TAB, uberonDescIds, CRT))

    return 0

# end parseMPFile() -------------------------------------

def parseUberonFile():
    # Purpose: parses Uberon input file into data structures, doing QC on the
    # input 
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    global uberonDict

    # values for parsing
    nameValue = 'name:'
    altIdValue = 'alt_id:'
    recordFound = 0
    uberonIdValue = 'id: UBERON:'
    emapaXrefValue = 'xref: EMAPA:'

    emapaList = [] # list of emapa Ids for a given uberon term
    uberonId = ''
    uberonName = ''
    isObsolete = 0
    firstRecord = 1
    lines = fpUin.readlines()
    if string.strip(lines[0]) != 'format-version: 1.2':
	# sanity #1
	fpLogCur.write('Uberon OBO file not in correct format, mp_emapload failed\n')
	sys.exit('Uberon OBO file not in correct format, mp_emapload failed')
    recordCt = 0	
    for line in lines:
	if string.find(line,'[Term]') == 0:
	    recordCt += 1
	    recordFound = 0 
	    if firstRecord != 1:
		rel = Relationship()
		rel.id1 = uberonId
		rel.term = uberonName
		rel.isObsolete = isObsolete
		rel.id2 = emapaList    
		uberonDict[uberonId] = rel
	    if firstRecord == 1:
		firstRecord =  0
	elif line[:11] == uberonIdValue:
	    recordFound = 1
	    emapaList = []
	    uberonId = line[4:-1]
	elif recordFound and line[:5] == nameValue:
	    uberonName = line[6:-1]
	    if uberonName.find('obsolete') == 0:
		isObsolete = 1
	    else:
		isObsolete = 0
	elif recordFound and line[:12] == emapaXrefValue:

	    emapaId = line[6:-1]
	    # one record like this: xref: EMAPA:35128 {source="MA"}
	    emapaId = emapaId.split(' ')[0]
	    if emapaId not in emapaList: # don't want duplicates
		emapaList.append(emapaId)
        elif recordFound and line[:7] == altIdValue:
	    altId = line[8:-1] 
	    uberonAltIdList.append(altId)
	else: # we don't care about this line, go to the next line
	    continue
    if len(uberonDict) < MIN_RECORDS:
	# sanity #2
	fpLogCur.write('UBERON File has less than the configured minimum records: %s, number of records: %s\n' % (MIN_RECORDS, len(uberonDict)))
	sys.exit('UBERON File has less than the configured minimum records: %s, number of records: %s' % (MIN_RECORDS, len(uberonDict)))
    if DEBUG:
	keys = uberonDict.keys()
	keys.sort()
	fpUtoE.write('Uberon ID%sUberon Term%sisObsolete%sEmapa IDs%s' % (TAB, TAB, TAB, CRT))
	for key in keys:
	    rel = uberonDict[key]
	    uberonID = rel.id1
	    uberonTerm = rel.term
	    isObsolete = rel.isObsolete
	    emapaIds = string.join(rel.id2, ', ')
	    fpUtoE.write('%s%s%s%s%s%s%s%s' % (uberonID, TAB, uberonTerm, TAB, isObsolete, TAB, emapaIds, CRT))

    return 0

# end parseUberonFile() -------------------------------------

def parseEmapaFile():
    # Purpose: parses EMAPA input file into data structures, doing QC on the
    # input
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    global emapaDict
    
    emapaIdValue = 'id: EMAPA:'
    altIdValue = 'alt_id:'

    recordFound = 0
    nameValue = 'name:'

    lines = fpEin.readlines()
    if string.strip(lines[0]) != 'format-version: 1.2':
	# sanity #1
	fpLogCur.write('EMAPA OBO file not in correct format, mp_emapload failed\n')
        sys.exit('EMAPA OBO file not in correct format, mp_emapaload failed')

    recordCt = 0
    for line in lines:
	if line == '[Term]':
	    recordCt += 1
	    recordFound = 0

	elif line[:10] == emapaIdValue:
	    emapaId = line[4:-1]
	    recordFound = 1
	elif recordFound and line[:7] == altIdValue:
            altId = line[8:-1]
            emapaAltIdList.append(altId)
	elif recordFound and line[:5] == nameValue:
	    emapaTerm = line[5:-1]
	    rel = Relationship()
	    rel.id1 = emapaId
	    rel.term = emapaTerm
	    termKey = 0
	    isObsolete = 0
	    preferred = 0
	    if emapaId in emapaLookup:
		termKey, isObsolete,preferred = emapaLookup[emapaId]

	    rel.termKey = termKey
	    rel.isObsolete = isObsolete
	    rel.preferred = preferred
	    emapaDict[emapaId] = rel
    if len(emapaDict) < MIN_RECORDS:
	# sanity #2
	fpLogCur.write('EMAPA File has less than the configured minimum records: %s, number of records: %s\n' % (MIN_RECORDS, len(emapaDict)))
        sys.exit('EMAPA File has less than the configured minimum records: %s, number of records: %s' % (MIN_RECORDS, len(emapaDict)))
    if DEBUG:
        keys = emapaDict.keys()
        keys.sort()
	fpEmapa.write('Emapa ID%spreferred%sEmapa Term%sisObsolete%s' % (TAB, TAB, TAB, CRT))
        for key in keys:
	    rel = emapaDict[key]
	    emapaId = rel.id1
	    preferred = rel.preferred
            emapaTerm = rel.term
            isObsolete = rel.isObsolete
	    fpEmapa.write('%s%s%s%s%s%s%s%s' % (emapaId, TAB, preferred, TAB, emapaTerm, TAB, isObsolete, CRT))
    
    return 0

# end parseEMapaFile() -------------------------------------

def findRelationships():
    # Purpose: finds the transitive relationship between MP and EMAPA via Uberon
    # and writes them to a bcp file
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    global loadedCt, nextRelationshipKey, distinctMpLoaded, distinctEmapaLoaded 
    # iterate thru the MP records and get their Uberon associations
    for mpId in mpDict:
	mpRel = mpDict[mpId]
	mpTerm = mpRel.term
	mpTermKey = mpRel.termKey
	uberonSomeList = mpRel.id2
	uberonDescList = mpRel.id3

	# The chosen list of uberon IDs from mpID
	uberonList = []

	if len(uberonDescList) and len(uberonSomeList):
	    uberonList = uberonDescList # if we have uberon for both choose Description
	    # report #10 - In cases where an mp stanza has both Description Uberon 
	    # ID(s) and someValueFrom Uberon ID(s), report the cases where the 
	    # Description Uberon ID(s) have no EMAPA association, but the 
	    # someValueFrom Uberon ID(s) do have EMAPA associations. In other words - 		 
	    # report the mappings to EMAPA that we are loosing by new requirements.

	    # check the chosen list of Uberon IDs for EMAPA associations
	    
	    # flag to report EMAPA ID mappings from uberonSomeList
	    reportSome = 0 # default
	    for ubId in uberonDescList:
		# If no emapa associations in uberonDescList - 
		# set flag to report EMAPA IDs from uberonSomeList
		if ubId in uberonDict and uberonDict[ubId].id2 == []:
		    reportSome = 1
		if reportSome:
		    for ubId in uberonSomeList:
			if ubId in uberonDict:
			    emapaList = uberonDict[ubId].id2
			    if emapaList != []:
				msg = '%s %s %s %s' % (mpId, mpTerm, ubId, string.join(emapaList, ', '))
				if msg not in someValuesFromLost:
				    someValuesFromLost.append(msg)

	elif len(uberonDescList):
	    uberonList = uberonDescList # if we have just Description
	else:
	    uberonList = uberonSomeList # otherwise choose someValuesFrom (which may be empty)
	if len(uberonList) > 1:
	    # report #8 and load
	    termList = []
	    for u in uberonList:
		uTerm = 'not in Uberon file' # default
		if u in uberonDict:
		    uTerm = uberonDict[u].term
		elif u in uberonAltIdList:
		    uTerm = 'altId'
		termList.append('%s (%s)' % (u, uTerm))
	    msg = '%s (%s) %s' % (mpId, mpTerm, string.join(termList, ', '))
	    if msg not in oneMpMultiUberon:
		oneMpMultiUberon.append(msg)
	# For each uberon ID get the EMAPA Ids it maps to, if any
	for ubId in uberonList:
	    ubTerm = 'not in Uberon file' # default
	    if ubId in uberonDict:
		ubTerm = uberonDict[ubId].term
	    if ubId in uberonAltIdList:
                # report and skip #6 ubid is alt_id
		msg = '%s (%s) ubId is alternate id: %s (%s)' % (mpId, mpTerm, ubId, ubTerm)
		if msg not in  obsAltUberonInMP:
		    obsAltUberonInMP.append(msg)
	    elif ubId in uberonDict and uberonDict[ubId].isObsolete:
                # report and skip #6 and skip, is obsolete
		msg = '%s (%s) ubId is obsolete: %s (%s)' % (mpId, mpTerm, ubId, ubTerm)
		if msg not in obsAltUberonInMP:
		    obsAltUberonInMP.append(msg)
	    elif ubId not in uberonDict:
	 	# report and skip #5 MP that don't map to emapa - no mp to uberon
		msg = '%s %s %s' % (mpId, mpTerm, ubId)
		if msg not in mpNoEmapa:
		    mpNoEmapa.append(msg)
	    else:
		uberonRel = uberonDict[ubId]
		ubTerm = uberonRel.term
		emapaList = uberonRel.id2
		# report and skip #5 no Uberon to EMAPA xref
		if emapaList == []:
		    msg = '%s %s %s %s' % (mpId, mpTerm, ubId, ubTerm)
		    if msg not in mpNoEmapa:
			mpNoEmapa.append(msg)
		if len(emapaList) > 1:
		    # report and load #9 report
		    termList = []
		    for e in emapaList:
			eTerm = 'not in EMAPA file'
			if e in emapaDict:
			    eTerm = emapaDict[e].term
		 	elif  e in emapaAltIdList:
			    eTerm = 'altId'
			termList.append('%s (%s)' % (e, eTerm))

		    msg = '%s (%s): %s' % (ubId, ubTerm, string.join(termList, ', '))
		    if msg not in oneUberonMultiEmapa:
			oneUberonMultiEmapa.append(msg)
		for emapaId in emapaList:
		    emapaTerm = 'altId' # default
		    if emapaId in emapaDict:
			emapaTerm =  emapaDict[emapaId].term

		    if emapaId in emapaAltIdList:
                        # report and skip #7 emapa is alt_id
			msg = '%s (%s) emapaId is alternate id: %s (%s)' % (ubId, ubTerm, emapaId, emapaTerm)
			if msg not in obsAltEmapaInUberon:
			    obsAltEmapaInUberon.append(msg)
		    elif emapaId in emapaDict and emapaDict[emapaId].isObsolete:
                        # report and skip #7 emapa is obsolete
			msg = '%s (%s) emapaId is obsolete: %s (%s)' % (ubId, ubTerm, emapaId, emapaTerm)
			if msg not in obsAltEmapaInUberon:
			    obsAltEmapaInUberon.append(msg)
		    elif emapaId not in emapaDict:
			# report and skip #5  uberon that don't map to emapa - no uberon to emapa
			msg = '%s (%s) %s (%s) %s' % (mpId, mpTerm, ubId, ubTerm, emapaId)
			if msg not in mpNoEmapa:
			    mpNoEmapa.append(msg)
		    else:
			# load this mpId/emapaId relationship
			loadedCt += 1
			emapaRel = emapaDict[emapaId]

			objKey1 = mpTermKey
			objKey2 = emapaRel.termKey

			# add ids to sets for distinct count reporting later
			distinctMpLoaded.add(mpId)
			distinctEmapaLoaded.add(emapaId)

			# MGI_Relationship
			fpRelationshipFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % \
			    (nextRelationshipKey, TAB, catKey, TAB, objKey1, TAB, objKey2, TAB, relTermKey, TAB, qualKey, TAB, evidKey, TAB, refsKey, TAB, userKey, TAB, userKey, TAB, DATE, TAB, DATE, CRT))

			nextRelationshipKey += 1

    return 0

# end findRelationships() ---------------------------------------------------

def writeCuratorLog():
    # Purpose: writes QC errors to the curator log
    # Returns: 1 if error, else 0
    # Assumes: file descriptors have been initialized
    # Effects: sets global variables, writes to the file system
    # Throws: Nothing

    fpLogCur.write('\n%s Relationships Loaded\n\n' % loadedCt)
    fpLogCur.write('%s Distinct MP Loaded\n' % len(distinctMpLoaded))
    fpLogCur.write('%s Distinct EMAPA Loaded\n\n' % len(distinctEmapaLoaded))
    # #3
    if mpNotInDatabase:
	fpLogCur.write('MP Terms in the MP OWL file not in the database\n')
	fpLogCur.write('-' * 60 + '\n')
	fpLogCur.write(string.join(mpNotInDatabase, CRT))
	fpLogCur.write('\nTotal: %s' % len(mpNotInDatabase))
	fpLogCur.write('%s%s' % (CRT, CRT))
    # #4 - obsoleted
    if emapaNotInDatabase:
	fpLogCur.write('EMAPA Terms in the EMAPA OBO file not in the database\n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(emapaNotInDatabase, CRT))
        fpLogCur.write('%s%s' % (CRT, CRT))
    # #5  
    if mpNoEmapa:
	fpLogCur.write('MP Terms that do not map to EMAPA\n')
	fpLogCur.write('When MP and Uberon ID reported, Uberon ID not in Uberon file\n')
	fpLogCur.write('When MP, Uberon ID and EMAPA reported, EMAPA not in EMAPA file or no EMAPA cross reference\n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(mpNoEmapa, CRT))
	fpLogCur.write('\nTotal: %s' % len(mpNoEmapa))
        fpLogCur.write('%s%s' % (CRT, CRT))
    # #6
    if obsAltUberonInMP:
        fpLogCur.write('Obsolete or Alt Uberon Terms in the MP File\n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(obsAltUberonInMP, CRT))
	fpLogCur.write('\nTotal: %s' % len(obsAltUberonInMP))
        fpLogCur.write('%s%s' % (CRT, CRT))
    # #7
    if obsAltEmapaInUberon:
        fpLogCur.write('Obsolete or Alt EMAPA Terms in the Uberon File\n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(obsAltEmapaInUberon, CRT))
	fpLogCur.write('\nTotal: %s' % len(obsAltEmapaInUberon))
        fpLogCur.write('%s%s' % (CRT, CRT))
    # #8
    if oneMpMultiUberon:
        fpLogCur.write('MP Terms that Map to Multiple Uberon Terms\n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(oneMpMultiUberon, CRT))
	fpLogCur.write('\nTotal: %s' % len(oneMpMultiUberon))
        fpLogCur.write('%s%s' % (CRT, CRT))
    # #9
    if oneUberonMultiEmapa:
        fpLogCur.write('Uberon Terms that Map to Multiple EMAPA Terms\n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(oneUberonMultiEmapa, CRT))
	fpLogCur.write('\nTotal: %s' % len(oneUberonMultiEmapa))
        fpLogCur.write('%s%s' % (CRT, CRT))
    # #10 
    if someValuesFromLost:
	fpLogCur.write('EMAPA from "someValuesFrom" lost because using EMAPA from Description \n')
        fpLogCur.write('-' * 60 + '\n')
        fpLogCur.write(string.join(someValuesFromLost, CRT))
        fpLogCur.write('\nTotal: %s' % len(someValuesFromLost))
        fpLogCur.write('%s%s' % (CRT, CRT))

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

# parse the three input files
print '%s' % mgi_utils.date()
print 'running parseFiles'
if parseFiles() != 0:
    print 'Parsing Files failed'
    closeFiles()
    sys.exit(1)

# find the transitive relationships between mp and emapa
print '%s' % mgi_utils.date()
print 'running findRelationships'
if findRelationships() != 0:
    print 'Finding Relationships failed'
    closeFiles()
    sys.exit(1)

# write QC
print '%s' % mgi_utils.date()
print 'running writeCuratorLog'
if writeCuratorLog() != 0:
    print 'Writing the Curator Log failed'
    closeFiles()
    sys.exit(1)

if QC_ONLY == 'false':
    # delete existing relationships
    print '%s' % mgi_utils.date()
    print 'running doDeletes()'
    if doDeletes() != 0:
	print 'Do Deletes failed'
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
