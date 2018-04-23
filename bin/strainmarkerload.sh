#!/bin/sh

#
# This script is a wrapper around the process that loads 
# MGP and MGI B6 Strain Marker Objects 
#
#
#     strainmarker.sh
#

cd `dirname $0`/..
CONFIG_LOAD=`pwd`/strainmarkerload.config

cd `dirname $0`
LOG=`pwd`/strainmarkerload.log
rm -rf ${LOG}

USAGE='Usage: strainmarkerload.sh'

#
#  Verify the argument(s) to the shell script.
#
if [ $# -ne 0 ]
then
    echo ${USAGE} | tee -a ${LOG}
    exit 1
fi

#
# verify & source the configuration file
#

if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}"
    exit 1
fi

. ${CONFIG_LOAD}

#
# Just a verification of where we are at
#

echo "MGD_DBSERVER: ${MGD_DBSERVER}"
echo "MGD_DBNAME: ${MGD_DBNAME}"

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# create 'cleansed' input files and put them in INPUTDIR - check for minimum size
#
DOFILES=1
if [ $DOFILES ]
then
    echo "Preprocessing input files" | tee -a ${LOG}
    cd ${INPUT_GFF_DIR}
    for dir in ${INPUT_DIR_LIST}
    do
	echo ${dir}
	# parse strain name and add to top of file
	strain=`zcat ${dir}/*.gz | grep genome-version | cut -d' ' -f2 | cut -d_ -f1,2`
	echo ${strain}
	if [ -z "$strain" ]
	then
	    echo "Load Failed: ${INPUTDIR}/${dir}.txt has missing 'genome-version'" | tee -a ${LOG_DIAG} ${LOG_CUR}
	    exit 1
	fi
	echo ${strain} > ${INPUTDIR}/${dir}.txt
	# add only gene lines
	zcat ${dir}/*.gz | grep 'ID=gene:' >> ${INPUTDIR}/${dir}.txt
	count=`cat ${INPUTDIR}/${dir}.txt | wc -l`
	count=$((count-1)) # remove the strain name line from the count
	echo "count: $count"
	echo "min_records: ${MIN_RECORDS}"
	# check that for the minimum record count
	if [ ${count} -lt ${MIN_RECORDS} ]
	then
	    echo "Load Failed: ${INPUTDIR}/${dir}.txt has $count records which is less than the required ${MIN_RECORDS}" | tee -a ${LOG_DIAG} ${LOG_CUR}
	    exit 1
	fi
    done
fi

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
#

preload ${OUTPUTDIR}

#
# rm all files/dirs from OUTPUTDIR
#
#cleanDir ${OUTPUTDIR}

#
# run the load
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Run strainmarkerload.py"  | tee -a ${LOG_DIAG}
${STRAINMARKERLOAD}/bin/strainmarkerload.py
STAT=$?
checkStatus ${STAT} "${STRAINMARKERLOAD}/bin/strainmarkerload.py"

# run postload cleanup and email logs

shutDown

