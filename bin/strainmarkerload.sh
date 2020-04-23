#!/bin/sh
#
#     strainmarker.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that loads 
#      MGP and MGI B6 Strain Genes 
#       
#      This script assumes the caller (straingenemodelload/bin/straingenemodelload.sh)
#
#  Usage:
#
#      strainmarker.sh
#
#  Env Vars:
#
#      See the configuration file: strainmarkerload.config
#
#  Outputs:
#
#       - Log file for the script initialization
#	- Diagnostic and curator logs
#	- MRK_StrainMarker bcp file
#	- ACC_Accession and ACC_AccessionReference bcp files
#	- input files for the straingenemodelload
#	    - b6 gene model file
#	    - b6 biotype file
#           - mgp gene model file
#           - mgp biotype file
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) unzip the MGI GFF3 file
#      2) if we are loading MGP (B6_ONLY=false),  preprocess the MGP files
#      3) run strainmarkerload.py
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
#  04/25/2018  sc  Initial development
#
###########################################################################

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
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
#

preload ${OUTPUTDIR} 

#
# Copy MGI.gff3 from public ftp site
#

if [ -f ${INPUT_MGI_GFF_FILE} ]
then
    echo "Removing MGI GFF File from input directory"
    rm  "${INPUT_MGI_GFF_FILE}" | tee -a ${LOG_DIAG}
fi

echo "Copying new MGI GFF File from FTP site" | tee -a ${LOG_DIAG}
echo "scp -p ${GFF3_SERVER}:${INPUT_MGI_GFF} ${INPUTDIR}"
scp -p "${GFF3_SERVER}:${INPUT_MGI_GFF}" ${INPUTDIR}

echo "Unzipping MGI GFF FILE" | tee -a ${LOG_DIAG}
gunzip ${INPUT_MGI_GFF_FILE}.gz

#
# create 'cleansed' MGP input files and put them in INPUTDIR - check for minimum size
#
if [ "${B6_ONLY}" = "false" ]
then
    echo "Preprocessing MGP input files in  ${INPUT_MGP_GFF_DIR}" | tee -a ${LOG_DIAG}
    cd ${INPUT_MGP_GFF_DIR}
    for dir in ${INPUT_MGP_DIR_LIST}
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
# run the load
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Run strainmarkerload.py"  | tee -a ${LOG_DIAG}
${PYTHON} ${STRAINMARKERLOAD}/bin/strainmarkerload.py
STAT=$?
checkStatus ${STAT} "${STRAINMARKERLOAD}/bin/strainmarkerload.py" >> ${LOG_DIAG} 2>&1

#
# Archive a copy of the input file, adding a timestamp suffix.
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Archive input file" >> ${LOG_DIAG}
TIMESTAMP=`date '+%Y%m%d.%H%M'`
ARC_FILE=`basename ${INPUT_MGI_GFF_FILE}`.${TIMESTAMP}
cp -p ${INPUT_MGI_GFF_FILE} ${ARCHIVEDIR}/${ARC_FILE}

#
# run postload cleanup and email logs
#
shutDown
