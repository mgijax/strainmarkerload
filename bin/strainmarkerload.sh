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
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" >> ${LOG_DIAG} 2>&1
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." >> ${LOG_DIAG} 2>&1
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
#date >> ${LOG_DIAG} 2>&1
#echo "Removing MGI GFF File from input directory"
#rm -rf "${INPUT_MGI_GFF_FILE}" >> ${LOG_DIAG} 2>&1
#echo "Copying new MGI GFF File from FTP site" >> ${LOG_DIAG} 2>&1
#echo "scp -p ${GFF3_SERVER}:${INPUT_MGI_GFF} ${INPUTDIR}" >> ${LOG_DIAG} 2>&1
#scp -p "${GFF3_SERVER}:${INPUT_MGI_GFF}" ${INPUTDIR} >> ${LOG_DIAG} 2>&1
#echo "Unzipping MGI GFF FILE" >> ${LOG_DIAG} 2>&1
#gunzip ${INPUT_MGI_GFF_FILE}.g>> ${LOG_DIAG} 2>&1

#
# run the load
#
date >> ${LOG_DIAG} 2>&1
echo "" >> ${LOG_DIAG}
echo "Run strainmarkerload.py"  >> ${LOG_DIAG} 2>&1
${PYTHON} -W ignore::SyntaxWarning ${STRAINMARKERLOAD}/bin/strainmarkerload.py >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "${STRAINMARKERLOAD}/bin/strainmarkerload.py" >> ${LOG_DIAG} 2>&1

#
# Archive a copy of the input file, adding a timestamp suffix.
#
#echo "" >> ${LOG_DIAG} 2>&1
#date >> ${LOG_DIAG} 2>&1
#echo "Archive input file" >> ${LOG_DIAG} 2>&1
#TIMESTAMP=`date '+%Y%m%d.%H%M'`
#ARC_FILE=`basename ${INPUT_MGI_GFF_FILE}`.${TIMESTAMP}
#cp -p ${INPUT_MGI_GFF_FILE} ${ARCHIVEDIR}/${ARC_FILE}

#
# run postload cleanup and email logs
#
shutDown
