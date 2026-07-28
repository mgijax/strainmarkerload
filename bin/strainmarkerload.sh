#!/bin/sh
#
#     strainmarker.sh
###########################################################################
#
#  Purpose:
#
#      This script is called from straingenemodelload/bin/straingenemodelload.sh.
#      This script is a wrapper around the process that loads MGP and MGI B6 Strain Genes 
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
#  Implementation:
#
#      This script will perform following steps:
#
#      1) run ensembl gff3 116 patching : 'patched' folder
#      2) copy patched files from 'patched' folder to 'input' folder and unzip
#      3) copy & unzip the MGI gff3 file from '/export/???/ftp/pub/mgigff3' folder to 'input' folder
#      4) run strainmarkerload.py
#
###########################################################################
#
#  Modification History:
#
#  Date        SE   Change Description
#  ----------  ---  -------------------------------------------------------
#  07/28/2026  lec  sprt-26/Strainmarker load for new Strain Genes
#  07/28/2026  jer  sprt-383/Strain Gene GFF3 file patching
#  04/25/2018  sc   Initial development
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
# Run the Ensembl Gff3 Patching
#
#date >> ${LOG_DIAG} 2>&1
#echo "Running patchEnsemblGff116.sh" >> ${LOG_DIAG} 2>&1
#pushd ${STRAINMARKERLOAD}/bin/patching >> ${LOG_DIAG} 2>&1
#./patchEnsemblGff116.sh ${INPUT_MGP_GFF_DIR}/*/*.gff3.gz >> ${LOG_DIAG} 2>&1
#popd >> ${LOG_DIAG} 2>&1

#
# Copy the Ensembl Gff3 Patched Files to Input Folder & Unzip
#
#date >> ${LOG_DIAG} 2>&1
#echo "Removing old Strain GFF3 Files from input directory" >> ${LOG_DIAG} 2>&1
#rm -rf ${INPUTDIR}/Mus*.gff3 >> ${LOG_DIAG} 2>&1
#echo "Copying new Strain GFF3 Files from patched directory & gzip" >> ${LOG_DIAG} 2>&1
#cp ${PATCH_ODIR}/*/*.gz ${INPUTDIR} >> ${LOG_DIAG} 2>&1
#cd ${INPUTDIR}
#for i in *.gz
#do
#gunzip $i
#done

#
# Copy MGI.gff3 from public ftp site
#
#date >> ${LOG_DIAG} 2>&1
#echo "Removing MGI GFF File from input directory" >> ${LOG_DIAG} 2>&1
#rm -rf ${INPUT_MGI_GFF_FILE} >> ${LOG_DIAG} 2>&1
#rm -rf ${INPUT_MGI_GFF_FILE}.gz >> ${LOG_DIAG} 2>&1
#echo "Copying new MGI GFF File from FTP site" >> ${LOG_DIAG} 2>&1
#cp ${INPUT_MGI_GFF} ${INPUTDIR} >> ${LOG_DIAG} 2>&1
#echo "Unzipping MGI GFF Files" >> ${LOG_DIAG} 2>&1
#gunzip ${INPUT_MGI_GFF_FILE}.gz >> ${LOG_DIAG} 2>&1

#
# run the load
#
date >> ${LOG_DIAG} 2>&1
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
