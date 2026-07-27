#!/bin/bash

if [ -e "../../strainmarkerload.config" ]
then
    source ../../strainmarkerload.config
else
    echo "Cannot find configuration file."
    exit 1
fi


function log {
    echo $* >> ${PATCH_LOG}
}

rm -f ${PATCH_LOG}
log "PatchEnsemblGff116 started at:" `date`

while [ $# -gt 0 ]; do
    dir="${1%/*}"
    strain="${dir##*/}"
    fname="${1##*/}"
    lastExt="${fname##*.}"
    fnameNoExt="${fname%.*}"
    ofile="${PATCH_ODIR}/${strain}/${fnameNoExt}"
    if [ "${lastExt}" == "gz" ] ; then
	ofile="${ofile}.gz"
    fi
    log ""
    log `date`
    log "Current argument: $1"
    log "dir=" $dir
    log "strain=" $strain
    log "fname=" $fname
    log "lastExt=" $lastExt
    log "fnameNoExt=" $fnameNoExt
    log "ofile=" $ofile
    mkdir -p ${PATCH_ODIR}/${strain}
    if [ "${lastExt}" == "gz" ] ; then
	gunzip -c $1 | $PYTHON patchEnsemblGff116.py -L $PATCH_PPG_LIMIT 2>> ${PATCH_LOG} | gzip > ${ofile} 
    else
	$PYTHON patchEnsemblGff116.py -L $PATCH_PPG_LIMIT < "$1" > "${ofile}" 2>> ${PATCH_LOG}
    fi
    shift
done

log "PatchEnsemblGff116 ended at:" `date`
