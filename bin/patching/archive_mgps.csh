#!/bin/csh -f

#
# http://bhmgiwk01lp.jax.org/mediawiki/index.php?title=sw:Next_Scrum_Project#Ensembl_114_Strain_Gene_and_Gene_Model_update
#
# 1. _strain_key
# 2. _marker_key
# 3. strain
# 4. MGP ID
# 5. Ensembl ID (C57BL/6J)
# 6. Marker MGI ID
# 7. gene symbol
#


if ( ${?MGICONFIG} == 0 ) then
        setenv MGICONFIG /usr/local/mgi/live/mgiconfig
endif

source ${MGICONFIG}/master.config.csh

cd `dirname $0`

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date | tee -a $LOG
 
# https://useast.ensembl.org/Mus_musculus_C57BL_6NJ/Gene/Summary?db=core;g=MGP_C57BL6NJ_G0029802 

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG

--select * from map_coord_collection;

--select * from acc_mgitype where _mgitype_key in (19, 44);

select s._strain_key, m._marker_key, s.strain, a.accid, ens.accid, ma.accid, m.symbol
from acc_accession a, MRK_StrainMarker sm, PRB_Strain s, MRK_Marker m, ACC_Accession ens, ACC_Accession ma
where a._mgitype_key = 44
and a._object_key = sm._StrainMarker_key
and sm._strain_key = s._strain_key
and sm._marker_key = m._marker_key
and m._marker_key = ens._object_key
and ens._logicaldb_key = 60
and m._marker_key = ma._object_key
and ma._mgitype_key = 2
and ma._logicaldb_key = 1
and ma.preferred = 1
order by s.strain
;

EOSQL

date |tee -a $LOG

