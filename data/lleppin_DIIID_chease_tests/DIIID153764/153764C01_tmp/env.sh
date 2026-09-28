#!/bin/bash -l
module purge
module load ntcc
export PREACTDIR=/dev/null
export ADASDIR=/dev/null
export TRANSPGRID_SERVER=transpgrid.pppl.gov
export MDS_TRANSP_SERVER='atlas.gat.com'
export MDS_TRANSP_TREE='transp'
export TR_EMAIL='aon2113@columbia.edu'
echo $MDS_TRANSP_SERVER : $MDS_TRANSP_TREE
cat 153764C01CC.TMP
rm -f tr_start.log
eval `tset -s xterm`
