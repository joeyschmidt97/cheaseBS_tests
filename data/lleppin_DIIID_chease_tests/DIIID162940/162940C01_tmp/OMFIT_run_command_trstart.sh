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
cat 162940C01CC.TMP
rm -f tr_start.log
eval `tset -s xterm`
chmod +x answerTRANSP.py
python3 -u answerTRANSP.py 162940C01 D3D pshare 8 1 1 0
