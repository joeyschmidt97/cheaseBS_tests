#!/bin/bash -l
module purge
module load ntcc

export TRANSPGRID_SERVER=transpgrid.pppl.gov
tr_send 153764C01
