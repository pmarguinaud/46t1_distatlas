#!/bin/bash

set -x

ulimit -s unlimited

ls -l

REPO=$HOME/46t1_distatlas
RED=.20

data=$HOME

if [ ! -f "$data/SFX_databases$RED/ecoclimapI_covers_param.bin" ]
then
  echo "Cannot find $data/SFX_databases$RED/ecoclimapI_covers_param.bin"
  exit 1
fi

ln -sf $data/SFX_databases$RED/ecoclimapI_covers_param.bin ecoclimapI_covers_param.bin
ln -sf $data/SFX_databases$RED/orography.dir               SFX.ZS.dir
ln -sf $data/SFX_databases$RED/orography.hdr               SFX.ZS.hdr
ln -sf $data/SFX_databases$RED/ecoclimap.dir               SFX.COVER.dir 
ln -sf $data/SFX_databases$RED/ecoclimap.hdr               SFX.COVER.hdr 
ln -sf $data/SFX_databases$RED/SAND_HWSD_MOY.dir           SFX.SAND.dir
ln -sf $data/SFX_databases$RED/SAND_HWSD_MOY.hdr           SFX.SAND.hdr
ln -sf $data/SFX_databases$RED/CLAY_HWSD_MOY.dir           SFX.CLAY.dir
ln -sf $data/SFX_databases$RED/CLAY_HWSD_MOY.hdr           SFX.CLAY.hdr

rm -f PGD2*

OMP_NUM_THREADS=1 mpirun -np 4 $REPO/pgd


