#!/bin/bash
#SBATCH -N2
#SBATCH --exclusive

set -x

ulimit -s unlimited

NN=$SLURM_NNODES

if [ "x$NN" = "x" ]
then
  NN=1
else
  TMPDIR=$workdir/tmp/tmp.$SLURM_JOBID
  mkdir -p $TMPDIR
  cd $TMPDIR
fi

ls -l

REPO=$HOME/46t1_distatlas
RED=.20

data=$HOME
ln -sf $data/SFX_databases$RED/ecoclimapI_covers_param.bin ecoclimapI_covers_param.bin
ln -sf $data/SFX_databases$RED/orography.dir               SFX.ZS.dir
ln -sf $data/SFX_databases$RED/orography.hdr               SFX.ZS.hdr
ln -sf $data/SFX_databases$RED/ecoclimap.dir               SFX.COVER.dir 
ln -sf $data/SFX_databases$RED/ecoclimap.hdr               SFX.COVER.hdr 
ln -sf $data/SFX_databases$RED/SAND_HWSD_MOY.dir           SFX.SAND.dir
ln -sf $data/SFX_databases$RED/SAND_HWSD_MOY.hdr           SFX.SAND.hdr
ln -sf $data/SFX_databases$RED/CLAY_HWSD_MOY.dir           SFX.CLAY.dir
ln -sf $data/SFX_databases$RED/CLAY_HWSD_MOY.hdr           SFX.CLAY.hdr

rm PGD2*

if [ 0 -eq 1 ]
then

#xport ATLAS_TRACE=1
#xport ATLAS_TRACE_REPORT=1


~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 16 -openmp 8 -- $REPO/pgd
elif [ 0 -eq 1 ]
then
~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 2 -openmp 10 -- $REPO/pgd
else
#MP_NUM_THREADS=1 valgrind $REPO/pgd
#MP_NUM_THREADS=1 gdb --ex=run --args $REPO/pgd
OMP_NUM_THREADS=1 mpirun -np 4 $REPO/pgd
fi


