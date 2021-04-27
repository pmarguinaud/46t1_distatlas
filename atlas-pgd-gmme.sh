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

PACK=/home/gmap/mrpm/marguina/46t1_distatlas

if [ 0 -eq 1 ]
then
  cp -f $PACK/data/fort.4.3x3ll fort.4
  RED=.10
elif [ 0 -eq 1 ]
then
  cp -f $PACK/data/fort.4.100x50ll fort.4
  RED=.10
elif [ 0 -eq 1 ]
then
  cp -f $PACK/data/fort.4.180x91 fort.4
  RED=.20
elif [ 0 -eq 1 ]
then
  cp -f $PACK/data/fort.4.64x64 fort.4
  RED=.20
elif [ 1 -eq 1 ]
then
  cp -f $PACK/data/fort.4.t149 fort.4
  RED=.20
else
  cp -f $PACK/data/fort.4.t1798c2.2 fort.4
# cp -f $PACK/data/fort.4.t1798     fort.4
# cp -f $PACK/data/fort.4.t4000     fort.4
# cp -f $PACK/data/fort.4.t8000     fort.4
  RED=""
fi

ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimapI_covers_param.bin ecoclimapI_covers_param.bin
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.dir               SFX.ZS.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.hdr               SFX.ZS.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.dir               SFX.COVER.dir 
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.hdr               SFX.COVER.hdr 
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.dir           SFX.SAND.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.hdr           SFX.SAND.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.dir           SFX.CLAY.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.hdr           SFX.CLAY.hdr


rm PGD2*

if [ 0 -eq 1 ]
then

#xport ATLAS_TRACE=1
#xport ATLAS_TRACE_REPORT=1


~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 16 -openmp 8 -- $PACK/pgd
elif [ 0 -eq 1 ]
then
~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 2 -openmp 10 -- $PACK/pgd
else
OMP_NUM_THREADS=1 /home/ext/bull/trivinoc/tools/valgrind/3.16.1/bin/valgrind $PACK/pgd
#MP_NUM_THREADS=1 gdb --ex=run --args $PACK/pgd
#MP_NUM_THREADS=1 $PACK/pgd
fi


