#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH -p huge512,normal256,nmipt
#SBATCH --exclusive

set -x

ulimit -s unlimited
export DR_HOOK_NOT_MPI=1 
export DR_HOOK=0 

NN=$SLURM_NNODES

if [ "x$NN" = "x" ]
then
  NN=1
  TMPDIR=$workdir/tmp/tmp.$$
  TMPDIR=$workdir/tmp/tmp.$SLURM_JOBID
  mkdir -p $TMPDIR
else
  TMPDIR=$workdir/tmp/tmp.$SLURM_JOBID
  mkdir -p $TMPDIR
  cd $TMPDIR
fi

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlasLATEST.01.I185274INTELMPI185274.x

ls -l

cp -f $PACK/data/fort.4.t498 fort.4
RED=.20

ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimapI_covers_param.bin ecoclimapI_covers_param.bin
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.dir               SURFGEOPOTENTIEL.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.hdr               SURFGEOPOTENTIEL.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.dir               SFX.COVER.dir 
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.hdr               SFX.COVER.hdr 
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.dir           SURFPROP.SABLE.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.hdr           SURFPROP.SABLE.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.dir           SURFPROP.ARGILE.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.hdr           SURFPROP.ARGILE.hdr

ln -sf /scratch/work/marguina/clim/dir/*_GL .
ln -sf /scratch/work/marguina/clim/dir/*_GL.* .
ln -sf /scratch/work/marguina/clim/dir/abc_coef.* .



rm Const.Clim*

if [ 1 -eq 1 ]
then
OMP_NUM_THREADS=1 $PACK/bin/ATLAS_CLIM > ATLAS_CLIM.eo 2>&1
else

/opt/softs/intel/2018.04/impi/2018.5.274/intel64/bin/mpirun -np 2 -ppn 2 \
/home/gmap/mrpm/marguina/pack/46t1_distatlasLATEST.01.I185274INTELMPI185274.x/bin/ATLAS_CLIM  > ATLAS_CLIM 2>&1

fi

ls -lrt 

lfitools=$PACK/bin/lfitools

for p in Const.Clim Const.Clim.m01 Const.Clim.m02 Const.Clim.m03 Const.Clim.m04 Const.Clim.m05 Const.Clim.m06 \
  Const.Clim.m07 Const.Clim.m08 Const.Clim.m09 Const.Clim.m10 Const.Clim.m11 Const.Clim.m12
do
$lfitools lfi_alt_index --lfi-file-in $p.* --lfi-file-out $p
$lfitools lfi_alt_pack --lfi-file-in $p --lfi-file-out $p.pack
done


pwd


