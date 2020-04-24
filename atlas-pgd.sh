#!/bin/bash
#SBATCH -N4
#SBATCH --exclusive

set -x

ulimit -s unlimited
export DR_HOOK_NOT_MPI=1 
export DR_HOOK=0 

NN=$SLURM_NNODES

if [ "x$NN" = "x" ]
then
  NN=1
else
  TMPDIR=$workdir/tmp/tmp.$SLURM_JOBID
  mkdir -p $TMPDIR
  cd $TMPDIR
fi

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x

ls -l


cp -f $PACK/data/fort.4.t49 fort.4


RED=.20

ln -sf /scratch/work/marguina/SFX_databases$RED/orography.dir relief.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.hdr relief.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.dir ecoclimap.dir 
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.hdr ecoclimap.hdr 
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.dir SAND.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.hdr SAND.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.dir CLAY.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.hdr CLAY.hdr


rm PGD.fa.*

~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 4 -openmp 10 -- $PACK/bin/ATLAS_PGD

ls -lrt 

exit

lfitools=$PACK/bin/lfitools

$lfitools lfi_alt_index --lfi-file-in HALFDIFF1.fa* --lfi-file-out HALFDIFF1.fa
$lfitools lfi_alt_pack --lfi-file-in HALFDIFF1.fa --lfi-file-out HALFDIFF1.pack.fa


