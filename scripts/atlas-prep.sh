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

ln -sf /scratch/work/marguina/t500/pgd_arpege.t1798.01.fa    PGD1
ln -sf /scratch/work/marguina/t500/ICMSH1798INIT.sfx         PRE1
ln -sf /scratch/work/marguina/PGD.t149.fa                    PGD2

rm PRE2*

~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 4 -openmp 10 -- $PACK/bin/ATLAS_PREP

ls -lrt 

lfitools=$PACK/bin/lfitools

for p in PRE2
do
$lfitools lfi_alt_index --lfi-file-in $p.* --lfi-file-out $p
$lfitools lfi_alt_pack --lfi-file-in $p --lfi-file-out $p.pack
done



