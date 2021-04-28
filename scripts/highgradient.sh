#!/bin/bash
#SBATCH -N1
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

rm -rf XYZ1* stdeo.* 

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI185274.x

ls -l

~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 1 -openmp 1 -- $PACK/bin/ATLAS_HIGHGRADIENT \
  --grid $PACK/data/fort.4.64x64 --xy 2 --order 3

ls -lrt 

lfitools=$PACK/bin/lfitools

for f in XYZ1 XYZ1GR
do
$lfitools lfi_alt_index --lfi-file-in ${f}.fa* --lfi-file-out ${f}.fa
$lfitools lfi_alt_pack --lfi-file-in ${f}.fa --lfi-file-out ${f}.pack.fa
done

ls -lrt $PWD/*.fa


