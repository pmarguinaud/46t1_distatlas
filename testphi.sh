#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --export="NONE"

set -x

JOBID=$(perl -e ' print(time ())')
PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x

ulimit -s unlimited
ulimit -l unlimited

export OMP_STACKSIZE=4G

set -e


NN=$SLURM_NNODES


if [ "x$NN" = "x" ]
then
  NN=1
  TMPDIR=/scratch/work/$USER/tmp/slurm-$$
else
  TMPDIR=/scratch/work/$USER/tmp/slurm-$JOBID
fi

mkdir -p $TMPDIR
cd $TMPDIR

NNP=4
let "NP=$NN*$NNP"

set +e
/opt/softs/mpiauto/mpiauto \
  --wrap --wrap-stdeo-pack --wrap-stdeo  \
  --prefix-mpirun "/usr/bin/time -f 'real=%e'" \
  --verbose -nn $NN -nnp $NNP -openmp 10 -- $PACK/bin/TESTPHI
set -e

ln -sf $PWD/stdeo.0 $PACK/stdeo.0.TESTPHI

for i in 0 1 2 3
do
  grep GREP stdeo.$i > $PACK/stdeo.GREP.TESTPHI.$i
done
