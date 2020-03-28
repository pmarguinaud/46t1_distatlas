#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --export="NONE"

set -x


rm bin/TESTFLO stdeo.0.TESTFLO

set -e
CXX="/home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/INTELMPI184274MT/mpic++ /home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/I185274/icpc -std=c++11"

$CXX -g -c -I/home/gmap/mrpm/marguina/atlas-intel/install/include -Isrc/local/atlas/include -o src/local/atlas/interpolation/interpolationA.o src/local/atlas/interpolation/interpolationA.cc

$CXX -g -o bin/TESTFLO src/local/atlas/programs/testflo.cc \
  src/local/atlas/interpolation/interpolationA.o \
  -I/home/gmap/mrpm/marguina/atlas-intel/install/include \
  -Isrc/local/atlas/include \
  -L/home/gmap/mrpm/marguina/atlas-intel/install/lib -Wl,-rpath,/home/gmap/mrpm/marguina/atlas-intel/install/lib \
  -latlas -leckit -leckit_geometry -leckit_mpi
set +e



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
  --verbose -nn $NN -nnp $NNP -openmp 1 -- $PACK/bin/TESTFLO
set -e

ln -sf $PWD/stdeo.0 $PACK/stdeo.0.TESTFLO

for i in 0 1 2 3
do
  grep GREP stdeo.$i > $PACK/stdeo.GREP.TESTFLO.$i
done