#!/bin/bash
#SBATCH -N4
#SBATCH --exclusive

set -x

ulimit -s unlimited
export DR_HOOK_NOT_MPI=1 
export DR_HOOK=0 

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x

cd $TMPDIR

$PACK/atlas_arpege_f-interp.pl --small --large --trace --interp4 --interpA

