#!/bin/bash
#SBATCH -N1
#SBATCH --time=00:02:00

set -x
ulimit -s unlimited

cd /home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI185274.x

~/SAVE/mpiauto/mpiauto -np 4 --verbose --wrap --wrap-stdeo -- ./bin/TESTPHI


