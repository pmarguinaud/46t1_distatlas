#!/bin/bash

set -x

\rm testflo.x

/home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/INTELMPI184274MT/mpic++ /home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/I185274/icpc -std=c++11 \
  -march=core-avx2 -no-fma -fp-model source -fp-model precise -fimf-use-svml=true -g -traceback -fPIC -ip -ftz -malign-double -qopenmp \
  -o testflo.x src/local/atlas/programs/testflo.cc -I /home/gmap/mrpm/marguina/atlas-intel/install/include -limf -lintlc -lsvml \
  -L/home/gmap/mrpm/marguina/atlas-intel/install/lib -latlas -leckit -leckit_mpi -Wl,-rpath,/home/gmap/mrpm/marguina/atlas-intel/install/lib


~/SAVE/mpiauto/mpiauto -verbose -np 4 --wrap --wrap-stdeo -- ./testflo.x
