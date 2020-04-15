#!/bin/bash

set -x

\rm testphi.x

/home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/INTELMPI184274MT/mpic++ /home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/I185274/icpc -std=c++11 \
  -march=core-avx2 -no-fma -fp-model source -fp-model precise -fimf-use-svml=true -g -traceback -fPIC -ip -ftz -malign-double -qopenmp \
  -o testphi.x src/local/atlas/programs/testphi.cc -I src/local/atlas/include -limf -lintlc -lsvml \
  -I /home/gmap/mrpm/marguina/atlas-intel18/install/include \
  -L/home/gmap/mrpm/marguina/atlas-intel18/install/lib -latlas -leckit -leckit_mpi -Wl,-rpath,/home/gmap/mrpm/marguina/atlas-intel18/install/lib

./testphi.x 
