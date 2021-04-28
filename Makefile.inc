
ATLAS_PREFIX=$(HOME)/atlas/install

F90=gfortran -DADDRESS64 -fconvert=swap -fopenmp -g -O2 -I$(ATLAS_PREFIX)/module/atlas -I$(ATLAS_PREFIX)/module/fckit -I$(ATLAS_PREFIX)/include 
CXX=g++ -fopenmp -g -O2 -std=c++11 -I$(ATLAS_PREFIX)/include
CC=gcc -fopenmp -g -O2

LIBS=-Wl,-rpath,$(ATLAS_PREFIX)/lib64 -L$(ATLAS_PREFIX)/lib64 -latlas_f -latlas -lfckit -leckit_mpi -leckit -lstdc++
LD=mpif90 -fopenmp -g 


