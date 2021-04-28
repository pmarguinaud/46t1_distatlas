#!/bin/bash
#SBATCH -N1
#SBATCH --exclusive

set -x

ulimit -s unlimited
export DR_HOOK_NOT_MPI=1 
export DR_HOOK=0 
export OMP_STACK_SIZE=2G

NN=$SLURM_NNODES

if [ "x$NN" = "x" ]
then
  NN=1
else
  TMPDIR=$workdir/tmp/tmp.$SLURM_JOBID
  mkdir -p $TMPDIR
  cd $TMPDIR
fi

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI185274.x

ls -l


ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimapI_covers_param.bin ecoclimapI_covers_param.bin
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.dir               SFX.ZS.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.hdr               SFX.ZS.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.dir               SFX.COVER.dir 
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.hdr               SFX.COVER.hdr 
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.dir           SFX.SAND.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.hdr           SFX.SAND.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.dir           SFX.CLAY.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.hdr           SFX.CLAY.hdr


rm PGD2*

export ECKIT_MPI_INIT_THREAD=MPI_THREAD_MULTIPLE


if [ 0 -eq 1 ]
then

#xport ATLAS_TRACE=1
#xport ATLAS_TRACE_REPORT=1


~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 16 -openmp 8 -- $PACK/bin/ATLAS_PGD_POINTS
elif [ 1 -eq 1 ]
then

#export MPIAUTOCONFIG=mpiauto.XGDB1.conf

~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 2 -openmp 10 -- $PACK/bin/ATLAS_PGD_POINTS
else

OMP_NUM_THREADS=10 gdb --ex=run --args $PACK/bin/ATLAS_PGD_POINTS

exit

rm log.*.txt.? log.*.txt

OMP_NUM_THREADS=1 $PACK/bin/ATLAS_PGD_POINTS > stdeo.0.1 2>&1

for l in log.*.txt
do
 mv $l $l.1
done

OMP_NUM_THREADS=2 $PACK/bin/ATLAS_PGD_POINTS > stdeo.0.2 2>&1

for l in log.*.txt
do
 mv $l $l.2
done

#MP_NUM_THREADS=2 gdb --ex=run --args $PACK/bin/ATLAS_PGD_POINTS
#MP_NUM_THREADS=1 $PACK/bin/ATLAS_PGD_POINTS
fi


ls -lrt 

