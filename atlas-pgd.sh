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


#cp -f $PACK/data/fort.4.t49  fort.4
 cp -f $PACK/data/fort.4.t149 fort.4
#cp -f $PACK/data/fort.4.t149c2.2 fort.4
#cp -f $PACK/data/fort.4.t32c2.4 fort.4


RED=.20

ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimapI_covers_param.bin ecoclimapI_covers_param.bin
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.dir               SFX.ZS.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/orography.hdr               SFX.ZS.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.dir               SFX.COVER.dir 
ln -sf /scratch/work/marguina/SFX_databases$RED/ecoclimap.hdr               SFX.COVER.hdr 
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.dir           SFX.SAND.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.hdr           SFX.SAND.hdr
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.dir           SFX.CLAY.dir
ln -sf /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.hdr           SFX.CLAY.hdr


rm PGD1.fa.* PGD2.fa.*

if [ 1 -eq 1 ]
then
~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 4 -openmp 10 -- $PACK/bin/ATLAS_PGD
elif [ 0 -eq 1 ]
then
~marguina/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 1 -openmp 10 -- $PACK/bin/ATLAS_PGD
else
OMP_NUM_THREADS=1 gdb --ex=run --args $PACK/bin/ATLAS_PGD
exit
fi

ls -lrt 

lfitools=$PACK/bin/lfitools

for p in PGD1 PGD2
do
$lfitools lfi_alt_index --lfi-file-in $p.fa.* --lfi-file-out $p.fa
$lfitools lfi_alt_pack --lfi-file-in $p.fa --lfi-file-out $p.pack.fa
done



