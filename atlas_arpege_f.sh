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

rm -rf XYZ1* stdeo.* XYZ2I*

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x

#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 1 -- $PACK/bin/ATLAS_ARPEGE_F --grid2 fort.4.t32
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid2 t31/ICMSHARPEINIT
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid2 fort.4.64x64  --dist2 checkerboard
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid2 ICMSHAROMINIT --dist2 checkerboard
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid2 ICMSHAROMINIT --dist2 checkerboard
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid2 fort.4.t32 --dist2 equal_regions
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 S80x40 --dist1 checkerboard --grid2 fort.4.t32 --dist2 equal_regions
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 fort.4.t32 --dist1 equal_regions --grid2 fort.4.64x64 --dist2 checkerboard
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 fort.4.64x64 --dist1 checkerboard --grid2 fort.4.32x32 --dist2 checkerboard --write1
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 L80x40 --dist1 checkerboard --grid2 fort.4.32x32 --dist2 checkerboard 
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 fort.4.3x3 --dist1 checkerboard --grid2 fort.4.32x32 --dist2 checkerboard
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 1 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 fort.4.t32 --dist1 equal_regions --grid2 fort.4.32x32 --dist2 checkerboard
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 X80x40 --dist1 checkerboard --grid2 fort.4.t32 --dist2 equal_regions

#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 L160x80 --dist1 checkerboard --grid2 fort.4.t32 --dist2 equal_regions --interpA
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 L40x20 --dist1 checkerboard --grid2 fort.4.t32 --dist2 equal_regions --interpA


#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 2 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 L40x20 --dist1 checkerboard --grid2 fort.4.64x64 --dist2 checkerboard --interpA
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 N320 --dist1 equal_regions --grid2 fort.4.64x64 --dist2 checkerboard --interpA --write1
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 fort.4.64x64 --dist1 checkerboard --grid2 fort.4.32x32_100km --dist2 checkerboard --interpA --write1
#~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 1 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 fort.4.64x64 --dist1 checkerboard --grid2 fort.4.32x32_100km --dist2 checkerboard --interpA --write1

 ~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F \
     --grid1 L400x200 --grid2 PFFCSTEUROC25+0000 --interp4 --write2

if [ 0 -eq 1 ]
then

export ATLAS_TRACE=1
export ATLAS_TRACE_REPORT=1

# --prefix-command '/home/gmap/mrpm/marguina/bin/perf_wrap --call-graph fp' \
~/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 8 -openmp 16 -- $PACK/bin/ATLAS_ARPEGE_F \
  --grid1 N1024        --dist1 equal_regions --interp4 \
  --grid2 L40000x20000 --dist2 checkerboard  --block2 40000 --light2

elif [ 0 -eq 1 ]
then

export ATLAS_TRACE=1
export ATLAS_TRACE_REPORT=1

# --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
~/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/home/gmap/mrpm/marguina/bin/perf_wrap --call-graph fp' \
  --wrap --wrap-stdeo -nn $NN -nnp 4 -openmp 32 -- $PACK/bin/ATLAS_ARPEGE_F \
  --grid1 L40000x20000 --dist1 checkerboard  --block1 40000 --light1 \
  --grid2 N1024        --dist2 equal_regions --interpA 

elif [ 0 -eq 1 ]
then

~/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 4 -openmp 32 -- $PACK/bin/ATLAS_ARPEGE_F \
  --grid1 L4000x2000 --dist1 checkerboard  --block1 4000 --light1 \
  --grid2 N512       --dist2 equal_regions --interpA 

elif [ 0 -eq 1 ]
then

~/SAVE/mpiauto/mpiauto \
  --prefix-mpirun '/usr/bin/time -f "time=%es"' \
  --prefix-command '/usr/bin/time -f "mem=%Mkb"' \
  --wrap --wrap-stdeo -nn $NN -nnp 1 -openmp 1 -- $PACK/bin/ATLAS_ARPEGE_F \
  --grid1 L400x200 --dist1 checkerboard  --block1 400 --light1 \
  --grid2 N16      --dist2 equal_regions --interpA 

fi

ls XYZ1.fa.* > /dev/null 2>&1

if [ $? -eq 0 ]
then
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ1.fa.*  --lfi-file-out XYZ1.fa > index.eo 2>&1
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ1.fa --lfi-file-out XYZ1.pack.fa
fi

ls XYZ2I4* > /dev/null 2>&1

if [ $? -eq 0 ]
then
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ2I4.fa.*  --lfi-file-out XYZ2I4.fa > index.eo 2>&1
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ2I4.fa --lfi-file-out XYZ2I4.pack.fa
fi


ls XYZ2IA* > /dev/null 2>&1

if [ $? -eq 0 ]
then
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ2IA.fa.*  --lfi-file-out XYZ2IA.fa > index.eo 2>&1
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ2IA.fa --lfi-file-out XYZ2IA.pack.fa
fi


if [ "x$SLURM_NNODES" != "x" ]
then

ls -lrt $PWD/*.pack.fa

\rm -f $PACK/XYZ*
for f in $PWD/*.pack.fa
do
  ln -sf $f $PACK/
done

fi


