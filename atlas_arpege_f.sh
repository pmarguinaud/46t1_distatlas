#!/bin/bash

set -x

ulimit -s unlimited
export DR_HOOK_NOT_MPI=1 
export DR_HOOK=0 

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


 ~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 1 -- $PACK/bin/ATLAS_ARPEGE_F --grid1 X80x40 --dist1 checkerboard --grid2 fort.4.t32 --dist2 equal_regions



ls XYZ1.fa.* > /dev/null 2>&1

if [ $? -eq 0 ]
then
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ1.fa.*  --lfi-file-out XYZ1.fa > index.eo 2>&1
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ1.fa --lfi-file-out XYZ1.pack.fa
fi

$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ2I.fa.*  --lfi-file-out XYZ2I.fa > index.eo 2>&1
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ2I.fa --lfi-file-out XYZ2I.pack.fa



