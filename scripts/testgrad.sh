#!/bin/bash

set -x

ulimit -s unlimited
export DR_HOOK=0
export DR_HOOK_NOT_MPI=1


PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x


\rm -f XYZ1*

#~/SAVE/mpiauto/mpiauto -verbose -np 4 --wrap --wrap-stdeo -- $PACK/bin/ATLAS_ARPEGE_F --grid1 N16 --gradient --write1
#~/SAVE/mpiauto/mpiauto -verbose -np 4 --wrap --wrap-stdeo -- $PACK/bin/ATLAS_ARPEGE_F --grid1 L80x40 --gradient --write1
#~/SAVE/mpiauto/mpiauto -verbose -np 4 --wrap --wrap-stdeo -- $PACK/bin/ATLAS_ARPEGE_F --grid1 data/fort.4.t32c2.4 --gradient --write1 
 ~/SAVE/mpiauto/mpiauto -verbose -np 4 --wrap --wrap-stdeo -- $PACK/bin/ATLAS_ARPEGE_F --grid1 data/fort.4.64x64 --gradient --write1


ls XYZ1* > /dev/null 2>&1

if [ $? -eq 0 ]
then
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ1.fa.* --lfi-file-out XYZ1.fa
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ1.fa --lfi-file-out XYZ1.pack.fa
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ1GT.fa.* --lfi-file-out XYZ1GT.fa
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ1GT.fa --lfi-file-out XYZ1GT.pack.fa
$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ1GC.fa.* --lfi-file-out XYZ1GC.fa
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ1GC.fa --lfi-file-out XYZ1GC.pack.fa
fi
