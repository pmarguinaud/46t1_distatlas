#!/bin/bash

set -x

ulimit -s unlimited
export DR_HOOK=0
export DR_HOOK_NOT_MPI=1


PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x


\rm -f UV1*

#~/SAVE/mpiauto/mpiauto -verbose -np 1 --wrap --wrap-stdeo -- $PACK/bin/ATLAS_ARPEGE_F --grid1 N16 --rotate  --grid2 data/fort.4.t32c2.4
 ~/SAVE/mpiauto/mpiauto -verbose -np 1 --wrap --wrap-stdeo -- $PACK/bin/ATLAS_ARPEGE_F --grid2 N16 --rotate  --grid1 data/fort.4.t32c2.4


ls UV1* > /dev/null 2>&1

if [ $? -eq 0 ]
then
$PACK/bin/lfitools lfi_alt_index --lfi-file-in UV1.fa.* --lfi-file-out UV1.fa
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in UV1.fa --lfi-file-out UV1.pack.fa
$PACK/bin/lfitools lfi_alt_index --lfi-file-in UV1R.fa.* --lfi-file-out UV1R.fa
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in UV1R.fa --lfi-file-out UV1R.pack.fa
fi
