#!/bin/bash

set -x

ulimit -s unlimited
export DR_HOOK_NOT_MPI=1 
export DR_HOOK=0 

rm XYZ1.fa*

PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x

~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- $PACK/bin/ATLAS_ARPEGE_F 

$PACK/bin/lfitools lfi_alt_index --lfi-file-in XYZ1.fa.*  --lfi-file-out XYZ1.fa > index.eo 2>&1
$PACK/bin/lfitools lfi_alt_pack --lfi-file-in XYZ1.fa --lfi-file-out XYZ1.pack.fa



