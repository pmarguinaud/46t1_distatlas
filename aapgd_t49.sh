#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --export="NONE"

set -x

JOBID=$(perl -e ' print(time ())')
PACK=/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x

#exec > aatestprog-$JOBID.eo 2>&1


ulimit -s unlimited
ulimit -l unlimited

unset DR_HOOK_NOT_MPI
export DR_HOOK=1
export DR_HOOK_OPT=prof
export DR_HOOK_SILENT=1
export EC_MEMINFO=0
export DR_HOOK_SHOW_PROCESS_OPTIONS=0
export DR_HOOK_IGNORE_SIGNALS=-1
export MPL_MBX_SIZE=2000000000
export OMP_STACKSIZE=4G

set -e


NN=$SLURM_NNODES


RED=".20"

if [ "x$NN" = "x" ]
then
  NN=1
  TMPDIR=/scratch/work/$USER/tmp/slurm-$$
else
  TMPDIR=/scratch/work/$USER/tmp/slurm-$JOBID
fi

mkdir -p $TMPDIR
cd $TMPDIR

#for TR in t1198 t1798 t2048 t4000 t6000 t8000
for TR in t49
#for TR in t149qf
do

mkdir -p $TR
cd $TR

ln -s /scratch/work/marguina/SFX_databases$RED/ecoclimapI_covers_param.bin .
ln -s /scratch/work/marguina/SFX_databases$RED/orography.dir relief.dir
ln -s /scratch/work/marguina/SFX_databases$RED/orography.hdr relief.hdr
#n -s /scratch/work/marguina/SFX_databases$RED/gtopo30.dir relief.dir
#n -s /scratch/work/marguina/SFX_databases$RED/gtopo30.hdr relief.hdr
ln -s /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.dir SAND.dir
ln -s /scratch/work/marguina/SFX_databases$RED/SAND_HWSD_MOY.hdr SAND.hdr
ln -s /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.dir CLAY.dir
ln -s /scratch/work/marguina/SFX_databases$RED/CLAY_HWSD_MOY.hdr CLAY.hdr
ln -s /scratch/work/marguina/SFX_databases$RED/ecoclimap.dir ecoclimap.dir 
ln -s /scratch/work/marguina/SFX_databases$RED/ecoclimap.hdr ecoclimap.hdr 

\rm -f drhook.prof.* *.msh stdeo.*

NNP=4
let "NP=$NN*$NNP"



cat -> fort.4 << EOF
&NAMPGD
  LGARDEN=.FALSE.,
  LWATER_TO_NATURE=.FALSE.,
  LTOWN_TO_ROCK=.TRUE.,
/
&NAMGRID
$(cat $PACK/fort.4.$TR)
/
EOF

#xport MPIAUTOCONFIG=mpiauto.DDT.conf
/opt/softs/mpiauto/mpiauto \
  --wrap --wrap-stdeo-pack --wrap-stdeo  \
  --prefix-mpirun "/usr/bin/time -f 'real=%e'" \
  --verbose -nn $NN -nnp $NNP -openmp 10 -- $PACK/bin/AAPGD

ln -sf $PWD/stdeo.0 $PACK/stdeo.0

for i in 0 1 2 3
do
  grep GREP stdeo.$i > $PACK/stdeo.GREP.AAPGD.$i
done


(
  export OMP_NUM_THREADS=1
  export DR_HOOK_NOT_MPI=1
  export DR_HOOK=0

  for b in ZXYZ2I ZXYZ2 ZXYZ1
  do
  $PACK/bin/lfitools lfi_alt_index --lfi-file-in $b.fa.*  --lfi-file-out $b.fa > index.eo 2>&1
  $PACK/bin/lfitools lfi_alt_pack --lfi-file-in $b.fa --lfi-file-out $b.pack.fa
  done

)

# scp *.pack.fa phil@s-entraider.net:


cd $TMPDIR

done




