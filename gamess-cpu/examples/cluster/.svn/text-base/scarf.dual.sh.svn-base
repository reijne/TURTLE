#!/bin/bash
# this script is used to submit a job under the
# control of the bsub command, this
#    bsub -a mpich_gm -n 32 < scarf.dual.sh
# here running a 32 CPU job with 2 CPUs per node
#BSUB -n 32
#BSUB -W 60
#BSUB -o par_44.out.32.16X2
#BSUB -J GAMESS
#BSUB -R "span[ptile=2]"

#export GMPI_SHMEM=0

export LD_LIBRARY_PATH=/opt/gm/lib:${LD_LIBRARY_PATH}

env

cd  /home/cse/cse0002/release-7-0/GAMESS-UK/examples/cluster

cp par_44.in datain

mpirun.lsf_pgi -np 32 /home/cse/cse0002/release-7-0/GAMESS-UK/bin/gamess-uk

