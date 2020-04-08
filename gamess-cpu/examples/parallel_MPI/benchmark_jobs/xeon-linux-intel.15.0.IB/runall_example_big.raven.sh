#!/bin/bash 
#PBS -l select=16:ncpus=16:mpiprocs=16
#PBS -l place=scatter:excl
#PBS -j oe
#PBS -N GAMESS-UK
#PBS -q workq
#PBS -l walltime=8:00:00
#PBS -P PR39

# latest intel compilers, mkl and intel-mpi

module purge
module load compiler/intel-15.0
module load mpi/intel.4.1.0

export OMP_NUM_THREADS=1

NCPUS=`cat ${PBS_NODEFILE} | wc -l`
PPN=`uniq -c ${PBS_NODEFILE} | awk '{print $1}' | uniq`

echo Number of Processing Elements is $NCPUS
echo Number of mpiprocs per node is $PPN

# Main GUK directory - change for your system
root=$HOME/MFG/GAMESS-UK

# The location of the binary to be validated
binary=$root/bin/gamess-uk

# Where the input files can be found
inputdir=$root/examples/parallel_MPI/input_files_benchmarks

# Where the output files can be found
outputdir=$root/examples/parallel_MPI/benchmark_jobs/xeon-linux-intel.15.0.IB

rundir=/scratch/$USER/GAMESS-UK.$PBS_JOBID

# The name of the file holding the list of jobs
jobs_list=$inputdir/jobs.big.list

# Loop over all the jobs listed in the jobs.list file and submit them:
# switch to where the jobs will be run
rm -rf ${rundir}
mkdir ${rundir}
cd ${rundir}

for job in `cat $jobs_list | grep -v "#"`
do
 echo running TEST=${job} NCPUs=$NCPUS PPN=$PPN
 cp $inputdir/${job}.in datain
 time mpirun -np $NCPUS ${binary} > ${outputdir}/${job}.n$NCPUS.PPN=$PPN.out.$PBS_JOBID
done
