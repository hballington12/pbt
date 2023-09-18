#!/bin/sh
#PBS -l walltime=0:10:00
#PBS -l nodes=2
#PBS -k oe
#PBS -o /Multirot/test
#PBS -e /Multirot/test
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

eval `/usr/bin/modulecmd bash load gcc-13.1`
eval `/usr/bin/modulecmd bash load openmpi-4.0.5`

export OMP_NUM_THREADS=32
export OMP_STACKSIZE=100M
ulimit -s unlimited

cd /home/hballington/abt

mpiexec \
/home/hballington/abt/src/mpi/abt \
-lambda 0.532 \
-rbi 1.3117 \
-ibi 0 \
-rec 8 \
-cmethod cc_hex \
-cft obj \
-cfn big_hex.obj \
-rot multi 40 \
-jobname test \
-cc_hex_l 20 \
-cc_hex_hr 10 \
-cc_hex_nfhr 4 \
-cc_hex_pfl 10 \
-cc_hex_nfpl 4 \
-cc_hex_pher 1 \
-cc_hex_pper 2 \
-cc_hex_nscales 1 \
-cc_hex_cls 1 \
-cc_hex_sds 0 \
-theta 0 0.1 1 0.5 5 2 179 0.25 180 \
-phi 0 2 360 \
-mt \
> log_mpi