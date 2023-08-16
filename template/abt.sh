#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
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



export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64_lin
export OMP_NUM_THREADS=32
export OMP_STACKSIZE=100M
ulimit -s unlimited

cd /home/hballington/abt

/soft/intel/impi/2019.8.254/intel64/bin/mpiexec \
/home/hballington/abt/src/seq/abt \
-lambda 0.532 \
-rbi 1.3117 \
-ibi 0 \
-rec 8 \
-cmethod read \
-cft obj \
-cfn test_sphere.obj \
-afn test_sphere_apertures.dat \
-rot euler 0.01 90 0.01 \
-jobname testing_a_sphere \
-cc_hex_l 10 \
-cc_hex_hr 5 \
-cc_hex_nfhr 4 \
-cc_hex_pfl 10 \
-cc_hex_nfpl 8 \
-cc_hex_pher 1 \
-cc_hex_pper 2 \
-cc_hex_nscales 1 \
-cc_hex_cls 1 \
-cc_hex_sds 0 \
-theta 0 1 180 \
-phi 0 2 360 \
-mt \
> log2
