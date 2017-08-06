#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -q shortq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=3:00:00
#PBS -J 1-6105
#PBS -N W1st_1-5k
#PBS -e /aurora_nobackup/sunglass/teifler/output2/
#PBS -o /aurora_nobackup/sunglass/teifler/output2/

cd $PBS_O_WORKDIR
/home/teifler/Dropbox/cosmolike/top-level/WFIRST/./compute_covariances_fourier9 $PBS_ARRAY_INDEX >& /aurora_nobackup/sunglass/teifler/job_output.log


