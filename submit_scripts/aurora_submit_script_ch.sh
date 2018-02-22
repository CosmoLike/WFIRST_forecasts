#!/bin/bash
#PBS -V
#PBS -q shortq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=3:00:00
#PBS -J 1-6105
#PBS -N W1st_1-5k
#PBS -e /aurora_nobackup/cosmos/heinrich/WFIRST_forecasts/output2/
#PBS -o /aurora_nobackup/cosmos/heinrich/WFIRST_forecasts/output2/

cd $PBS_O_WORKDIR
/home/teifler/Dropbox/cosmolike/top-level/WFIRST/./compute_covariances_fourier $PBS_ARRAY_INDEX >& /aurora_nobackup/cosmos/heinrich/WFIRST_forecasts/job_output.log

