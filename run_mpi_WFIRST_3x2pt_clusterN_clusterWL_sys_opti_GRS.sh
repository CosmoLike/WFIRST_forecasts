#!/bin/bash
#PBS -S /bin/bash
#PBS -W group_list=cosmo
#PBS -q qualified
### Set the number of nodes,cores and memory that will be used for this job
### select=3 is the node count, ncpus=28 are the cores in each node,
### mem=168gb is memory per node, pcmem=6gb is the memory per core - optional
#PBS -l select=40:ncpus=28:mem=168GB
#PBS -l place=free:shared
#PBS -l cput=28000:00:00
#PBS -l walltime=25:00:00
#PBS -N WF_all_GRS_sy_o
#PBS -e /home/u17/timeifler/output/
#PBS -o /home/u17/timeifler/output/


cd $PBS_O_WORKDIR

module load python/2
module load mpich/ge/gcc/64/3.2.1
module load openmpi

### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 1120 python runWFIRST_3x2pt_clusterN_clusterWL_allsys_opti_GRS.py
date