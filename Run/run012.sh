#! /bin/bash

#SBATCH -o run012.sh.log-%j-%a
#SBATCH -a 1-60
#SBATCH -c 4

# Pass to MATLAB
matlab -nodisplay -r "DSC_Modular_run_Sequential(${SLURM_ARRAY_TASK_ID}, 1, 3, '012'); exit"
