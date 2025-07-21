#! /bin/bash

#SBATCH -o run014.sh.log-%j-%a
#SBATCH -a 1-96
#SBATCH -c 4

# Pass to MATLAB
matlab -nodisplay -r "DSC_Modular_run_Sequential(${SLURM_ARRAY_TASK_ID}, 1, '014', 'Jupiter'); exit"
