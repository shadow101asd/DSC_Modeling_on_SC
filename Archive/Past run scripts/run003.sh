#! /bin/bash

#SBATCH -o run003.sh.log-%j-%a
#SBATCH -a 1-500
#SBATCH -c 2

# Run the job array
matlab -nodisplay -r "DSC4_SC_run(${SLURM_ARRAY_TASK_ID},'003'); exit" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT