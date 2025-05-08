#! /bin/bash

#SBATCH -o run005.sh.log-%j-%a
#SBATCH -a 1-240
#SBATCH -c 2

# Run the job array
matlab -nodisplay -r "DSC5_SC_run(${SLURM_ARRAY_TASK_ID},'005'); exit" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT