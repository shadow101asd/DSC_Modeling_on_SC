#! /bin/bash

#SBATCH -o run001.sh.log-%j-%a
#SBATCH -a 1-195

# Run the job array
matlab -nodisplay -r "DSC4_SC_run001(${SLURM_ARRAY_TASK_ID}); exit" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT