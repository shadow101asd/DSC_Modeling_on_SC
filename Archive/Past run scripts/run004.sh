#! /bin/bash

#SBATCH -o run004.sh.log-%j-%a
#SBATCH -a 1-240

# Run the job array
matlab -nodisplay -r "DSC4_SC_runN(${SLURM_ARRAY_TASK_ID}, 4, 240, '004'); exit" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT