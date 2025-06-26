#! /bin/bash

#SBATCH -o run008.sh.log-%j-%a
#SBATCH -a 1-240
#SBATCH -c 2

# Run the job array
matlab -nodisplay -r "DSC_Simple_SC_run(${SLURM_ARRAY_TASK_ID}, '008', {'B', 'C1', 'C2', 'D1', 'D2'}); exit" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT