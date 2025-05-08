#! /bin/bash

#SBATCH -o run006.sh.log-%j-%a
#SBATCH -a 1-120
#SBATCH -c 2

# Run the job array
matlab -nodisplay -r "DSC_SC_runN(${SLURM_ARRAY_TASK_ID}, 2, 120, '006', 6); exit" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT