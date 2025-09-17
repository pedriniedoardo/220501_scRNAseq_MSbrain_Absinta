#!/bin/bash
#SBATCH --job-name=integ
#SBATCH --account pedrini.edoardo
#SBATCH --mem=370GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=20  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="22_integration_SoupX.err"
#SBATCH --output="22_integration_SoupX.out"

echo "my job strart now" > 22_integration_SoupX.log;

date >> 22_integration_SoupX.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_renv;

Rscript /home/pedrini.edoardo/scr/project_edoardo/220501_scRNAseq_MSbrain_Absinta/R_analysis/scr/22_integration_CXandWM_sbatch.R

date >> 22_integration_SoupX.log;
echo "all done!!" >> 22_integration_SoupX.log