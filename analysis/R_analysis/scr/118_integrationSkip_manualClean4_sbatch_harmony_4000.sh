#!/bin/bash
#SBATCH --job-name=integ
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="118_integrationSkip_manualClean4_sbatch_harmony_4000.err"
#SBATCH --output="118_integrationSkip_manualClean4_sbatch_harmony_4000.out"

echo "my job strart now" > 118_integrationSkip_manualClean4_sbatch_harmony_4000.log;

date >> 118_integrationSkip_manualClean4_sbatch_harmony_4000.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_recover2;

Rscript /home/pedrini.edoardo/scr/project_edoardo/220501_scRNAseq_MSbrain_Absinta/R_analysis/scr/118_integrationSkip_manualClean4_sbatch_harmony_4000.R

date >> 118_integrationSkip_manualClean4_sbatch_harmony_4000.log;
echo "all done!!" >> 118_integrationSkip_manualClean4_sbatch_harmony_4000.log