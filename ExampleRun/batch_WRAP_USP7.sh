#!/bin/bash
#PBS -l walltime=48:00:00                       # 48 hour runtime
#PBS -l nodes=1:ppn=1,pmem=64gb,mem=64gb        # Request 4 node, 2 processor per node, and 64gb of memory
#PBS -t 1-10                                    # Run 10 processes at once


export PATH=$PATH:/share/dept_bbsr/Projects/Shaw_Timothy/3352_Splicing_Pipeline_2021/scripts/src/STAR-2.7.6a/source:/share/dept_bbsr/Projects/Shaw_Timothy/3352_Splicing_Pipeline_2021/scripts/DRPPM/DRPPM-master/export/:/share/dept_bbsr/Projects/Shaw_Timothy/3352_Splicing_Pipeline_2021/scripts/src/RSeQC-2.6.4/install/share/apps/python-2.7.9/bin/:/share/Lab_Shaw/software/bin/

# Establish starting directory
cd /share/Lab_Shaw/projects/ShawLab/USP7_shRNA_A_Project/USP7_rerun2/

# Assigns an ArrayID (1-10) to each line (p) in the execute_everything.sh script, where each line is a seperate shell script
line=`sed -n "${PBS_ARRAYID}p" execute_everything.sh`

# Essential modules for scripts
module load gcc/5.5.0
#module load glibc/2.14.1
module load samtools/1.1
module load fastqc/0.11.7
module load python/2.7.9
#module load python/3.7.6
module load R/3.5.1
module load bedtools2/2.27.1

# Each line/shell script of the excute_everything.sh script if run through here as a seperate but grouped job
${line}
