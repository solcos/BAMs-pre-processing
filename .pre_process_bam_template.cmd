#!/bin/bash

#SBATCH --job-name=samplenumber
#SBATCH --output=logs/filename_logs.out
#SBATCH --error=logs/filename_error.out
#SBATCH --ntasks=1


hostname >&2
printf START >&2; uptime >&2
date >&2

## Load modules
module load samtools/latest java/12.0.2 gatk/4.1.8.1 

## Run script
time bash ./src/pre_process_bam_parallel_jobs.sh -i template

printf END >&2; uptime >&2
