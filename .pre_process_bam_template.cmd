#!/bin/bash

#SBATCH --job-name=samplenumber
#SBATCH --output=logs/filename_logs.out
#SBATCH --error=logs/filename_error.out
#SBATCH --ntasks=1

hostname >&2
printf START >&2; uptime >&2
date >&2

## Dates variable for log files name
date_var=$(date +%Y%m%d%H%M%S)

## Change log file name with date time
mv ./logs/filename_logs.out ./logs/filename_logs_${date_var}.out
mv ./logs/filename_error.out ./logs/filename_error_${date_var}.out

## Remove the job file itself if it is cancelled before launched
exit_abnormal() {
    rm -f "$0"
}

trap 'exit_abnormal' EXIT
trap 'exit_abnormal' ERR
trap 'exit_abnormal' INT
trap 'exit_abnormal' TERM

## Load modules
module load samtools/latest java/12.0.2 gatk/4.1.8.1 

## Run script
time bash ./src/pre_process_bam_parallel_jobs.sh -i template

printf END >&2; uptime >&2
