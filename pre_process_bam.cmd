#!/bin/sh

#SBATCH --job-name=pre_process
#SBATCH --output=./logs/pre_process_logs.out
#SBATCH --error=./logs/pre_process_error.out
#SBATCH --ntasks=1

hostname >&2
printf START >&2; uptime >&2
date >&2

## Load modules
module load samtools/latest java/12.0.2 gatk/4.1.8.1 

## Debug
##set -x 

## Set the color variable
green='\033[0;32m'
red='\033[0;31m'
## Clear the color after that
clear='\033[0m'

## Create slurm files for each sample/job
echo -e "${green}INFO${clear}    Creating parallel jobs"

## Setting variables values file
conf_file='./conf.yaml'

## Check if 'conf.yaml' exists in the path
if [ -z $conf.yaml ]; then
    echo -e "${red}ERROR${clear}    ${date_var} 'conf.yaml' file is missing."
    exit_abnormal
fi

## Path to yq program folder
yq_path="/slgpfs/projects/slc00/sv_data/bin"

## Check if 'yq' exists in the path
if [ -z ${yq_path}/yq ]; then
    echo -e "${red}ERROR${clear}    ${date_var} 'yq' binary is missing. It must be in ${yq_path}, the PATH variable or you can edit the 'pre_process_bam.sh' script to change the path."
    exit_abnormal
fi

## REQUIRED SETTINGS IN CONF.YAML
SAMPLES=$(${yq_path}/yq '.REQUIRED_SETTINGS.SAMPLES' $conf_file) ## Path to the reference genome of the input BAM file

## Check if the 'SAMPLES' file exists and if it has any problem
## Set full dir
SAMPLES=$(echo $(readlink -f "${SAMPLES}"))
if ! [ -e $SAMPLES ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Parameter 'SAMPLES' in 'conf.yaml' file is required or set incorrectly."
    exit 1
else
   
    ## Remove empty lines in samples.txt file
    sed '/^$/d' $SAMPLES | grep -v "^#" > .edited_samples.txt

    if ! [ -s .edited_samples.txt ]; then ## Check if the variable is empty
        echo -e "${red}ERROR${clear}    ${date_var} Parameter 'SAMPLES' in 'conf.yaml' has no samples in it or are malformed configure."
    else

        ## While loop to iterate among the list of samples
        while read -r line; do

            # Check if the 'line' in the 'samples.csv' file is a path
            if ! [ -f "$line" ]; then
                echo -e "${red}ERROR${clear}    ${date_var} Input sample '${line}' does not have a proper path or it's invalid."
                rm -f .edited_samples.txt
                exit 1
            fi        
 
            ## Input sample
            input_sample=$(echo $(readlink -f "$line"))
            
            ## Check if 'input_sample' (sample) parameter was passed to the script
            if ! test -f "$input_sample"; then

                ## Check if the file is a '.bam'
                if [ "${input_sample:-4}" = ".bam" ];
                then

                    header=$(samtools view --threads 35 -H ${input_sample} || :)
                    
                    if ! [ -z "${header}" ];
                    then
                        echo -e "${red}ERROR${clear}    ${date_var} Input file (${input_sample}) does not have a valid header and/or body."
                        rm -f .edited_samples.txt
                        exit 1
                    fi
                else
                    echo -e "${red}ERROR${clear}    ${date_var} Input file (${input_sample}) does not have the proper extension (.bam) or it can't be read properly."
                    rm -f .edited_samples.txt
                    exit 1
                fi
                echo -e "${red}ERROR${clear}    ${date_var} Can't start the pipeline for ${input_sample}, input file or path are not correct."
                rm -f .edited_samples.txt
                exit 1
            fi

        done < .edited_samples.txt ## End while loop
    fi
fi

counter=1
while IFS= read -r bam_file;
do
 
  # Extract the file names from the absolute paths and count their occurrences
  duplicate_files=$(awk -F/ '{print $NF}' .edited_samples.txt | sort | uniq -cd)

  # Check if there are any duplicate file names
  if [[ -n $duplicate_files ]]; then
      
      # Adding a tab in front of each line (to look good in the output)
      duplicate_files=$(echo "$duplicate_files" | sed 's/^/\t/')
      echo -e "${red}ERROR${clear}   Duplicate file names found:"
      echo "$duplicate_files"
      rm -f .edited_samples.txt
      exit 1
  fi

  filename=$(basename ${bam_file} .bam)
  cp .pre_process_bam_template.cmd .job_${filename}.cmd
  sed -i "s+template+$bam_file+g" .job_${filename}.cmd ## Sed can use any separator rather than '/'
  sed -i "s/filename/${filename}/g" .job_${filename}.cmd
  sed -i "s/number/${counter}/g" .job_${filename}.cmd
  let counter++

done < .edited_samples.txt

rm -f .edited_samples.txt

## Launch all jobs
for job in .job*; do sbatch $job; done

echo -e "${green}INFO${clear}    Done!"

printf END >&2; uptime >&2
