#!/bin/bash

##############################################################

## STRUCTURAL VARIANTS PRE-PROCESS PIPELINE ##

Author="Arnau Soler Costa"
Creation_date="09/02/2023"
Update_date="02/03/2023"
Version="v0.0.1"

##############################################################

# Script to pre-process the BAM files. 

# This pipeline follows the GATK (gatk-4.2.6.0) best practices. The execution of this script allows to obtain a final BAM file: sorted (samtools v.1.3.1), marked the duplicates (GATK, gatk-4.2.6.0) and perform an optional step called 'Base Quality Score Recalibration' (BQSR from GATK). For the 'BQSR' step we are using the information on gnomAD (Broad Institute) structural variants VCF as the known sites of real variants (https://gnomad.broadinstitute.org/downloads/#v2-structural-variants).

# The only requirements for the input BAM file are the need of a header with at least the '@RG' tag in it defined.

# The reference genome need to have in the same folder specified the '.fasta', '.fasta.dict' and '.fasta.fai' files.

# Indexed gnomad SV known sites VCF needed.

# Calculate statistics and depth coverage additional file scripts needed in order to perform the final statistics.

# Paths to gnomad VCF and calculate final statistics files, as well as path to yq program, must be specified in this script. 'Configuration variables and paths' section.

# The pipeline has 'logic coding' skipping the steps of sorting and marking the duplicates if the initial BAM file is already sorted and marked and if it's wanted. Moreover, it's possible to select if you want to sort and/or mark the duplicates again, or force the script to do it if there is a problem with the initial BAM. It's possible to choose if the 'BQSR' step is run or not, and when running this step, it's possible to choose if removing all intermediate files or keeping the BAM file sorted and marked the duplicates while the one sorted, marked the duplicates and 'BQSR', among other parameters.

# The scipt allows single and multi-sample BAM files.

# The script has a set of commands to return a set of files with information and result statistics of the initial and final BAMs. It's also possible to choose to no return these final informative files. 

# The script creates a new and final directory where is possible to find all the results.

# No read is trimmed or removed in this pipeline.

# Use the 'conf.yaml' file to set all the parameters.

# Please, do not use large path directories to specify where are the files or where to store the results.

# Note: This script need the modules: GATK (gatk-4.2.6.0), samtools (v.1.3.1), java (v.12.0.2) and yq (v4.33.2). For the GATK, samtools and java modules the program needs to be in the '$PATH' variable and for the yq program is possible to set the path in the '$yq_path' variable a few lines below (feel free to change/remove the path of the yq program in this script and add it to the '$PATH' variable).


##############################################################
##############################################################

# Global variables

# Set the color variable
red='\033[0;31m'
green='\033[0;32m'
yellow='\033[0;33m'
# Clear the color after that
clear='\033[0m'

# Function to print a help message.
usage() {
  echo -e "\nUsage: $0 [OPTIONS]

            Use a 'conf.yaml' file to set the parameters.  
            
            REQUIRED SETTINGS:

                INPUT_BAM_FILE: Path to the bam to be processed (ie. filename path). No default.
            
                REFERENCE_GENOME: Path to reference genome of the input bam file (ie. b37, hg19, GRCh38,...). No default.
            
                TMP_DIR: Path to the directory for the temporary files. (ie. directory path). No default. 

 
            PROGRAM SETTINGS:

                OUTPUT_DIR: Path to the directory for the output files. (ie. directory path). Same directory as default.

                SORT_THREADS: Set the number of additional threads to use for the 'samtools sort' function (-@/--threads) (ie. 4,8,16). Default 8 threads.
   
                SORT_MEM: Set the maximum memory per thread to use in the 'samtools sort' function (-m) (ie. 750K, 500M, 1G). Default 500M.
 

            OPTIONAL SETTINGS:
   
                BQSR: If TRUE, BQSR will be performed (ie. true/false). Default 'false'.
            
                FORCE_SORT_COORDINATES: Force sort the bam file by coordinates (ie. true/false). Default 'false'.
            
                FORCE_NOT_MARK_DUPLICATES: Force not performing mark the duplicates of the bam file. The user has to be sure that the BAM file has the duplicates already marked by the MarkDuplicates (GATK) tool (version gatk-4.2.6.0). Default 'false'. 
            
                RM_SORT_MARKED_FILE: Remove the sorted and marked bam file when you perform the BQSR as a final step (ie. true/false). Default 'false'.
            
                INITIAL_INFO: Show the initial information of the input bam file. Default true. (ie. true/false). Default 'true'.
            
                FINAL_STATISTICS: Perform and show some final statistics from the last obtained bam file and comparisons with the input bam file. Default true. (ie. true/false). Default 'true'.
            
                RM_TMP_FILES: Delete the temporary files. Default true. (ie. true/false). Default 'true'.
            
                INDEX_INPUT_FILE: If TRUE, indexing of the input bam file will be performed (ie. true/false). Default 'false'.

                CLEAN_SAM: To clean the SAM/BAM file by the function CleanSam (Picard - GATK) (ie. true/false). Default 'false'.
            
                VALIDATE_SAM_FILE: To validate the SAM/BAM file by the function ValidateSamFile (Picard - GATK) (ie. true/false). Default 'false'.

            INFORMATION (not in 'conf.yaml' file):

                -h/--help --- Show usage of the script.
            
                -v/--version --- Show version of the script.\n" 1>&2
}

# Function: Exit with error
exit_abnormal() {
  usage
  # Remove tmp folder for sample
  rm -r -f ${TMP_DIR_SAMPLE}

  # Remove the job wrapper
  rm -f ${TMP_DIR}/.job_${filename}.cmd 
  exit 1
}

##############################################################

###################
# Configuration variables and paths
###################

# NOTE:  All the paths must be relative to where the 'pre_process_bam.cmd' launcher is.

# Declare date variable
date_var=$(date +"%F %H:%M:%S")

# Setting variables values file
conf_file='./conf.yaml' 

# Check if 'conf.yaml' exists in the path
if [ -z $conf.yaml ]; then
    echo -e "${red}ERROR${clear}    ${date_var} 'conf.yaml' file is missing."
    exit_abnormal
fi

# Path to Structural Variants known-sites
known_sites="./data/vcf/gnomad_v2.1_sv.sites.vcf.gz"

# Path to calculate statistics and depth coverage programs
calculate_statistics="./src/calculate_statistics.sh"
calculate_depth_coverage="./src/calculate_depth_coverage.sh"

# Path to yq program folder
yq_path="./bin"

# Check if 'yq' exists in the path
if [ -z ${yq_path}/yq ]; then
    echo -e "${red}ERROR${clear}    ${date_var} 'yq' binary is missing. It must be in ${yq_path}, the PATH variable or you can edit the 'pre_process_bam.sh' script to change the path."
    exit_abnormal
fi

# Location where the script runs
initial_path=$(pwd)

# REQUIRED SETTINGS IN CONF.YAML
REFERENCE_GENOME=$(${yq_path}/yq '.REQUIRED_SETTINGS.REFERENCE_GENOME' $conf_file) # Path to the reference genome of the input BAM file
TMP_DIR=$(${yq_path}/yq '.REQUIRED_SETTINGS.TMP_DIR' $conf_file) # Set the path to the temprorary files directory.

# PROGRAM SETTINGS
OUTPUT_DIR=$(${yq_path}/yq '.PROGRAM_SETTINGS.OUTPUT_DIR' $conf_file); OUTPUT_DIR=${OUTPUT_DIR:-${PWD}}
SORT_THREADS=$(${yq_path}/yq '.PROGRAM_SETTINGS.SORT_THREADS' $conf_file); SORT_THREADS=${SORT_THREADS:-8}
SORT_MEM=$(${yq_path}/yq '.PROGRAM_SETTINGS.SORT_MEM' $conf_file); SORT_MEM=${SORT_MEM:-500M}

# OPRIONAL SETTINGS
BQSR=$(${yq_path}/yq '.OPTIONAL_SETTINGS.BQSR' $conf_file); BQSR=${BQSR:-true} # True = Execute BQSR function, False = no execution BQSR.
FORCE_SORT=$(${yq_path}/yq '.OPTIONAL_SETTINGS.FORCE_SORT_COORDINATES' $conf_file); FORCE_SORT=${FORCE_SORT:-false} # True = Force sort bam file by coordinates.
FORCE_NOT_MARK_DUPLICATES=$(${yq_path}/yq '.OPTIONAL_SETTINGS.FORCE_NOT_MARK_DUPLICATES' $conf_file); FORCE_NOT_MARK_DUPLICATES=${FORCE_NOT_MARK_DUPLICATES:-false} # True = Force not to perform mark duplicates the bam file.
RM_SORT_MARKED_FILE=$(${yq_path}/yq '.OPTIONAL_SETTINGS.RM_SORT_MARKED_FILE' $conf_file); RM_SORT_MARKED_FILE=${RM_SORT_MARKED_FILE:-true} # True = Remove the intermidate file 'sort.marked_duplicates'. Only set True if 'BQSR' is also set True.
INITIAL_INFO=$(${yq_path}/yq '.OPTIONAL_SETTINGS.INITIAL_INFO' $conf_file); INIITIAL_INFO=${INITIAL_INFO:-true} # True = Get initial information of the initial bam file (True).
FINAL_STATISTICS=$(${yq_path}/yq '.OPTIONAL_SETTINGS.FINAL_STATISTICS' $conf_file); FINAL_STATISTICS=${FINAL_STATISTICS:-true} # Get final statistics of the final bam file, and comparison with the initial bam file (True).
INDEX_INPUT=$(${yq_path}/yq '.OPTIONAL_SETTINGS.INDEX_INPUT_FILE' $conf_file); INDEX_INPUT=${INDEX_INPUT:-false} # Set 'true' to index the input bam file 
CLEAN_SAM=$(${yq_path}/yq '.OPTIONAL_SETTINGS.CLEAN_SAM' $conf_file); CLEAN_SAM=${CLEAN_SAM:-false} # To execute the CleanSam fuction (true). Default false.
VALIDATE_SAM=$(${yq_path}/yq '.OPTIONAL_SETTINGS.VALIDATE_SAM_FILE' $conf_file); VALIDATE_SAM=${VALIDATE_SAM:-false} # To execute the ValidateSamFile fuction (true). Default false.

##############################################################

# ERROR CONTROL

# Control errors (Add -e to exit with non-zero errors. Add -x to debug. Add -o pipefail, the returned error code of a function will be used as the return code of the whole pipeline. Add -u to traces the unset variables. Add -h to save the location where the command has been used.)
set -uo pipefail

# Debug
#set -x

# Unset all variables
unset filename
unset fullname
unset TIME
unset enddir
unset LOG_SPLIT_FILE
unset LOG_FILE 

# Sets the locale environment variable to "C" in the Bash shell
export LC_ALL=C

##############################################################

# Control variables
while [[ "$#" -gt 0 ]]; do    
                              
  case $1 in                    
    
    -h|--help)
      exit_abnormal
      shift
      ;;
    -v|--version)
      echo "Author: ${Author}"
      echo "Creation date: ${Creation_date}"
      echo "Update date: ${Update_date}"
      echo "Version: ${Version}"
      exit 0
      shift
      ;;
    -i|--input) BAM_FILE="$2";
      if ! [ -f $BAM_FILE ]; then
          echo -e "${red}ERROR${clear}    ${date_var} Sample '$BAM_FILE' in 'samples.txt' file is required or set incorrectly."
          exit_abnormal
      else
          # Set full dir
          fullname=$(echo $(readlink -f "${BAM_FILE}"))
    
          # Check if 'input_sample' (BAM_FILE) parameter was passed to the script properly
          if ! test -f "$fullname"; then

              # Check if the file is a '.bam'
              if [ "${fullname:-4}" = ".bam" ];
              then
                  # Extract header
                  header=$(samtools view --threads 35 -H ${fullname} || :)
                 
                  if ! [ -z "${header}" ]; # Check if the header is good
                  then
                      echo -e "${red}ERROR${clear}    ${date_var} Input file (${fullname}) does not have a valid header and/or body."
                      exit_abnormal
                  fi
              else
                  echo -e "${red}ERROR${clear}    ${date_var} Input file (${fullname}) does not have the proper extension (.bam) or it can't be read properly."
                  exit_abnormal
              fi
        
              echo -e "${red}ERROR${clear}    ${date_var} Can't start the pipeline for ${fullname}, input file or path are not correct."
              exit_abnormal
          
          else
              # Basename of the file
              filename=$(basename $BAM_FILE .bam)
          fi
         
      fi
      shift      
      ;;    
    :)
      echo -e "${red}ERROR${clear}    ${date_var} The parameter $1 is wrong.";
      exit_abnormal
      ;;
    *) # If unknown (any other) option:
      echo -e "${red}ERROR${clear}    ${date_var} One or more parameters are wrong.";
      exit_abnormal # Exit abnormally.
      ;;

  esac
  shift
  
done

# Check if required modules are in the '$PATH' variable
module_samtools=$(echo $PATH | grep -i samtools || :)
module_gatk=$(echo $PATH | grep -i gatk || :)
module_java=$(echo $PATH | grep -iE "(java|JDK)" || :)

if ! [ -n "${module_samtools}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The module SAMTOOLS is not loaded or in the 'PATH'."
    exit 1
elif ! [ -n "${module_gatk}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The module GATK is not loaded or in the 'PATH'."
    exit 1
elif ! [ -n "${module_java}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The module JAVA/JDK is not loaded or in the 'PATH'."
    exit 1
fi

# Check if needed data and aditional code exist
if ! [ -e "${known_sites}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The SV known sites VCF is not in the path specified in this script."
    exit 1
elif ! [ -e "${calculate_statistics}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The 'calculate_statistics.sh' file is not in the path specified in this script."
    exit 1
elif ! [ -e "${calculate_depth_coverage}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The 'calculate_depth_coverage.sh' file is not in the path specified in this script."
    exit 1
else
    # Set full dir
    known_sites=$(echo $(readlink -f "$known_sites"))
    calculate_statistics=$(echo $(readlink -f "$calculate_statistics"))
    calculate_depth_coverage=$(echo $(readlink -f "$calculate_depth_coverage"))
fi

# Check if 'TMP_DIR' parameter was passed to the script
if ! [ -d $TMP_DIR ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Parameter 'TMP_DIR' in 'conf.yaml' is required or set incorrectly."
    exit_abnormal
elif ! [ -w $TMP_DIR ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The 'TMP' folder does not have write permissions."
    exit_abnormal
else
    # Set full dir
    TMP_DIR=$(echo $(readlink -f "$TMP_DIR"))
    # Create tmp folder for sample
    mkdir -p ${TMP_DIR}/${filename}_tmp
    TMP_DIR_SAMPLE="${TMP_DIR}/${filename}_tmp"
fi

# Check 'OUTPUT_DIR' parameter was passed to the script
if ! [ -d $OUTPUT_DIR ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Parameter 'OUTPUT_DIR' is required or set incorrectly."
    exit_abnormal
elif ! [ -w $OUTPUT_DIR ]; then
    echo -e "${red}ERROR${clear}    ${date_var} The 'OUTPUT' folder does not have write permissions."
    exit_abnormal
else
    # Set full output dir
    OUTPUT_DIR=$(echo $(readlink -f "$OUTPUT_DIR"))
fi

# Check 'REFERENCE_GENOME' parameter was passed to the script
if ! [ -e $REFERENCE_GENOME ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Parameter 'REFERENCE_GENOME' is required or set incorrectly."
    exit_abnormal
else
    # Set full dir
    REFERENCE_GENOME=$(echo $(readlink -f "$REFERENCE_GENOME"))
fi

# Check if 'REFERENCE_GENOME' has .fasta, .dict and .fasta.fai
dict_ref="$(dirname "$REFERENCE_GENOME")/$(basename $REFERENCE_GENOME .fasta)"
if ! [ -e "${dict_ref}".fasta.dict ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Missing reference_genome.fasta.dict in 'REFERENCE_GENOME' path folder."
    echo "You can run:
              gatk CreateSequenceDictionary -R reference.fasta -O reference.fasta.dict"
    exit 1
elif ! [ -e ${REFERENCE_GENOME}.fai ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Missing reference_genome.fasta.fai in 'REFERENCE_GENOME' path folder."
    echo "You can run:
              samtools faidx reference.fasta"
    exit 1 
fi

# Check INT parameters
re='^[0-9]+$'

# SORT_THREADS
if ! [[ ${SORT_THREADS} =~ $re ]] ; then # if $SORT_THREADS is an INT
   echo -e "${red}ERROR${clear}    ${date_var} 'SORT_THREADS' is not a number or INT" 
   exit_abnormal
fi

# SORT_MEM
if ! [[ ${SORT_MEM::-1} =~ $re ]] ; then # if $SORT_MEM is an INT
   echo -e "${red}ERROR${clear}    ${date_var} 'SORT_MEM' is not a number or INT"
   exit_abnormal
elif [[ "${SORT_MEM: -1}" != "K" ]] && [[ "${SORT_MEM: -1}" != "M" ]] && [[ "${SORT_MEM: -1}" != "G" ]]; then
   echo -e "${red}ERROR${clear}    ${date_var} 'SORT_MEM' is not ended with K/M/G."
   exit_abnormal
fi

# Check all bolean variables

# BQSR
if ! [[ "${BQSR,,}" == "true" ]] && ! [[ "${BQSR,,}" == "false" ]]; then   # if $BQSR is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'BQSR' must be a boolean, true/false."
    exit_abnormal
fi

# FORCE_SORT_COORDINATES
if ! [[ "${FORCE_SORT,,}" == "true" ]] && ! [[ "${FORCE_SORT,,}" == "false" ]]; then   # if $FORCE_SORT is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'FORCE_SORT_COORDINATES' must be a boolean, true/false."
    exit_abnormal
fi

# FORCE_NOT_MARK_DUPLICATES
if ! [[ "${FORCE_NOT_MARK_DUPLICATES,,}" == "true" ]] && ! [[ "${FORCE_NOT_MARK_DUPLICATES,,}" == "false" ]]; then   # if $FORCE_NOT_MARK_DUPLICATES is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'FORCE_NOT_MARK_DUPLICATES' must be a boolean, true/false."
    exit_abnormal
fi

# RM_SORT_MARKED_FILE
if ! [[ "${RM_SORT_MARKED_FILE,,}" == "true" ]] && ! [[ "${RM_SORT_MARKED_FILE,,}" == "false" ]]; then   # if $RM_SORT_MARKED_FILE is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'RM_SORT_MARKED_FILE' must be a boolean, true/false."
    exit_abnormal
fi

# INITIAL_INFO
if ! [[ "${INITIAL_INFO,,}" == "true" ]] && ! [[ "${INITIAL_INFO,,}" == "false" ]]; then   # if $INITIAL_INFO is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'INITIAL_INFO' must be a boolean, true/false."
    exit_abnormal
fi

# FINAL_STATISTICS
if ! [[ "${FINAL_STATISTICS,,}" == "true" ]] && ! [[ "${FINAL_STATISTICS,,}" == "false" ]]; then   # if $FINAL_STATISTICS is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'FINAL_STATISTICS' must be a boolean, true/false."
    exit_abnormal
fi

# INDEX_INPUT
if ! [[ "${INDEX_INPUT,,}" == "true" ]] && ! [[ "${INDEX_INPUT,,}" == "false" ]]; then # if $INDEX_INPUT is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'INDEX_INPUT' must be a boolean, true/false."
    exit_abnormal
fi

# CLEAN_SAM
if ! [[ "${CLEAN_SAM,,}" == "true" ]] && ! [[ "${CLEAN_SAM,,}" == "false" ]]; then # if $CLEAN_SAM is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'CLEAN_SAM' must be a boolean, true/false."
    exit_abnormal
fi

# VALIDATE_SAM
if ! [[ "${VALIDATE_SAM,,}" == "true" ]] && ! [[ "${VALIDATE_SAM,,}" == "false" ]]; then # if $VALIDATE_SAM is not boolean:
    echo -e "${red}ERROR${clear}    ${date_var} 'VALIDATE_SAM_FILE' must be a boolean, true/false."
    exit_abnormal
fi

##################################################################

##############################
# Start of the pipeline
##############################

echo -e "${green}##### START OF THE PIPELINE #####${clear}" 

# Set names and path directories

# Need to be 'sudo' to run the script
#[[ $(id -u) -eq 0 ]] || { echo >&2 "Must be root to run the script properly"; exit 1; }

# Set variable control exit
set -e

# Debug option
#set -x

# Config parameters for logs (info. https://serverfault.com/questions/103501/how-can-i-fully-log-all-bash-scripts-actions)
# Save original stdout and stderr file descriptors
exec 3>&1 4>&2

# Set up trap to restore original file descriptors on signal
trap 'exec 2>&4 1>&3' 0 1 2 3

# Set up signal handler for INT signal (Ctrl + C). Exits immediatly when pressing (Ctrl + C). If not, thanks to the while loop, the user would need to press it as many times as samples in samples.txt.
trap "echo; exit" INT

############
     
# Set some variables and create needed files/folders

# Create results folder (end directory)
TIME=$(date +%Y%m%d%H%M%S)
enddir=${OUTPUT_DIR}/${filename}_${TIME} #.${random_number}
mkdir $enddir

# Create error logs split file
touch ${enddir}/${filename}_split_sample.log
LOG_SPLIT_FILE="${enddir}/${filename}_split_sample.log"

# Function to display error message and exit with non-zero status code
handle_error() {
  
# Save status of the error
status=$?
  
# Check if the split file log is empty
if [ -s "${LOG_SPLIT_FILE}" ]; then
    echo -e "${red}ERROR${clear}    ${date_var} Command failed with status ${status}. See '${LOG_SPLIT_FILE}'" 
else
    echo -e "${red}ERROR${clear}    ${date_var} Command failed with status ${status}"
    rm -f ${LOG_SPLIT_FILE}
fi
  
# Check if the created directory is empty
if ! [ "$(ls -A ${enddir})" ]; then
    # Folder is empty. Remove it
    rm -r -f ${enddir}
fi

# If we have more than one sample to be processed
if [ `cat ${enddir}/list_of_samples.txt | wc -l` -gt 1 ]
then
    # Continue the while loop (next sample)
    continue
else
    #rm -f ${OUTPUT_DIR}/.samples_tmp.txt
    rm -r -f ${TMP_DIR_SAMPLE}
    rm -f ${TMP_DIR}/.job_${filename}.cmd # Remove the job wrapper
    exit 1
fi
}

# Set up trap to call handle_error function on error
trap 'handle_error' ERR
    
# Trap the final message to only echo if the script has been run successfully
    
exit_normal() {
  
    # Remove temporary files
    N_TMP_FILES=$(ls ${TMP_DIR} | wc -l) # Set number of temporary files 
    echo -e "${green}INFO${clear}    ${date_var} Removing ${N_TMP_FILES} temporary file(s)" 
    
    rm -r -f ${TMP_DIR_SAMPLE} 
        
    rm -f ${TMP_DIR}/.job_${filename}.cmd # Remove the job wrapper

    echo -e "${green}INFO${clear}    ${date_var} Check the results in '${enddir}' folder."
    echo -e "${green}##### END OF THE PIPELINE #####${clear}" >&3
}

trap 'exit_normal' EXIT

# Echo the processing file
echo -e "${green}INFO${clear}    ${date_var} PROCESSING FILE: ${filename}.bam" >&3

###############################################################

####################
## Check multi-sample bam file
####################
    
# Find if the bam is multi sample
multi_sample=$(samtools view -H ${fullname} | grep "^@RG" | grep "SM" | sed 's/^.*\(SM:.*\).*$/\1/' | cut -f1 | sort | uniq | wc -l) || multi_sample=1 # The '||' is a OR condition to assing 1 if the first statement is not accomplished or is empty assign 1 to the variable

# If it's multi-sample bam, split it
if [[ $multi_sample -gt 1 ]];then
        
    echo -e "${green}INFO${clear}    ${date_var} Split samples logs file for '${filename}.bam'" &>> ${enddir}/${filename}_split_sample.log
    echo -e "${green}INFO${clear}    ${date_var} The input BAM file is a multi-sample file" &>> ${enddir}/${filename}_split_sample.log
    echo -e "${green}INFO${clear}    ${date_var} Splitting the bam file by samples (according to the '@RG' tag and 'SM' field in the BAM header)" &>> ${enddir}/${filename}_split_sample.log

    # Split reads function
    echo -e "${green}INFO${clear}    ${date_var} BAM file has more than one sample" >&3
    echo -e "${green}INFO${clear}    ${date_var} Splitting BAM file" >&3
        
    gatk SplitReads -I ${fullname} -O ${enddir} --split-sample true --create-output-bam-index false --read-validation-stringency LENIENT --tmp-dir ${TMP_DIR_SAMPLE} --add-output-sam-program-record true &>> ${enddir}/${filename}_split_sample.log
        
    echo -e "${green}INFO${clear}    ${date_var} Splitting of BAM file done. See '${enddir}/${filename}_split_sample.log'" >&3
        
    # List of samples
    find ${enddir} -type f -name "*bam" > ${enddir}/list_of_samples.txt
        
    # Remove @RG tag from other samples in each new BAM file
    for bad_head_bam in `cat ${enddir}/list_of_samples.txt`
    do
        # Found the name of each SM in @RG
        sample_name=$(samtools view -H ${bad_head_bam} | grep "^@RG" | grep "SM" | sed 's/^.*\(SM:.*\).*$/\1/' | cut -f1 | cut -d ":" -f2 | sort | uniq)

        for sample_tag in ${sample_name}
        do
            filename_string=$(basename ${bad_head_bam} .bam)

            # If the sample_tag is in the name of the file, do not remove it from the header
            if ! [[ "$filename_string" == *"$sample_tag"* ]];
            then
                # Extract the original header, remove the lines where the sample name that we dont want is found and change the header of the original file
                samtools view -H ${bad_head_bam} > ${TMP_DIR_SAMPLE}/.header.sam

                sed  "/^@RG.*SM:${sample_tag}/ d" ${TMP_DIR_SAMPLE}/.header.sam > ${TMP_DIR_SAMPLE}/.header_corrected.sam

                samtools reheader ${TMP_DIR_SAMPLE}/.header_corrected.sam ${bad_head_bam} > ${TMP_DIR_SAMPLE}/.tmp_reheader.bam

                # Change the name to the original filename
                mv ${TMP_DIR_SAMPLE}/.tmp_reheader.bam ${bad_head_bam}

                # Remove intermediate files
                rm -f ${TMP_DIR_SAMPLE}/.header.sam ${TMP_DIR_SAMPLE}/.header_corrected.sam
            fi
        done
    done

else
    # The input BAM file is a single sample
    echo ${fullname} > ${enddir}/list_of_samples.txt	
fi

# Remove split_sample.log file if it's empty (single sample input BAM file)
if ! [ -s "${enddir}/${filename}_split_sample.log" ]; then
    rm -f ${enddir}/${filename}_split_sample.log
fi

###############################################################
###############################################################

## Start of pipeline

for sample in `cat ${enddir}/list_of_samples.txt`
do

    # Basename of the file
    filename=$(basename $sample .bam)
      
    # Full path to file
    fullname=$(echo $(readlink -f "$sample"))
   
    # Index file name
    INDEX="${fullname}.bai"

    # Create error logs file
    LOG_FILE="${enddir}/${filename}.log"

    # Redirect stdout to file.log and stderr to stdout
    exec 1>${LOG_FILE} 2>&1

    # Function to display error message and exit with non-zero status code
    handle_error() {
        status=$?
        echo -e "${red}ERROR${clear}    ${date_var} Command failed with status ${status}. See '${LOG_FILE}'" >&3
        rm -f ${TMP_DIR}/.job_${filename}.cmd # Remove the job wrapper
        exit 1
    }

    # Set up trap to call handle_error function on error
    trap 'handle_error' ERR
        
    # Echo the processing file
    if [[ $multi_sample -gt 1 ]];then
        echo -e "${green}INFO${clear}    ${date_var} PROCESSING FILE: ${filename}.bam" >&3
    fi

    # Title for the log file
    echo -e "${green}INFO${clear} LOGS FILE FOR '${filename}.bam'"

    # Everything below will go to the file '$sample.log':
    # Only work inside the loop (because is inside it)?
    # Adding '>&3' to the end of the command, the output goes to the console, not the 'log' file 
    # Other way: Adding '&>> ${enddir}/${filename}.log' to the end of each command

    # Print all the parameters values/settings.
    echo -e "${green}INFO${clear} CONFIGURATION:

                          REQUIRED SETTINGS:

           	            INPUT_BAM_FILE --- ${fullname}
                    
                            REFERENCE_GENOME --- ${REFERENCE_GENOME}
                      
                            TMP_DIR --- ${TMP_DIR} 


               	          PROGRAM SETTINGS:

                            OUTPUT_DIR --- ${OUTPUT_DIR}

                            SORT_THREADS --- ${SORT_THREADS}

                            SORT_MEM --- ${SORT_MEM}

         
                          OPTIONAL SETTINGS:

                            BQSR --- ${BQSR}
                    
                            FORCE_SORT_COORDINATES --- ${FORCE_SORT} 
                    
                            FORCE_NOT_MARK_DUPLICATES --- ${FORCE_NOT_MARK_DUPLICATES} 
                    
                            RM_SORT_MARKED_FILE --- ${RM_SORT_MARKED_FILE}
                    
                            INITIAL_INFO --- ${INITIAL_INFO}
                    
                            FINAL_STATISTICS --- ${FINAL_STATISTICS}
                    
                            INDEX_INPUT_FILE --- ${INDEX_INPUT}

                            CLEAN_SAM --- ${CLEAN_SAM}
                    
                            VALIDATE_SAM_FILE --- ${VALIDATE_SAM}\n"

    #########################

    # Read Groups tag check + creating it if necessary
    rg_tags=$(samtools view -H ${fullname} | grep "@RG" || :)

    if ! [ -n "${rg_tags}" ]; then
        echo -e "${red}ERROR${clear}    ${date_var} There is not Read Group tags (RG tags) in the BAM file, which are required for the next steps. Please add the @RG tags to the BAM file and rerun the script."
        echo -e "${green}INFO${clear}    ${date_var} You can use the next command as an example:"
        echo "          gatk AddOrReplaceReadGroups -I input.bam -O output.bam -RGLB lib1 -RGPL platform1 -RGPU unit1 -RGSM sample"
        handle_error
    fi

    # Clean sam 
    if [ "${CLEAN_SAM,,}" = "true" ];
    then
        echo -e "${green}INFO${clear}    ${date_var} Start CleanSam"
             
        gatk CleanSam --INPUT ${fullname} --OUTPUT ${enddir}/${filename}_clean.bam --VALIDATION_STRINGENCY LENIENT --TMP_DIR ${TMP_DIR_SAMPLE}
        mv ${enddir}/${filename}_clean.bam ${enddir}/${filename}.bam # Change the file name to the initial one
            
        echo -e "${green}INFO${clear}    ${date_var} CleanSam done"

    fi

    # Check if the input BAM file is in the 'enddir' create due to being edited or changed by the CLEAN_SAM function
    if [ -f ${enddir}/${filename}.bam ] && [ "${CLEAN_SAM,,}" = "true" ];
    then
        fullname=${enddir}/${filename}.bam
    fi

    #########################

    # Validate Sam file before or after cleaning
    if [ "${VALIDATE_SAM,,}" = "true" ];
    then
            #echo "Start ValidateSamFile - MODE=SUMMARY"
            #/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk ValidateSamFile -I ${fullname} --MODE SUMMARY --TMP_DIR ${TMP_DIR
            
        echo -e "${green}INFO${clear}    ${date_var} Example provided to perform the Validation of the BAM file by the user using ValidateSamFile (Picard - GATK) tool and be able to fix the errors/warnings if needed:
                  
                  SUMMARY MODE:

                          ValidateSamFile -I input.bam --MODE SUMMARY

                  VERBOSE MODE (seeing only errors):

        	          ValidateSamFile -I input.bam --MODE VERBOSE --IGNORE_WARNINGS true
               
        	  VERBOSE MODE (seeing only warnings):

                          ValidateSamFile -I input.bam --MODE VERBOSE --IGNORE type
        	
        	  https://gatk.broadinstitute.org/hc/en-us/articles/360036854731-ValidateSamFile-Picard-
                  https://sites.google.com/a/broadinstitute.org/legacy-gatk-documentation/solutions-to-problems/7571-Errors-in-SAMBAM-files-can-be-diagnosed-with-ValidateSamFile" 

    fi

    # Index the input BAM file if necessary
    # Index file if its not already indexed
    if test -f "${INDEX}"; then
        echo -e "${yellow}WARNING${clear} ${date_var} '${filename}.bam.bai' already exists. Skipping the index of the input file"
    else
        if [ "${INDEX_INPUT,,}" = "true" ];
        then

            echo -e "${green}INFO${clear}    ${date_var} Start indexing '${filename}.bam'" 

            samtools index ${fullname}

            echo -e "${green}INFO${clear}    ${date_var} End index"
 
        else
            echo -e "${yellow}WARNING${clear} ${date_var} No indexed file of input bam in the folder"
        fi

    fi

    ##############################################################
    ##############################################################

    ### INFORMATION OF THE BAM ###

    if [ "${INITIAL_INFO,,}" = "true" ];
    then
            
        # Get all PG tags
        pg_tags=$(samtools view -H ${fullname} | grep '^@PG' || :) # The '|| :' is a OR condition to add an empty line if the first statement is not accomplished or is empty

        if [ -n "${pg_tags}" ]; then
            echo "${pg_tags}" > ${enddir}/${filename}.pg_tags
        else
            echo -e "${red}ERROR${clear} There is not Program tags (PG tags) in the BAM file" > ${enddir}/${filename}.pg_tags
        fi 
             
        # Get Mapper
        # Mapper bwa tool
        echo -e "${green}INFO${clear} BWA mapper" > ${enddir}/${filename}.mapper
            
        # Mapper
        mapper_bwa=$(samtools view -H ${fullname} | grep '^@PG' | grep -wo bwa | head -n1 || :) 
           
        if [ -n "${mapper_bwa}" ]; then
            echo "${mapper_bwa}" >> ${enddir}/${filename}.mapper
        else
            echo -e "${red}ERROR${clear} The mapper tool is not BWA or is not in the BAM file" >> ${enddir}/${filename}.mapper
        fi

        # Mapper method
        mapper_method=$(samtools view -H ${fullname} | grep '^@PG' | grep bwa | grep -woE 'aln|mem' | head -n1 || :) 
            
        if [ -n "${mapper_method}" ]; then
            echo "${mapper_method}" >> ${enddir}/${filename}.mapper
        else
            echo -e "${red}ERROR${clear} There is not the BWA mapper method (aln/mem) in the BAM file" >> ${enddir}/${filename}.mapper
        fi

        # Mapper version
        mapper_version=$(samtools view -H ${fullname} | grep '^@PG' | grep bwa | grep -E 'aln|mem' | grep -o "VN:.*" | awk '{print $1}' | cut -c4- | head -n1 || :) 
            
        if [ -n "${mapper_version}" ]; then
            echo "${mapper_version}" >> ${enddir}/${filename}.mapper
        else
            echo -e "${red}ERROR${clear} There is not the BWA mapper version in the BAM file" >> ${enddir}/${filename}.mapper
        fi

        # Read groups
        if [ -n "${rg_tags}" ]; then
            echo "${rg_tags}" > ${enddir}/${filename}.rg_tags
        else
            echo -e "${red}ERROR${clear} There is not Read Group tags (RG tags) in the BAM file: A new should be created for the pipeline" > ${enddir}/${filename}.rg_tags
        fi

        # Get quality threshold
        quality_threshold=$(samtools view -H ${fullname} | grep '^@PG' | grep bwa | grep -o "\-q.*" | cut -c4- | awk '{print $1}' | head -n1 || :)

        if [ -n "${quality_threshold}" ]; then
            echo "${quality_threshold}" > ${enddir}/${filename}.quality_threshold
        else
            echo -e "${red}ERROR${clear} There is not the quality threshold used in the BAM file" > ${enddir}/${filename}.quality_threshold
        fi

        echo -e "${green}INFO${clear}    ${date_var} INITIAL INFORMATION DONE" 

    fi

    ###############################################################
    ###############################################################

    # START OF THE PIPELINE

    echo -e "${green}INFO${clear}    ${date_var} START OF THE PIPELINE" 

    #####################
    ### Sort by coordinate (samtools)
    #####################

    # Getting if the bam file is sorted by coordinates

    SORT_COORDINATE=$(timeout 5 samtools view -H ${fullname} | grep "^@HD" | grep "coordinate" || :)

    if ! [ -n "$SORT_COORDINATE" ] || [ "${FORCE_SORT,,}" = "true" ];
    then

        echo -e "${green}INFO${clear}    ${date_var} START SORT BY COORDINATES" 
            
        # SAMTOOLS SORT
        samtools sort ${fullname} -m ${SORT_MEM} --threads ${SORT_THREADS} -T ${TMP_DIR_SAMPLE} -o ${enddir}/${filename}.sort.bam

        # GATK SORT
        #/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk SortSam INPUT=${fullname} OUTPUT=${enddir}/${filename}sort.bam SORT_ORDER=coordinate --TMP_DIR /vault/bio-scratch/arnau/structural_variants/tmp_files --VALIDATION_STRINGENCY LENIENT 

        echo -e "${green}INFO${clear}    ${date_var} SORT BY COORDINATES DONE" 

    else
        echo -e "${yellow}WARNING${clear} ${date_var} Assuming file is already sorted by coordinates: Sorting is not performed" 
            
        echo -e "${green}INFO${clear}    ${date_var} Copying file to the correct folder"
 
        cp ${fullname} ${enddir}/${filename}.sort.bam

        echo -e "${green}INFO${clear}    ${date_var} File is copyied correctly" 

    fi

    ###############################################################

    #####################
    ### Mark Duplicates (bam file ordered by coordinate) and remove duplicates. Use MarkDuplicatesSpark if you want to run it in parallel.
    #####################

    # Always marked the duplicates if the user is not sure that the BAm file has the duplicates marked by MarkDuplicates GATK tool (version gatk-4.2.6.0) 

    if [[ "${FORCE_NOT_MARK_DUPLICATES,,}" == "false" ]];
    then

        echo -e "${green}INFO${clear}    ${date_var} START MARK DUPLICATES" 
            
        gatk MarkDuplicates -I ${enddir}/${filename}.sort.bam -O ${enddir}/${filename}.sort.marked_duplicates.bam -M ${enddir}/${filename}.sort.marked_duplicates_metrics.txt --TMP_DIR ${TMP_DIR_SAMPLE} --VALIDATION_STRINGENCY LENIENT --ADD_PG_TAG_TO_READS true --ASSUME_SORTED true 

        # Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. For memory issues (ie. --java-options "-Xmx5g -Djava.io.tmpdir=${TMP_DIR}")

        echo -e "${green}INFO${clear}    ${date_var} MARK DUPLICATES DONE" 

    else

        echo -e "${yellow}WARNING${clear} ${date_var} Assuming file is already marked by duplicates: MarkDuplicates (GATK) is not applied" 
            
        echo -e "${green}INFO${clear}    ${date_var} Copying file to the correct folder" 
        cp ${enddir}/${filename}.sort.bam ${enddir}/${filename}.sort.marked_duplicates.bam
        echo -e "${green}INFO${clear}    ${date_var} File is copyied correctly" 

    fi
            
    #######################
    ## Create Index for sort.marked_duplicates file
    #######################

    if [ -e ${enddir}/${filename}.sort.marked_duplicates.bam ];
    then
        echo -e "${green}INFO${clear}    ${date_var} Start create index for 'sort.marked_duplicates' file"

        # Index
        samtools index ${enddir}/${filename}.sort.marked_duplicates.bam
          
        echo -e "${green}INFO${clear}    ${date_var} Indexing done" 
    fi

    ################################################################

    ######################
    ### BaseRecalibrator before
    ######################

    if [ "${BQSR,,}" = "true" ];
    then

        echo -e "${green}INFO${clear}    ${date_var} START BASE RECALIBRATOR"
         
        gatk BaseRecalibrator -I ${enddir}/${filename}.sort.marked_duplicates.bam -R ${REFERENCE_GENOME} --known-sites ${known_sites} -O ${enddir}/${filename}.recal_data.table --tmp-dir ${TMP_DIR_SAMPLE} --read-validation-stringency LENIENT --add-output-sam-program-record true 

        echo -e "${green}INFO${clear}    ${date_var} BASE RECALIBRATOR DONE" 

    #####################
    ## ApplyBQSR
    #####################

        echo -e "${green}INFO${clear}    ${date_var} START APPLY BQSR"
 
        gatk ApplyBQSR -R ${REFERENCE_GENOME} -I ${enddir}/${filename}.sort.marked_duplicates.bam --bqsr-recal-file ${enddir}/${filename}.recal_data.table -O ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam  --add-output-sam-program-record true --create-output-bam-index false --read-validation-stringency LENIENT --tmp-dir ${TMP_DIR_SAMPLE} 

    #--create-output-bam-index true #--preserve-qscores-less-than 10 

        echo -e "${green}INFO${clear}    ${date_var} APPLY BQSR DONE" 

    else
        echo -e "${yellow}WARNING${clear} ${date_var} BQSR is not applied" 
    fi

    ###############################################################

    #######################
    ## Create Index for BQSR
    #######################

    if [ -e ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam ];
    then
        echo -e "${green}INFO${clear}    ${date_var} Start create index for 'BQSR' file" 

        samtools index ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam 

        echo -e "${green}INFO${clear}    ${date_var} Indexing done" 
    fi

    ###############################################################

    #####################
    ## Remove intermediate files
    #####################

    echo -e "${green}INFO${clear}    ${date_var} Removing intermediate files"

    rm -f ${enddir}/${filename}.sort.bam

    if [ "${RM_SORT_MARKED_FILE,,}" = "true" ] && [ "${BQSR,,}" = "true" ];
    then

        rm -f ${enddir}/${filename}.sort.marked_duplicates.bam

    fi

    ###############################################################

    ######################
    # FINAL STATISTICS
    ######################

    if [ "${FINAL_STATISTICS,,}" = "true" ];
    then
            
        echo -e "${green}INFO${clear}    ${date_var} Calculating final statistics"
         
        cd ${enddir} # Get inside the final directory to create a new folder called statistics

        mkdir -p statistics # Create a new folder called statistics if does not exist
            
        if [ "${BQSR,,}" = "true" ] && [ "${RM_SORT_MARKED_FILE,,}" = "false" ];
        then

            # Calculate statistics and depth coverage from the two files

            bash ${calculate_statistics} ${filename}.sort.marked_duplicates 

            bash ${calculate_depth_coverage} ${filename}.sort.marked_duplicates

            bash ${calculate_statistics} ${filename}.sort.marked_duplicates.BQSR 

            bash ${calculate_depth_coverage} ${filename}.sort.marked_duplicates.BQSR 

        elif [ "${BQSR,,}" = "true" ] && [ "${RM_SORT_MARKED_FILE,,}" = "true" ];
        then

            # Calculate statistics and depth coverage from the BQSR file

            bash ${calculate_statistics} ${filename}.sort.marked_duplicates.BQSR 

            bash ${calculate_depth_coverage} ${filename}.sort.marked_duplicates.BQSR
         
        else

            # Calculate statistics and depth coverage from the sort and arked duplicates file
               
            bash ${calculate_statistics} ${filename}.sort.marked_duplicates

            bash ${calculate_depth_coverage} ${filename}.sort.marked_duplicates

        fi

        # Retrun to the original directory
        cd ${initial_path}

         echo -e "${green}INFO${clear}    ${date_var} Final statistics done"
 
    fi

    ###############################################################
    ###############################################################

    # Change file permissions
    #chmod ugo-x ${enddir}/${filename}.sort.marked_duplicates.bam

    echo -e "${green}INFO${clear}    ${date_var} File '${filename}.bam' processing is finished"

    # Restart the set -e command
    set +e 
        
    # End message to final workflow of the sample
    echo -e "${green}INFO${clear}    ${date_var} Processing of '${filename}.bam' done!" >&3

done # End for loop for multiple sample
