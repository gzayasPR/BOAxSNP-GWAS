#!/bin/bash
# Create a function for waiting job to finish (avoid to use all memory of the group in the server)
# made by Giovanni Ladeira, July 12, 2024
# Thank you Giovanni
function wait_for_job_completion() {
    # Crete a local variable (job_id) to store the first argument passed to the function
    local job_id=$1
    local time_min=$2
    echo "          Waiting for job $job_id to complete..."
    # While loop only stops when break or exit commands are executed
    while true; do
        # Check if the job is still in the queue
        job_status=$(squeue -j $job_id -h -o %T 2>/dev/null)
        if [[ -z $job_status ]]; then
            echo ""
            break
        elif [[ $job_status == "FAILED" ]]; then
            echo "          Job $job_id failed."
            exit 1
        else
            # Job is still running or pending
            echo "          Job $job_id is $job_status. Checking again in ${time_min} minutes..."
            sleep ${time_min}m
        fi
    done
}

function wait_for_job_completion_less_output() {
    # Crete a local variable (job_id) to store the first argument passed to the function
    local job_id=$1
    local time_min=$2
    # While loop only stops when break or exit commands are executed
    while true; do
        # Check if the job is still in the queue
        job_status=$(squeue -j $job_id -h -o %T 2>/dev/null)
        if [[ -z $job_status ]]; then
            break
        elif [[ $job_status == "FAILED" ]]; then
            echo "          Job $job_id failed."
            exit 1
        else
            # Job is still running or pending
            echo "          Job $job_id is $job_status. Checking again in ${time_min} minutes..."
            sleep ${time_min}
        fi
    done
}

# Function to calculate heritability
calculate_heritability() {
    local var_g=$1
    local var_e=$2
    echo "scale=4; $var_g / ($var_g + $var_e)" | bc
}


function wait_job_check_progress_GWAS() {
    # Input arguments
    local job_file=$1
    local time_min=$2
    local start_time_gwas=$3
    local result_files=$4

    # Main loop to track job statuses
    while true; do
        # Read job IDs from the file
        local job_ids=($(cat "$job_file"))
        echo "Checking jobs: ${job_ids[@]}"
        
        # Track incomplete jobs
        local remaining_jobs=()
        
        for job_id in "${job_ids[@]}"; do
            # Check if the job is still in the queue
            job_status=$(squeue -j "$job_id" -h -o %T 2>/dev/null)
            if [[ -n $job_status ]]; then
                if [[ $job_status == "FAILED" ]]; then
                    echo "          Job $job_id failed."
                    exit 1
                else
                    remaining_jobs+=("$job_id")
                fi
            fi
        done

        # If there are no remaining jobs, break the loop
        if [[ ${#remaining_jobs[@]} -eq 0 ]]; then
            echo "All jobs completed."
            break
        fi

        # Update job file with remaining jobs
        echo "${remaining_jobs[@]}" > "$job_file"
        echo "Jobs still running: ${remaining_jobs[@]}"
        
        # Capture the intermediate time after job check
        check_time_gwas=$(date +%s)
        
        # Calculate time difference from the start
        elapsed_seconds=$((check_time_gwas - start_time_gwas))
        
        # Convert seconds to days, hours, minutes, and seconds
        days=$((elapsed_seconds / 86400))
        hours=$(( (elapsed_seconds % 86400) / 3600 ))
        minutes=$(( (elapsed_seconds % 3600) / 60 ))
        seconds=$((elapsed_seconds % 60))
        
        # Calculate the total number of lines in all files within the directory
        num_lines=$(wc -l "${result_files}"/* 2>/dev/null | tail -1 | awk '{print $1}')
        
        # Count the number of files in the directory
        num_files=$(ls -1 "${result_files}" 2>/dev/null | wc -l)
        
        # Calculate the number of SNPs by subtracting the number of files from the total number of lines
        num_snps=$((num_lines - num_files))
        
        # Display the number of SNPs processed so far and the elapsed time
        echo "$num_snps SNPs done in ${days} days, ${hours} hours, ${minutes} minutes, ${seconds} seconds"
        
        # Wait before checking again
        echo "Checking again in ${time_min} minutes..."
        sleep "${time_min}m"
    done
}