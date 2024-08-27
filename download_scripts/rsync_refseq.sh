#!/bin/bash
# Synchronize with the NCBI RefSeq databases
# the base URL is the same, however, for curl we need https://, while rsync has rsync://
# files to download: '.faa.gz' '.gff.gz' '_assembly_report.txt'

# Usage:
# Start the script in terminal on worker-schwede.scicore.unibas.ch
# cd /scicore/home/schwede/GROUP/gcsnap_db
# nohup ./download_scripts/rsync_refseq.sh >/dev/null 2>&1 &

export WORK_DIR="/scicore/home/schwede/GROUP/gcsnap_db/refseq"
# Define the folder for rsync files
export DESTINATION="$WORK_DIR/data"
# prefix of Database
export PREFIX="GCF"
# summary name
SUMMARY_NAME="assembly_summary_refseq.txt"
# URL for summary
SUMMARY_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/$SUMMARY_NAME"

# Output file to store progress, only the last url to keept that file small
export PROGRESS_FILE="$WORK_DIR/logging/urls_success.txt"
> ${PROGRESS_FILE}
# final output file
export FINAL_FILE="$WORK_DIR/logging/final_summary.txt"
# Failure file to store progress
export FAILURE_FILE="$WORK_DIR/logging/urls_failed.txt"
# clear failure file
> ${FAILURE_FILE}
# files with the URLS
export ALL_URLS="$WORK_DIR/logging/urls_to_sync.txt"

# Number of Cores as paralle GNU
NUM_CORES=10

# Part 1: Get assembly summary.
# -------
echo "Downloading summary ..."
# download
wget -c ${SUMMARY_URL} -P ${WORK_DIR}
# Get all URLs out of this file, this is done once, 
awk -F'\t' '{if(NR>2) print $20}' ${WORK_DIR}/${SUMMARY_NAME} | sort >  ${ALL_URLS}
# remove the na, put it to tmp file and rename file. 
# > directly does not work as awk does not support inplace operation
awk '!/^na$/' "${ALL_URLS}" > temp && mv temp "${ALL_URLS}"

# Record the start time
start_time=$SECONDS

# Part 2: parallel routine
# --------
# load modules
ml parallel
echo "Performing Rsync ..."
# testing with only the first 5099 (head -n 5099 ${ALL_URLS}), for all use cat rsync_urls.txt
cat ${ALL_URLS} | parallel -j ${NUM_CORES} '
    url={};
    if [[ "$url" == "na" ]]; then
        echo "Skipping na URL";
        continue;
    fi
    max_retries=10;
    retries=0;
    success=0;
    # Replace https with rsync in the URL
    rsync_url="${url/https:\/\//rsync://}"
    while [ $retries -lt $max_retries ] && [ $success -eq 0 ]; do
        rsync --timeout=60 --no-relative --partial --include="*.faa.gz" --include="*.gff.gz" --include="*_assembly_report.txt" --exclude="*" "$rsync_url/*" "$DESTINATION" -q  && success=1 || retries=$((retries+1));
        if [ $success -eq 0 ]; then
            echo "Rsync failure for $rsync_url. Retrying $((retries+1))/$max_retries...";
            sleep 10;
        fi
    done;
    if [ $success -eq 0 ]; then
        echo "$rsync_url" | tee -a ${FAILURE_FILE};
    else
		echo "$rsync_url" | tee -a ${PROGRESS_FILE};
    fi
'

# remove the artefacts
mkdir -p ${DESTINATION}/delete_temp
mv ${DESTINATION}/.${PREFIX}* ${DESTINATION}/delete_temp/
rm -rf ${DESTINATION}/delete_temp/

elapsed_time=$(($SECONDS - start_time))

# Convert elapsed time to HH:MM:SS
hours=$(($elapsed_time / 3600))
minutes=$((($elapsed_time % 3600) / 60))
seconds=$(($elapsed_time % 60))

# Format the time to HH:MM:SS
formatted_time=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)

echo "Done with $NUM_CORES cores in $formatted_time" > ${FINAL_FILE}

# count number of assemblies done
success_count=$(grep -c "rsync://ftp" ${PROGRESS_FILE})
echo "Total successful URL processes: $success_count" >> ${FINAL_FILE}

