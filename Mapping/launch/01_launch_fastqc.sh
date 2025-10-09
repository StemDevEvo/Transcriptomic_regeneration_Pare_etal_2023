#!/usr/bin/bash

. ../config.txt


fastqc_output="$output"/fastqc

# Determine if samples are paired ended based on the sample file.
pairedend=true
nbr_of_col=$( head -n 1 "$sample_file" | tr -s " " | tr " " "\t" | tr "\t" "\n" | wc -l )

if [ $nbr_of_col = 3 ]; then
    pairedend=false
fi

if [ $nbr_of_col -gt 4 ]; then
    echo "There is an error with the sample file, please verify that there is no spaces in the file names and there are only 3 or 4 columns."
    exit
fi

# Create output folder if not exist
[ ! -d $output ] && mkdir $output

[ ! -d $fastqc_output ] && mkdir $fastqc_output

# Get the directory where the reads are
read_path=$( head -n 1 $sample_file | tr -s ' ' | tr " " "\t" | cut -f3 | xargs dirname)

# Launch fastqc, using sbatch system
#sbatch -A $project ../src/01_fastqc.sh $fastqc_output $read_path &
bash ../src/01_fastqc.sh $fastqc_output $read_path &