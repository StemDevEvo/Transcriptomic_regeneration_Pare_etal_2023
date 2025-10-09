#!/usr/bin/bash

. ../config.txt

# Determine if reads are paired ended based on the sample file.
pairedend=true
nbr_of_col=$( head -n 1 "$sample_file" | tr -s " " | tr " " "\t" | tr "\t" "\n" | wc -l )

if [ $nbr_of_col = 3 ]; then
    pairedend=false
fi

if [ $nbr_of_col -gt 4 ]; then
    echo "There is an error with the sample file, please verify that there is no spaces in the file names and there are only 3 or 4 columns."
    exit
fi

# get number of line, wc -l not working as intended.
nbr_of_samples=$( sed -n '$=' $sample_file )

# Create output directories
[ ! -d ${output}/fastp_trimmed/ ] && mkdir ${output}/fastp_trimmed/

[ ! -d ${output}/fastp_reports/ ] && mkdir ${output}/fastp_reports/

# Create new sample file with the trimmed reads location.
> ${output}/fastp_trimmed/sample_file.txt 

#For each, either pair of sample or unique sample, perform fastp.
for i in $(eval echo "{1..${nbr_of_samples}}");  
do
    if [ "$pairedend" = true ]; then
        read_1=$( head -n $i $sample_file | tail -n 1 | tr -s ' ' | tr " " "\t" | cut  -f 3 )
        read_2=$( head -n $i $sample_file | tail -n 1 | tr -s ' ' | tr " " "\t" | cut  -f 4 )

        bash ../src/02_fastp.sh \
        ${output}/fastp_trimmed/ \
        ${output}/fastp_reports/$( echo ${read_1::-7} | xargs basename) \
        $read_1 \
        $read_2 \
        ${read_1::-7} \
        "$fastp_parameters" \
        $pairedend &

        # update sample file by adding the pair of read currently being processed.
        # It does not need to wait the fastp script to finish.
        bash ../src/update_sample_file.sh "$sample_file" "$output" "$i" &

    else
        read=$( head -n $i $sample_file | tail -n 1 | tr -s ' ' | tr " " "\t" | cut  -f 3 )

        bash ../src/02_fastp.sh \
        ${output}/fastp_trimmed/ \
        ${output}/fastp_reports/$( echo ${read::-7} | xargs basename ) \
        $read \
        $read \
        ${read::-7} \
        $fastp_parameters \
        $pairedend &
        # update sample file by adding the read currently being processed.
        # It does not need to wait the fastp script to finish.
        bash ../src/update_sample_file.sh "$sample_file" "$output" "$i" &
    fi
done