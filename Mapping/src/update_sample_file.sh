# This script take the old sample file and which sample it needs to update and write this information to the new sample file.
# Since it is only read and write to some files, it does not need to be launched with sbatch for exemple.

sample_file=$1
output=$2
read_to_add_i=$3

new_sample_file="$output"/fastp_trimmed/sample_file.txt


line_to_add=$( head -n "$read_to_add_i" "$sample_file" | tail -n 1 | tr -s ' ' | tr " " "\t" | cut -f 1)
line_to_add="$line_to_add"\\t$( head -n "$read_to_add_i" "$sample_file" | tail -n 1 | tr -s ' ' | tr " " "\t" | cut -f 2)

# Determine if paired ended based on original sample file.
pairedend=true
nbr_of_col=$( head -n 1 "$sample_file" | tr -s " " | tr " " "\t" | tr "\t" "\n" | wc -l )

if [ $nbr_of_col = 3 ]; then
    pairedend=false
fi

if [ $nbr_of_col -gt 4 ]; then
    echo "There is an error with the sample file, please verify that there is no spaces in the file names and there are only 3 or 4 columns."
    exit
fi


#If paired ended, 4 lines per samples, else 3 lines per samples
if [ "$pairedend" = true ]; then
    read_1="$output"/fastp_trimmed/$(  head -n "$read_to_add_i" "$sample_file" | tail -n 1 | tr -s ' ' | tr " " "\t" | cut -f 3 | xargs basename )
    read_2="$output"/fastp_trimmed/$(  head -n "$read_to_add_i" "$sample_file" | tail -n 1 | tr -s ' ' | tr " " "\t" | cut -f 4 | xargs basename )
    line_to_add="$line_to_add"\\t"$read_1"\\t"$read_2"

else
    read_1="$output"/fastp_trimmed/$(  head -n "$read_to_add_i" "$sample_file" | tail -n 1 | tr -s ' ' | tr " " "\t" | cut -f 3 | xargs basename )
    line_to_add="$line_to_add"\\t"$read_1"
fi

# Save result to new sample file. (it does not destroy the file but only append the new line)
echo -e $line_to_add >> "$output"/fastp_trimmed/sample_file.txt