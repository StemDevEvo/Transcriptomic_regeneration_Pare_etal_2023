transcriptome=$1
output="$3"
# take the new sample file with the trimmed reads location.
samples_file="${output}"/fastp_trimmed/sample_file.txt 
thread=$4
gene_trans_map=$5
fragment_length=$6

# Determine if paired ended based on the sample file.
pairedend=true
nbr_of_col=$( head -n 1 "$samples_file" | tr -s " " | tr " " "\t" | tr "\t" "\n" | wc -l )

if [ $nbr_of_col = 3 ]; then
    pairedend=false
fi


if [ $nbr_of_col -gt 4 ]; then
    echo "There is an error with the sample file, please verify that there is no spaces in the file names and there are only 3 or 4 columns."
    exit
fi

# Perform alignment 
# The only difference between paired end and single end
# is that for single end, fragment length and standard deviation must be given.

if [ "$pairedend" = true ]; then

    align_and_estimate_abundance.pl --transcripts $transcriptome \
    --prep_reference \
    --seqType fq \
    --samples_file "$samples_file" \
    --est_method kallisto \
    --output_dir "$output" \
    --thread_count $thread \
    --gene_trans_map $gene_trans_map

else
    align_and_estimate_abundance.pl --transcripts $transcriptome \
    --seqType fq \
    --samples_file "$samples_file" \
    --est_method kallisto \
    --output_dir "$output" \
    --thread_count $thread \
    --gene_trans_map $gene_trans_map \
    --fragment_length $fragment_length \
    --fragment_std 20 \
    --prep_reference
fi


# Perform conversion to matrix
nbr_of_sample=$( sed -n '$=' "$samples_file" )


# Prepare a variable containing location of all abundance.tsv file for each sample.
# Mandatory for the abundance estimates to matrix script.
for i in $(eval echo "{1..${nbr_of_sample}}");
do
    samp_rep=$( head -n $i "$samples_file" | tail -n 1 | tr -s " " | tr " " "\t"| cut -f 2 )
    list_file=""$list_file" "$samp_rep"/abundance.tsv"
done


abundance_estimates_to_matrix.pl --est_method kallisto \
--gene_trans_map $gene_trans_map \
--out_prefix kallisto \
--name_sample_by_basedir $list_file
[ ! -d "${output}"/matrix/ ] && mkdir "${output}"/matrix/

# Rearrange output folder

# Move output from alignment to the given output in the parameters
# Note that for some  reasons, the output is saved on the current directory
# instead to the output parameters. Moving results at the end is a "good" fix.
mv kallisto* "$output"/matrix/

[ ! -d "${output}"/alignment/ ] && mkdir "${output}"/alignment/

for i in $list_file;
do
    directory=$( dirname $i )
    [ ! -d "${output}"/alignment/${directory} ] && mkdir "${output}"/alignment/${directory}
    mv $directory "$output"/alignment/$directory
done