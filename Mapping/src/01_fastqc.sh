output_directory=$1
path_to_reads=$2

all_reads=$(ls ${path_to_reads}/*.gz)

fastqc -t 6 -o $output_directory -f fastq $all_reads