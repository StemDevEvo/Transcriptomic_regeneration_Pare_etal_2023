#!/usr/bin/bash

. ../config.txt

[ ! -d "$output" ] && mkdir "$output"



bash ../src/03_alignment.sh $transcriptome $sample_file "$output" $threads $gene_trans_map $r_length