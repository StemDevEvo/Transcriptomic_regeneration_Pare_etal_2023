#!/usr/bin/env python

#########
# Launch Ex90N50 statistic computation with the chosen parameters for the snakemake pipeline. 
# Written by Louis PARE, mail: louis.pare@ijm.fr
#########


#Imports :
import os
import os.path
from os import path
import sys
import subprocess

# Functions :

# inputs : line = string containing words separated by spaces, nbr_word = which word to be returned, start at 0 to n-1 where n is number of words. 
# output : chosen word.
# return the chosen word from a line, words must be separated by spaces.
def get_word(line, nbr_word):
    return line.split()[nbr_word]


# input : file path to a sample.
# output : name of the sample.
# used to return the name of the sample.  
def trim_name(sample):
    return sample.split("/")[-1]

# input : path = path to the file containing information on samples.
# output : list containing path of each samples files. 
# Used to get the names of the paired sample from the sample file. return a 2 dimentional list with first dimension each replicate and second dimension the paired end samples.
def get_list_sample(path):
    #create empty list.
    list_sample = []
    # open file and parse line by line
    with open(path,"r+") as file:
        for line in file:
            list_sample.append((get_word(line,2),get_word(line,3)))
    return list_sample


#checking arg
if len(sys.argv) != 6:
    print("One argument is missing.")
    exit(1)

otp_cdhit = sys.argv[1]
otp_quant = sys.argv[2]
proc = sys.argv[3]
sample = sys.argv[4]
file_type = sys.argv[5]


# Get samples names and modify it to match the trim files.
cwd = os.getcwd()
os.chdir(otp_quant)
list_sample = get_list_sample(sample)
left_reads = []
right_reads = []
for i in list_sample:
    left_reads.append(i[0])
    right_reads.append(i[1])

# Sample files creation for the trimmed samples.
lines = []
with open(sample,"r+") as file:
    for line in file:
        lines.append(line)
with open(otp_quant + "/sample_file.txt","w") as file:
    for i in range(len(lines)):
        file.write(get_word(lines[i],0) + "\t" + get_word(lines[i],1) + "\t" + left_reads[i] + "\t" + right_reads[i])
        if i != len(lines)-1:
            file.write("\n")
with open(otp_quant + "/abundance_estimates_samples.txt","w") as file:
    for i in range(len(lines)):
        file.write(otp_quant +"/"+ get_word(lines[i],1) + "/RSEM.isoforms.results")
        if i != len(lines)-1:
            file.write("\n")

    # Launch RSEM.
subprocess.call("$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts " + otp_cdhit + " --seqType " + file_type + " --samples_file " + otp_quant + "/sample_file.txt" + " --est_method RSEM --output_dir " + otp_quant + " --aln_method bowtie2 --prep_reference --thread_count " + proc, shell = True)
subprocess.call("$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --gene_trans_map none --est_method RSEM --quant_files " + otp_quant + "/abundance_estimates_samples.txt --name_sample_by_basedir", shell = True)
subprocess.call("$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl " + otp_quant + "/RSEM.isoform.TMM.EXPR.matrix " + otp_cdhit +" | tee ExN50.stats", shell = True)
os.chdir(cwd)