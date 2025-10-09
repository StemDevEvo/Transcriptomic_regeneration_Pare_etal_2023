#!/usr/bin/env python

#########
# Launch bowtie2 with the chosen parameters for the snakemake pipeline. 
# Done in two steps: first an index is created then each sample is aligned with the transcriptome.
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
if len(sys.argv) != 5:
    exit(1)

otp_cdhit = sys.argv[1]
otp_bowtie = sys.argv[2]
proc = sys.argv[3]
sample = sys.argv[4]




# index creation.
subprocess.call("bowtie2-build " + otp_cdhit + " " + otp_bowtie + "/bowtie_build.fasta ", shell = True)
# align each paired sample with transcriptome.
list_sample = get_list_sample(sample)
for i in list_sample:
    sample_name = i[0].split("/")[-1]
    subprocess.call("bowtie2 -p " + proc + " -q --no-unal -k 20 -x " + otp_bowtie + "/bowtie_build.fasta -1 " + i[0] + " -2 " + i[1] + " 2>" + otp_bowtie + sample_name + ".txt| samtools view -@10 -Sb -o bowtie2.bam", shell = True)
