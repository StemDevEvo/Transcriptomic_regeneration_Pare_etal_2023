#!/usr/bin/env python

#########
# launch transrate with the chosen parameter for the snakemake pipeline.
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
if len(sys.argv) != 5:
    exit(1)

otp_cdhit = sys.argv[1]
otp_transrate = sys.argv[2]
proc = sys.argv[3]
sample = sys.argv[4]

# get trimmed samples name.
list_sample = get_list_sample(sample)
left_reads = ""
right_reads = ""
for i in list_sample:
    left_reads = left_reads + i[0] + ","
    right_reads = right_reads + i[1] + ","
left_reads = left_reads[:-1]
right_reads = right_reads[:-1]

# launch transrate.
subprocess.call("transrate --assembly " + otp_cdhit + " --left "+ left_reads +" --right " + right_reads + " --threads " + proc + " --output " + otp_transrate, shell = True)
