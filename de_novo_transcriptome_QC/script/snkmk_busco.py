#!/usr/bin/env python

#########
# Launch busco with the chosen parameters for the snakemake pipeline. 
# Written by Louis PARE, mail: louis.pare@ijm.fr
#########


#Imports :
import os
import os.path
from os import path
import sys
import subprocess

# main

otp_busco = sys.argv[1]
transcriptome = sys.argv[2]

# Select database type.
odb = otp_busco.split("/")[-1]

# Swiching current directory because busco can not take full path as output.
cwd = os.getcwd()
os.chdir(otp_busco)

#Launch busco
subprocess.call("busco -i "+transcriptome+" -l "+odb+"_odb10 -o "+odb+" -m transcriptome", shell = True)

os.chdir(cwd)
