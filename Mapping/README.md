# Mapping of the RNA-seq data.

## Description

This pipeline will align and create count matrices from RNAseq data using a transcriptome as a reference.
It will check quality of the reads, trim them and finally align and create count matrices for the given samples.

## What are the steps of the pipeline?

The pipeline is divided in 3 mains steps.
Firstly, the reads are quality checked using [fastqc](https://github.com/s-andrews/FastQC), which generates reports where the quality of the reads can be easily estimated based on several metrics.
Then, the pipeline use [fastp](https://github.com/OpenGene/fastp) in order to trim the reads.
Finally, the pipeline use the [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) implementation of the transcript quantification pipeline in order to generate count matrices.

## How to make it work?

The pipeline use a [conda](https://www.anaconda.com/products/distribution) environment to gather all the dependencies needed. 
After the installation of conda, the environment can be created using the environment.yml file using the following command: 

`conda env create -f environment.yml` 

Then the environment can be loaded using: 

`conda activate NAME_OF_THE_ENVIRONMENT` 


Once the conda environment is ready the configuration file must be completed. It is a text file where parameters can be written. 
Each parameters are in the form of: 

`name_of_the_parameter=value_of_the_parameter` 

Their must not be any spaces between the value of the parameters and the = sign. 

Once the parameter file has been completed, each step of the pipeline can be launched. In order to launch a step, the current directory must be changed to the `launch` folder, then each launch file can be started using the command: 

`bash name_of_the_launch_file.sh` 

Before launching the next step, make sure that the current step has finished. For exemple, if you launch the 01_launch_fastqc.sh step, which will perform quality control on the reads, wait for all the output to be written, then you can perform the second step of the pipeline.