# Transcriptomic_regeneration_Pare_etal_2023

This repository contain the code used to generate the analysis from Pare et al., 2023.

**Transcriptomic landscape of posterior regeneration in the annelid Platynereis dumerilii**
Louis Paré, Loïc Bideau, Loeiza baduel, Caroline Dalle, Médine Benchouaia, Stephan Q.Schneider, Lucie Laplane, Yves Clément, Michel Vervoort & Eve Gazave
BMC Genomics 2023; doi: https://doi.org/10.1186/s12864-023-09602-z 

The raw fastq files can be obtained on the ENA repository linked to the paper, at the following adress: https://www.ebi.ac.uk/ena/browser/view/PRJEB63219

It is composed of three parts, the `de_novo_transcriptome_QC` folder contain a snakemake pipeline to compute quality control of de novo generated transcriptomes.

The `Mapping` folder contain the script to obtain the count tables from the Platynereis reference transcriptome

The thrid folder `DEG_analysis` contains the code and results for the differentially expressed genes analysis.