# DENV_MinION_Script
Shell script used to run the ARTIC minion pipeline to analyze amplicon sequencing data generated using the pan-DENV primers sequenced on ONT platforms

This pipeline takes raw read data in the form of fastq files generated from Oxford Nanopore Technologies platforms, determines the most likely Dengue virus serotype by mapping initial reads to multiple reference genomes, then utilizes the ARTIC minion pipeline to align reads to appropriate reference genome, soft clip primer sequences, and generates a consensus sequence genome. 

This program will output multiple file types, include consensus sequence genomes in the form of a fasta file, a vcf file displaying where mutations occur in relation to the aligned reference genome, and multiple BAM files. 
To run this software, you must have the access to the ARTIC minion pipeline and all underlying programs. Visit their GitHub and/or readthedocs page to learn more about the ARTIC minion pipeline. 

ARTIC Network field-bioinformatics pipeline GitHub: https://github.com/artic-network/fieldbioinformatics

ARTIC Network read the docs: https://artic.readthedocs.io/en/latest/commands/

Running the script
This script is designed to be run in a UNIX environment. It requires users to set hard paths to an output directory, a directory containing the primer scheme (see the ARTIC minion documents for more information), a directory containing fastq files for each sample to be analyzed in subdirectories, and to a list file that contains the barcode number and sample name for each sample to be analyzed. 

Other options to edit are the minimum and maximum read length for filtering and the Medaka model to use. 
