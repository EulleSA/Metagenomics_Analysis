# Cactus: A pipeline for identification and visualization of microbial consortium interactions in metabolic pathways using shotgun metagenoma sequences.

# TOPIC 1
- Install Python3: https://www.python.org/downloads/release/python-391/

- Install snakemake according with documentation: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

- Download the repository using the program "git". If you don't have the program "git" installed in your computer, access the following link: https://git-scm.com/downloads

                                git clone https://github.com/EulleSA/Metagenomics_Analysis.git.

# TOPIC 2
After, go to the folder of the application downloaded and create the following folders:
        mkdir ...
Put the metagenome file you want to analyze in the folder **"data_raw/"** and execute the follow command line in the directory of the **Snakefile**:

                                                    snakemake -p --cores

# NOTES

The pipeline needs databases to perform the functional and taxonomic annotations. If you already have BASTA and its respective databases, in addition to NCBI-nr on your machine, you just create a symbolic link to the respective folders: 
