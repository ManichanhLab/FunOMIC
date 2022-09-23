# FunOMIC
![FunOMIC logo](https://manichanh.vhir.org/assets/img/funomic_logo.jpeg)
## Introduction
While analysis of the bacterial microbiome has become routine, that of the fungal microbiome is still hampered by the lack of robust databases and bioinformatic 
pipelines. Here, we present FunOMIC, a pipeline with built-in taxonomic (1.6 million marker genes) and functional (3.4 million non-redundant fungal proteins) 
databases for the identification of fungi. Applied to more than 2,600 human metagenomic samples, the tool revealed fungal species associated with geography, 
body sites, and diseases. Correlation network analysis provided new insights into inter-kingdom interactions. With this pipeline and two of the most comprehensive 
fungal databases, we foresee a fast-growing resource for mycobiome studies.
## Usage
1, Required third-party tools:

A full installation of FunOMIC requires Python 3.3+ (modules needed: pandas, math, subprocess, argparse), bowtie2, samtools, flash2, and diamond. 

Please make sure you successfully installed the dependancies.

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[samtools](http://www.htslib.org/)

[flash2](https://github.com/dstreett/FLASH2)

[diamond](https://github.com/bbuchfink/diamond/wiki)

2, Download FunOMIC pipeline:
```
git clone https://github.com/ManichanhLab/FunOMIC.git
```
3, Add the path of FunOMIC folder in your PATH variable, for example:
```
export PATH=$PATH:/home/ManichanhLab/Downloads/FunOMIC
```
You can also add this command line to your ~/.profile file

4, To use FunOMIC pipeline, make sure you have downloaded the required database from https://manichanh.vhir.org/funomic/.
```
FunOMIC.sh -1 read1.fastq -2 read2.fastq -p output_prefix -o output_folder \
-a bacterialDB_path -b FunOMIC-T_folder -c FunOMIC-P_folder -t #_threads
```

