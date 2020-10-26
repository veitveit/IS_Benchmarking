# Workflow to analyze label-free data with MaxQuant
This workflow is based on Nextflow, running with docker and SDRF implemented. 


## Getting started

This script assumes that Docker is available on your system and is targeting Debian/Ubuntu linux distributions.

This script is working with Nextflow,and needs the SDRF file for the project wanted to be run.

The SDRF can be found under annotated projects, and for the PXD001819, the file is added under data.
URL for SDRF files: https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects

## Run benchmarking data set

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow, given the following commands.

1) The rawfolder, and nextflow needs read and write access to the directory. The path needs to absoulte and not relative.
2) The SDRF.tsv file. Relative paths works fine.
3) The fasta file, and this needs a absolute file path, not relative. 


The SDRF for pxd001819 file is given in the data folder of this workflow:s https://github.com/veitveit/IS_Benchmarking/tree/master/MaxQuant/Nextflow/NextFlow/data

The outputfiles will be written into _RAWFOLDER_ and the copied into the results folder in the base directory


```
nextflow main.nf --raws "rawfile Path without ending "/" " --sdrf "SDRF file path .tsv" --fasta "fasta file Path .fasta"  -profile docker
```


