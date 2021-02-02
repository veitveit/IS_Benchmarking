# Workflow to analyze label-free data with MaxQuant
This workflow is based on Nextflow, running with singularity and SDRF implemented. Normalization and statistical comparisons using NormalyzerDE are conducted on the MaxQuant results.


## Getting started

This script assumes that Singularity is available on your system and is targeting Debian/Ubuntu linux distributions.

This script is working with Nextflow,and needs the SDRF file for the project wanted to be run.

The SDRF can be found under annotated projects, and for the PXD001819, the file is added under data.
URL for SDRF files: https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects

An experimental design file for the NormalyzerDE part can also be found in the data folder. 

## Run benchmarking data set

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow, giving the following parameters:

1) The rawfolder, and nextflow needs read and write access to the directory. The path needs to absolute and not relative.
2) The SDRF.tsv file. Relative paths works fine.
3) The fasta file, and this needs a absolute file path, not relative. 
4) The absolute path to the Normalyzer experimental design file
5) Which normalization method to use.
5) The group comparisons to perform in NormalyzerDE.

These settings and more for the PXD001819 project are given in the configuration file PXD001819.conf under data.
Just make sure to update the paths in the configuration file, and then run as

```
nextflow main.nf -c PXD001819.conf -profile singularity
```

The SDRF for PXD001819 file is given in the data folder of this workflow:s https://github.com/veitveit/IS_Benchmarking/tree/master/MaxQuant/Nextflow/NextFlow/data

The output files will be written into the combined directory of the specified MaxQuant raw folder and then copied into the specified output folder

The NormalyzerDE results will be copied into a folder named as the 'project' parameter

As alternative to running with a comnfiguration file, parameters can be set at the command line when running NextFlow:

```
nextflow main.nf --raws "rawfile Path without ending "/" " --sdrf "SDRF file path .tsv" --fasta "fasta file Path .fasta"  -profile singularity
```


