# Workflow to analyze label-free data with MaxQuant
This workflow is based on Nextflow, running with SDRF implemented. Normalization and statistical comparisons using NormalyzerDE are conducted on the MaxQuant results.


## Getting started

This script is working with Nextflow,and needs the SDRF file for the project wanted to be run, as well as a file 
for the experimental design

The SDRF can be found under annotated projects, and for the PXD001819, the file is added under data.
URL for SDRF files: https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects

An experimental design file for the NormalyzerDE part can also be found in the data folder. 

## Run benchmarking data set

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow, giving the following parameters:

1) The rawfolder, and nextflow needs read and write access to the directory. The path needs to absolute and not relative.
2) The SDRF.tsv file. Relative paths works fine.
3) The fasta file, and this needs a absolute file path, not relative. 
4) The absolute path to the experimental design file
5) (Optional) Which normalization method to use.
5) (Optional) The group comparisons to perform in NormalyzerDE. Without this parameter, the comparison will be versus 
the first condition

Just make sure to update the paths in the configuration file, and then run as

```
nextflow run main.nf --raws 'RAW FILES OR FOLDER WITH RAW FILES' --fasta 'SPECIFYFOLDER/yeast_UPS.fasta' --sdrf data/sdrf_UPS.tsv --experiment_design 'SPECIFYFOLDER/pxd001819.txt' 
```

Still missing: paramter for tolerances, enzymes and PTMs

