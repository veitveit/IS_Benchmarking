# Workflow to analyze label-free data with MaxQuant
This workflow is based on a docker image which has to be created beforehand: exaexa:maxquant

## Getting started

This script assumes that Docker is available on your system and is targeting Debian/Ubuntu linux distributions.

This workflow is not so far implemented in Nextflow as it does not use more than one tool

## Run benchmarking data set

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow with the following command and parameters after  
a) changing _RAWFOLDER_ to the folder where the raw files are located.   
b) placing the parameter file _mqpar_LFW_PXD001819.xml_ into in the current folder  
c) changíng the paths of both fasta file and the raw files in the maxquant parameter file. 

The parameter file is given in the Results folder of this workflow: https://github.com/veitveit/IS_Benchmarking/edit/master/MaxQuant/Results

The outputfiles will be written into _RAWFOLDER_


```
docker run -it --rm -v ${PWD}:/import -v RAWFOLDER:/input_raw veitveit/maxquant:dev bash  -c 'cd /import; maxquant  mqpar_LFQ_PXD001819.xml'
```


