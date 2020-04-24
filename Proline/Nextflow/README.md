# Workflow to analyze label-free data with the searchgui and proline
This workflow is based on a docker image which has to be created beforehand. A Dockerfile is available. 
The image needs to be named veitveit/prolineworkflow:dev

## Getting started

Here is a little script that can be used to test the workflow execution.
This script assumes that Docker is available on your system and is targeting Debian/Ubuntu linux distributions.

```
# Install Java
sudo apt-get install openjdk-8-jdk

# Fix docker permission issue (https://stackoverflow.com/questions/48957195/how-to-fix-docker-got-permission-denied-issue)
# Note: this allows to run docker as non-root user
sudo usermod -aG docker $USER
newgrp docker

# Fetch scripts and data
git clone https://github.com/veitveit/IS_Benchmarking.git
cd ./IS_Benchmarking/Proline/Nextflow/data
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/UPS1_500amol_R1.raw

# Nextflow
cd ..
curl -s https://get.nextflow.io | bash

# Worfklow
./nextflow run main.nf -profile docker,test

```

## Run benchmarking data set

You should have successfully tested the workflow using the procedure above.

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow with the following command and parameters after changing _RAWFOLDER_ to the folder where the raw files are located. You will also need to place the files _pxd001819.txt_, _yeast_UPS.fasta_ and _lfq_param_file_pxd001819.txt_ into in the current folder. These files are given in the Results folder of this workflow: https://github.com/veitveit/IS_Benchmarking/tree/master/Proline/Results

Also adjust the parameter values _max_cpus_ and _max_memory_ to the computing power you have available.
```
nextflow run main.nf --raws 'RAWFOLDER/*.raw' --fasta yeast_UPS.fasta --precursor_mass_tolerance 5 --fragment_mass_tolerance 0.8 --miscleavages 2 \
--variable_mods 'Oxidation of M, Acetylation of protein N-term' --experiment_design 'pxd001819.txt' --lfq_param \
'lfq_param_file.txt' --max_cpus 8 --max_memory 8GB -profile docker
```


