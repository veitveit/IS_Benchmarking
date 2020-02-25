# Workflow to analyze label-free data with the searchgui and proline
This workflow is based on a docker image which has to be created beforehand. A Dockerfile is available. 
The image needs to be named veitveit/prolineworkflow:dev

# Getting started

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
