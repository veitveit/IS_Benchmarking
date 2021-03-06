## ELIXIR Implementation Study "Comparison, benchmarking and dissemination of proteomics data analysis pipelines"
See also https://elixir-europe.org/about-us/commissioned-services/proteomics-pipelines

## IMPORTANT NOTE:
__We now moved to https://github.com/orgs/wombat-p/ __


The aim of this study is to create different pipelines for the analysis of label-free proteomics data and systematically compare their results. If possible, the workflows are implemented in [NextFlow](https://www.nextflow.io/) to allow reproducible results and scalable execution on different computer architectures.

### File structure of repository
For examples and simple templates containing e.g. a workflow for file conversion (Thermo raw -> mgf), see [Examples_and_templates](https://github.com/veitveit/IS_Benchmarking/tree/master/Examples_and_templates)

Each of the different workflows is represented by a folder containing the subfolders "Nextflow" and "Results"

#### Content of Nextflow folder
To allow fully reproducible results, each workflow should be executable using the "docker" option in Nextflow 
- [ ] __Nextflow script__ Full set of routines to run the workflow
- [ ] __Nextflow configuration file__ Default options and basic profiles to run the workflow
- [ ] __environment.yml__ Conda packages all software available via (bio)conda
- [ ] __Dockerfile__ Installations instructions for all tools with automatic inclusion of conda packages
- [ ] __data__ Folder with simple example data set to test the workflow 
- [ ] __conf__ Folder with configuration files for testing and execution
- [ ] __README.md__ Detailed description of parameter files and how to run the workflow on our data set.  Include details about source of data and database. Please also include other options (e.g. command line options/flags) necessary to run the workflow

#### Content of Result folder
- [ ] __Result files__ from each lab for identified and quantified peptides and proteins, preferably in the form of a table. Please change the name of the files to _WorkflowName_Lab_originalfilename_
- [ ] __Parameter files__ for running the workflow on the data set
- [ ] __Performance report__ (TODO: details) E.g. execution time on defined architecture, ...
