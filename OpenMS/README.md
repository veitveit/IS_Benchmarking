Info about used software (e.g. environment.yml, Dockerfile): https://github.com/nf-core/proteomicslfq

Instructions:

* Install Java 8+
* Install docker or singularity
* Install nextflow (`wget -qO- https://get.nextflow.io | bash`)
* ~/nextflow run nf-core/proteomicslfq -r dev -profile docker --sdrf https://raw.githubusercontent.com/bigbio/proteomics-metadata-standard/master/annotated-projects/PXD001819/sdrf.tsv --database https://raw.githubusercontent.com/nf-core/test-datasets/proteomicslfq/testdata-aws/uniprot\_yeast\_reviewed\_isoforms\_ups1\_crap.fasta\_td.fasta

Results:

* MSstats folder contains final results/quants post-processed by MSstats
* openms\_msstats-ready\_out.csv contains the raw aggregated intensities
