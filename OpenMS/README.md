Instructions:

* Install Java 8+
* Install docker or singularity
* Install nextflow (`wget -qO- https://get.nextflow.io | bash`)
* ~/nextflow run nf-core/proteomicslfq -r dev -profile docker --sdrf https://raw.githubusercontent.com/bigbio/proteomics-metadata-standard/master/annotated-projects/PXD001819/sdrf.tsv --database https://raw.githubusercontent.com/nf-core/test-datasets/proteomicslfq/testdata-aws/uniprot\_yeast\_reviewed\_isoforms\_ups1\_crap.fasta\_td.fasta 
