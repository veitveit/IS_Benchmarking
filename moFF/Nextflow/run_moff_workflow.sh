#!/bin/bash

#time docker build . -t "veitveit/moffworkflow:dev";

time nextflow run ./main.nf -with-report -with-trace -with-timeline -with-dag -profile test,docker

#nextflow run main.nf -profile docker
