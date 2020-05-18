#!/bin/bash

time docker build . -t "veitveit/moffworkflow:dev";

nextflow  run main.nf --raws 'data/*.raw' --txts 'data/*.txt' --inifile "data/configuration_iRT.ini" -profile docker 
