#!/bin/bash

time docker build . -t "veitveit/moffworkflow:dev";

time nextflow run main.nf \
	--raws "absence_peak_data/" \
	--raw_repo "absence_peak_data/" \
	--inifile "absence_peak_data/configuration_iRT.ini" \
	--loc_out "./absence_peak_data/results/" \
	-profile docker \
		-with-report -with-trace -with-timeline -with-dag;