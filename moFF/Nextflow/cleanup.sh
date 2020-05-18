#!/bin/bash

find "./" -maxdepth 1 \( \
	-name "trace.txt*" \
		-o \
	-name "timeline.html*" \
		-o \
	-name "timeline.html*" \
		-o \
	-name "report.html*" \
		-o \
	-name ".nextflow.log*" \
		-o \
	-name "dag.dot*" \
\) -delete;
