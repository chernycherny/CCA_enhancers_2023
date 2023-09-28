This directory contains:

process_chipseq_enhancers.sh: shell script containing commands we used to process the H3K27ac chip-seq data
generate_UCSC_trackshub.R: R script containing commands to generate track hubs for UCSC genome browser
analysis_cca_enhancers_0_preliminary.Rmd: R script containing commands to analyze consensus peaks, recount peaks, plot heatmap and PAC


Both files are not meant to be run from start to finish. It won't work because:
1. The directory doesn't include the various required data files (please download and process the data files as required, see Data Availability in the paper)
2. The scripts are cobbled together from various code written in the course of this project

Instead, the code is provided only as a reference / documentation of how we performed the analysis in this project. You may refer to this code should you want to write your own program to reproduce our results!
