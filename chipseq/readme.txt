This directory contains:

process_chipseq_enhancers.sh: shell script containing commands we used to process the H3K27ac chip-seq data, peak calling, generating signal tracks, and motif finding
generate_UCSC_trackshub.R: R script containing commands to generate track hubs for UCSC genome browser
analysis_cca_enhancers_0_preliminary.Rmd: R script containing commands to analyze consensus peaks, recount peaks, plot heatmap and PAC
analysis_cca_enhancers_0b_unsup_clustering.Rmd: R script containing commands for unsupervised clustering and consensus clustering 
analysis_cca_enhancers_1_targetmapping.R: R script containing commands to map enhancers to target genes
analysis_cca_enhancers_1_targetmapping.Rmd: R script containing commands to analyze enhancer-gene mapping
analysis_cca_enhancers_2b.2_diffbind_new.Rmd: R script containing commands to setup and run differential binding analyses
analysis_cca_enhancers_2h_other_analyses.Rmd: R script containing commands for motif analysis with Homer
analysis_cca_enhancers_3b_somatic_enhancers_groupspecific.Rmd: R script containing commands to analyze group-specific somatic enhancers & genesets
analysis_cca_enhancers_3c_estrogen.Rmd: R script containing commands to analyze the ESTRO enhancers (Figure 2)
analysis_cca_enhancers_3d_oxphos.Rmd: R script containing commands to analyze the OXPHO enhancers (Figure 3)
analysis_cca_enhancers_3e_immune.Rmd: R script containing commands to analyze the IMMUN enhancers (Figure 4)
analysis_cca_enhancers_functions.R: R script containing functions used by the above scripts 


The source codes are not meant to be run from start to finish. It won't work because:
1. The data files are not included here (please download and process the data files as required, see Data Availability in the paper)
2. The source codes are cobbled together from various codes written in the course of this project
3. The source codes are far from "production-ready" code. They are at R&D stage and are written ad-hoc to accomplish research objectives.
Instead, the source codes are provided only as reference / documentation of how we performed the analysis in this project. You may refer to this code should you want to write your own program to reproduce our results!
