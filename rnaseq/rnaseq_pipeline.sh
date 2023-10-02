


#########################################################################################################
##################### Build references: see rnaseq_pipeline_build_references.sh #########################
########## Using GRCh37.p13 genome, and refseq curated annotations (Apr 2018) ###########################
#########################################################################################################





####################################################################################
############################# Run STAR to align ####################################
####################################################################################


### STAR mapping
# use STAR.config, which uses mostly Encode alignment parameters
cd ~/CCA/CCA_RnaSeq/processed/align_STAR
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH1_H69-1_S39_L005_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH1_H69-1_S39_L005_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH2_H69-2_S40_L005_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH2_H69-2_S40_L005_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Egi-1-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH3_Egi-1-1_S41_L005_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH3_Egi-1-1_S41_L005_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Egi-1-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH4_Egi-1-2_S42_L005_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH4_Egi-1-2_S42_L005_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH5_S5-EV-1_S43_L006_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH5_S5-EV-1_S43_L006_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH6_S5-EV-2_S44_L006_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH6_S5-EV-2_S44_L006_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH7_S5-WT-1_S45_L006_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH7_S5-WT-1_S45_L006_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH8_S5-WT-2_S46_L006_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH8_S5-WT-2_S46_L006_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH9_S5-MUT-1_S47_L007_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH9_S5-MUT-1_S47_L007_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH10_S5-MUT-2_S48_L007_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH10_S5-MUT-2_S48_L007_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-KO-EV-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH11_H69_KO_EV-1_S51_L008_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH11_H69_KO_EV-1_S51_L008_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-KO-EV-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH12_H69_KO_EV-2_S52_L008_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH12_H69_KO_EV-2_S52_L008_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-KO-WT-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH13_H69_KO_WT-1_S53_L008_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH13_H69_KO_WT-1_S53_L008_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-KO-WT-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH14_H69_KO_WT-2_S54_L008_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH14_H69_KO_WT-2_S54_L008_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-KO-MUT-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH15_H69_KO_MUT-1_S49_L007_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH15_H69_KO_MUT-1_S49_L007_R2_001.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix H69-KO-MUT-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/JH16_H69_KO_MUT-2_S50_L007_R1_001.fastq.gz ~/CCA/CCA_RnaSeq/raw/JH16_H69_KO_MUT-2_S50_L007_R2_001.fastq.gz


# batch 2
cd ~/CCA/CCA_RnaSeq/processed/align_STAR_batch2
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 17231203_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/17231203T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/17231203T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 2000861_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/2000861T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/2000861T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 29150572_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/29150572T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/29150572T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 33587479_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/33587479T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/33587479T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 72600277_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/72600277T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/72600277T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A035_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/A035T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/A035T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C080_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/C080T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/C080T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C096_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/C096T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/C096T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix HUCCT_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/HUCCT_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/HUCCT_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix KKU100_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/KKU100_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/KKU100_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix M213_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/M213_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/M213_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 260418N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/TBT-CCA-260418_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/TBT-CCA-260418_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 807N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/TBT-CCA-807_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/TBT-CCA-807_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix W39_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/W39T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/W39T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Y140_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/Y140T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/Y140T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Y65_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/Y65T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/Y65T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z3722_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/Z3722T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch2_jan2019/Z3722T_2.fastq.gz


# batch 3
cd ~/CCA/CCA_RnaSeq/processed/align_STAR_batch3
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 1202_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/1202T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/1202T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 2000123_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/2000123T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/2000123T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 2143_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/2143T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/2143T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 23474504_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/23474504T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/23474504T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 3001_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/3001T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/3001T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 30131247_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/30131247T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/30131247T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 45782570_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/45782570T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/45782570T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 824_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/824T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/824T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 86014838_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/86014838T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/86014838T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 90866096_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/90866096T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/90866096T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A153_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/A153T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/A153T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix B083_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/B083T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/B083T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C150_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/C150T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/C150T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 3011118N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/TBTCCA011118N_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/TBTCCA011118N_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 4081118N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/TBTCCA081118N_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/TBTCCA081118N_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z2403_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/Z2403T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/Z2403T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z3585_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/Z3585T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/Z3585T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z639_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/Z639T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch3_mar2019/Z639T_2.fastq.gz


# batch 4
cd ~/CCA/CCA_RnaSeq/processed/align_STAR_batch4

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 001G_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/001GT_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/001GT_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 21914187_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/21914187T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/21914187T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 31215286_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/31215286T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/31215286T_2.fastq.gz

# note that A003T is actually 1965T (see jingyi email)
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 1965_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/1965T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/1965T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A042_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A042T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A042T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A059_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A059T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A059T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A096_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A096T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A096T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A100_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A100T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/A100T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix R149_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/R149T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/R149T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix W40_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/W40T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/W40T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Y74_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/Y74T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch4_apr2019/Y74T_2.fastq.gz


# batch 5
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 26814273_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/26814273_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/26814273_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 61271150_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/61271150_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/61271150_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_684_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/684_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/684_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 93199056_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/93199056_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/93199056_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A074_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A074_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A074_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A142_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A142_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A142_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A157_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A157_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A157_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix A169_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A169_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/A169_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix B085_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/B085_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/B085_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C008_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C008_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C008_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C078_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C078_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C078_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C105_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C105_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C105_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C144_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C144_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C144_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix C176_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C176_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/C176_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix D017_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/D017_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/D017_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_33_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN33_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN33_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_39_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN39_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN39_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_43_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN43_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN43_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_46_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN46_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/TAIWAN46_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Y002_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Y002_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Y002_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Y091_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Y091_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Y091_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z12244N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z12244N_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z12244N_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z12265N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z12265N_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z12265N_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z12267N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z12267N_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z12267N_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z321_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z321_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z321_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z508_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z508_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch5_apr2019/Z508_2.fastq.gz


# batch 6
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-CON-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVCONR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVCONR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-CON-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVCONR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVCONR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-OLA-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVOLAR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVOLAR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-OLA-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVOLAR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVOLAR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-SP-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVSPR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVSPR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-EV-SP-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVSPR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5EVSPR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-CON-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTCONR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTCONR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-CON-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTCONR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTCONR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-OLA-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTOLAR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTOLAR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-OLA-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTOLAR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTOLAR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-SP-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTSPR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTSPR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-MUT-SP-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTSPR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5MUTSPR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-CON-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTCONR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTCONR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-CON-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTCONR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTCONR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-OLA-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTOLAR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTOLAR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-OLA-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTOLAR2_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTOLAR2_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-SP-1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTSPR1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTSPR1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5-WT-SP-2_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTSPR2R_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch6_apr2019_ribozerogold/S5WTSPR2R_2.fastq.gz


# batch 7
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 03B5663_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/03B5663T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/03B5663T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 09B02141_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/09B02141T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/09B02141T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 19286759_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/19286759T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/19286759T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 21914187N_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/21914187N_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/21914187N_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 39707259_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/39707259T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/39707259T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 77071507_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/77071507_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/77071507_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 97203223_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/97203223T_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/97203223T_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix B011_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/B011_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/B011_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_14_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/TAIWAN14_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/TAIWAN14_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TAIWAN_38_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/TAIWAN38_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/TAIWAN38_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z2778_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/Z2778_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch7_jun2019/Z2778_2.fastq.gz


# batch 8
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix 2KU452_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/2KU452_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/2KU452_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix M055_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/M055_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/M055_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix M156_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/M156_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/M156_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix MEC_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/MEC_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/MEC_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix SNU1196_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/SNU1196_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/SNU1196_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix SNU245_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/SNU245_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/SNU245_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix S5_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/TBTCCAS5_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/TBTCCAS5_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix TFK1_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/TFK1_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/TFK1_2.fastq.gz

qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --outFileNamePrefix Z5606_ --parametersFiles STAR.config --readFilesIn ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/Z5606_1.fastq.gz ~/CCA/CCA_RnaSeq/raw/batch8_aug2019/Z5606_2.fastq.gz


##################################################
#### Check mapping rate to mitochondria genome (check mapped.chrM.numreads / mapped.numreads)
for i in *_Aligned.sortedByCoord.out.bam
do
   sampname=$(basename _Aligned.sortedByCoord.out.bam)
   qsub samtools index $i
done

for i in *_Aligned.sortedByCoord.out.bam
do
   sampname=$(basename _Aligned.sortedByCoord.out.bam)
   qsub -o ${i}.mapped.numreads samtools view -F 4 -c $i
done

for i in *_Aligned.sortedByCoord.out.bam
do
   sampname=$(basename _Aligned.sortedByCoord.out.bam)
   qsub -o ${i}.mapped.chrM.numreads samtools view -F 4 -c $i chrM
done






####################################################################################
########################### RSEM quantification ####################################
####################################################################################

# check bam format
cd ~/CCA/CCA_RnaSeq/processed/quantify_RSEM/batch8
for i in ../../align_STAR_batch8/*_Aligned.toTranscriptome.out.bam
do
  qsub rsem-sam-validator $i
done



# quantify with RSEM
cd ~/CCA/CCA_RnaSeq/processed/quantify_RSEM/batch8
for i in ../../align_STAR_batch8/*_Aligned.toTranscriptome.out.bam
do
  sample_id=$(basename $i _Aligned.toTranscriptome.out.bam)
  echo $sample_id
  echo $i
  qsub -pe smp 8 rsem-calculate-expression -p 8 --bam --paired-end --no-bam-output --seed 123456789 --forward-prob 0.5 $i ~/CCA/genome/RSEM_GRCh37.p13_refseq/GRCh37.p13_refseq $sample_id
done










####################################################################################
########################### run PicardMetrics ######################################
####################################################################################
cd ~/CCA/CCA_RnaSeq/processed/picardmetrics
for i in ../align_STAR_batch2/*_Aligned.sortedByCoord.out.bam
do
  qsub picardmetrics run -f ~/picardmetrics/picardmetrics.conf -k -r -o ~/CCA/CCA_RnaSeq/processed/picardmetrics $i
done
