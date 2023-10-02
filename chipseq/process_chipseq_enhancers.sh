

mkdir preprocessed

### QC reads
mkdir ~/CCA/CCA_ChipSeq/processed/preprocessed/fastqc
cd ~/CCA/CCA_ChipSeq/chipseq_raw/batch9_feb2020
for i in *.fastq.gz
do
  qsub ~/fastqc/FastQC/fastqc $i --outdir=/home/gmsv1178/CCA/CCA_ChipSeq/processed/preprocessed/fastqc/
done



#### if FASTQC reports adapter sequences (so far only batch 3 has), trim adapters first
# why not use trimmomatic with quality & adapter trimming for all samples? because chipseq pipeline assumes all reads start from same position in the fragment (when extending reads)
mkdir preprocessed/trim_adapter
for i in ../chipseq_raw/batch3_sep2018/*.fastq.gz
do
  basefilename=$(basename $i)
  qsub -pe smp 8 java -jar ~/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 8 -phred33 $i preprocessed/trim_adapter/$basefilename ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:50
done

# paired end data with adapter (batch 8 & 9). Remove adapter, then just take read 1
# Note: The fragment sizes seem to be around 160 (just take revcomp of read2, overlap with read 1)!!!
# That means many reads will have adapter readthrough --> Trimmomatic discards read2. Reads with no readthrough --> no adapter content --> not trimmed.

mkdir ~/CCA/CCA_ChipSeq/processed/preprocessed/trim_adapter
cd ~/CCA/CCA_ChipSeq/processed
for i in ../chipseq_raw/batch9_feb2020/*_1.fq.gz
do
  samplename=$(basename $i _1.fq.gz)
  qsub -pe smp 8 java -jar ~/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 8 -phred33 -basein $i -baseout preprocessed/trim_adapter/${samplename}.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 MINLEN:50
done
# combine paired and unpaired read 1
cd preprocessed/trim_adapter
rm *_2P.fq.gz
rm *_2U.fq.gz
for i in *_1P.fq.gz
do
  samplename=$(basename $i _1P.fq.gz)
  qsub gunzip ${samplename}_1P.fq.gz
  qsub gunzip ${samplename}_1U.fq.gz
done
for i in *_1P.fq
do
  samplename=$(basename $i _1P.fq)
  mv ${samplename}_1P.fq ${samplename}_1PU.adaptertrim.fastq
  cat ${samplename}_1U.fq >> ${samplename}_1PU.adaptertrim.fastq
done
for i in *_1PU.adaptertrim.fastq
do
  qsub gzip $i
done



#### trim front and back 10bp (-Q33 just specifies illumina format score).
# Preferable to trim fixed amount rather adaptive trimming, because chipseq pipeline assumes all reads start from same position in the fragment, when extending reads
# gzip -c -d $1 | fastx_trimmer -Q33 -f 11 | fastx_trimmer -Q33 -t 10 | gzip --fast -c > $2
mkdir ~/CCA/CCA_ChipSeq/processed/preprocessed/trim
cd ~/CCA/CCA_ChipSeq/chipseq_raw/batch4_oct2018
for i in *.fastq.gz
do
  basefilename=$(basename $i .fastq.gz)
  qsub ~/CCA/CCA_ChipSeq/processed/chipseq_pipeline_1_trim.sh $i ~/CCA/CCA_ChipSeq/processed/preprocessed/trim/$basefilename.trim.fastq.gz
done


# chipseq_pipeline_1_trim.sh code:
gzip -c -d $1 | fastx_trimmer -Q33 -f 11 | fastx_trimmer -Q33 -t 10 | gzip --fast -c > $2



# count reads
cd ~/CCA/CCA_ChipSeq/processed/preprocessed/trim
for i in *.trim.fastq.gz
do
  zcat $i | grep "^@" | wc -l > $i.numreads
done

for i in {171..194}
do
  echo CCA${i} >> ../numreads_1_trim.txt
  cat CCA${i}*.numreads | awk '{sum+= $1}; END {print sum}' >> ../numreads_1_trim.txt
done

# for signomax batch (filenames don't start with CCA)
for i in *.numreads
do
  echo $i >> ../numreads_1_trim.txt
  cat $i | awk '{sum+= $1}; END {print sum}' >> ../numreads_1_trim.txt
done




# QC trimmed reads
cd ~/CCA/CCA_ChipSeq/processed/preprocessed/trim
mkdir ~/CCA/CCA_ChipSeq/processed/preprocessed/fastqc_trimmed
for i in *.trim.fastq.gz
do
  qsub ~/fastqc/FastQC/fastqc $i --outdir=../fastqc_trimmed/
done



# build bwa index
cd ~/genome
qsub ~/bwa-0.7.10/bwa index hs37d5.chr.fa


# align
cd ~/CCA/CCA_ChipSeq/processed
mkdir preprocessed/align
for i in ./preprocessed/trim/CCA00*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA01*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA02*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA03*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA04*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA05*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA06*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA07*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA08*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA09*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done

for i in ./preprocessed/trim/CCA10*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done


for i in ./preprocessed/trim/CCA*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/genome/hs37d5.chr.fa $i
done


for i in ./preprocessed/trim/*.trim.fastq.gz
do
  basefilename=$(basename $i .trim.fastq.gz)
  qsub -pe smp 8 -o preprocessed/align/$basefilename.sam ~/bwa-0.7.10/bwa mem -M -t 8 ~/CCA/genome/hs37d5.chr.fa $i
done







# convert to bam
cd ~/CCA/CCA_ChipSeq/processed
mkdir preprocessed/bam
for i in ./preprocessed/align/*.sam
do
  basefilename=$(basename $i .sam)
  qsub -o preprocessed/bam/$basefilename.bam samtools view -bS $i
done


# sort bam
mkdir preprocessed/sort
for i in ./preprocessed/bam/*.bam
do
  basefilename=$(basename $i .bam)
  qsub samtools sort -o preprocessed/sort/$basefilename.sorted.bam $i
done


# merge multiple bams
mkdir preprocessed/merge
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187T_K4me3.bam ./preprocessed/sort/CCA001*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187T_K4me1.bam ./preprocessed/sort/CCA002*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187T_K27ac.bam ./preprocessed/sort/CCA003*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187T_input.bam ./preprocessed/sort/CCA004*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187N_K4me3.bam ./preprocessed/sort/CCA005*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187N_K4me1.bam ./preprocessed/sort/CCA006*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187N_K27ac.bam ./preprocessed/sort/CCA007*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/21914187N_input.bam ./preprocessed/sort/CCA008*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572T_K4me3.bam ./preprocessed/sort/CCA009*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572T_K4me1.bam ./preprocessed/sort/CCA010*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572T_K27ac.bam ./preprocessed/sort/CCA011*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572T_input.bam ./preprocessed/sort/CCA012*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572N_K4me3.bam ./preprocessed/sort/CCA013*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572N_K4me1.bam ./preprocessed/sort/CCA014*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572N_K27ac.bam ./preprocessed/sort/CCA015*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/29150572N_input.bam ./preprocessed/sort/CCA016*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100T_K4me3.bam ./preprocessed/sort/CCA017*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100T_K4me1.bam ./preprocessed/sort/CCA018*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100T_K27ac.bam ./preprocessed/sort/CCA019*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100T_input.bam ./preprocessed/sort/CCA020*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100N_K4me3.bam ./preprocessed/sort/CCA021*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100N_K4me1.bam ./preprocessed/sort/CCA022*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100N_K27ac.bam ./preprocessed/sort/CCA023*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/A100N_input.bam ./preprocessed/sort/CCA024*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585T_K4me3.bam ./preprocessed/sort/CCA025*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585T_K4me1.bam ./preprocessed/sort/CCA026*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585T_K27ac.bam ./preprocessed/sort/CCA027*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585T_input.bam ./preprocessed/sort/CCA028*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585N_K4me3.bam ./preprocessed/sort/CCA029*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585N_K4me1.bam ./preprocessed/sort/CCA030*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585N_K27ac.bam ./preprocessed/sort/CCA031*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Z3585N_input.bam ./preprocessed/sort/CCA032*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_BAP1.bam ./preprocessed/sort/CCA033*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_inputTF.bam ./preprocessed/sort/CCA034*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/TBTS5_EV_BAP1.bam ./preprocessed/sort/CCA035*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/TBTS5_EV_inputTF.bam ./preprocessed/sort/CCA036*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/TBTS5_WT_BAP1.bam ./preprocessed/sort/CCA037*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/TBTS5_WT_inputTF.bam ./preprocessed/sort/CCA038*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/TBTS5_MUT_BAP1.bam ./preprocessed/sort/CCA039*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/TBTS5_MUT_inputTF.bam ./preprocessed/sort/CCA040*.sorted.bam

qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_K27ac.bam ./preprocessed/sort/CCA041*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_input.bam ./preprocessed/sort/CCA042*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_K27me3H.bam ./preprocessed/sort/CCA043*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_K27me3T.bam ./preprocessed/sort/CCA044*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_K119ubH.bam ./preprocessed/sort/CCA045*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/H69_K119ubT.bam ./preprocessed/sort/CCA046*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Egi1_K27ac.bam ./preprocessed/sort/CCA047*.sorted.bam
qsub ~/samtools-0.1.19/samtools merge ./preprocessed/merge/Egi1_input.bam ./preprocessed/sort/CCA048*.sorted.bam

# batch 3
mv ./preprocessed/sort/CCA049*.sorted.bam ./preprocessed/merge/H69_DNMT3A.bam
mv ./preprocessed/sort/CCA050*.sorted.bam ./preprocessed/merge/H69_DNMT3B.bam
mv ./preprocessed/sort/CCA051*.sorted.bam ./preprocessed/merge/H69_RING1A.bam
mv ./preprocessed/sort/CCA052*.sorted.bam ./preprocessed/merge/H69_RING1B.bam
mv ./preprocessed/sort/CCA053*.sorted.bam ./preprocessed/merge/H69_inputTF.bam
mv ./preprocessed/sort/CCA054*.sorted.bam ./preprocessed/merge/H69_K27me3.bam
mv ./preprocessed/sort/CCA055*.sorted.bam ./preprocessed/merge/H69_K119ub.bam
mv ./preprocessed/sort/CCA056*.sorted.bam ./preprocessed/merge/H69_input.bam
mv ./preprocessed/sort/CCA057*.sorted.bam ./preprocessed/merge/H69_KO_EV_K27ac.bam
mv ./preprocessed/sort/CCA058*.sorted.bam ./preprocessed/merge/H69_KO_EV_input.bam
mv ./preprocessed/sort/CCA059*.sorted.bam ./preprocessed/merge/H69_KO_WT_K27ac.bam
mv ./preprocessed/sort/CCA060*.sorted.bam ./preprocessed/merge/H69_KO_WT_input.bam
mv ./preprocessed/sort/CCA061*.sorted.bam ./preprocessed/merge/H69_KO_MUT_K27ac.bam
mv ./preprocessed/sort/CCA062*.sorted.bam ./preprocessed/merge/H69_KO_MUT_input.bam
mv ./preprocessed/sort/CCA063*.sorted.bam ./preprocessed/merge/TBTS5_EV_K27ac.bam
mv ./preprocessed/sort/CCA064*.sorted.bam ./preprocessed/merge/TBTS5_EV_input.bam
mv ./preprocessed/sort/CCA065*.sorted.bam ./preprocessed/merge/TBTS5_WT_K27ac.bam
mv ./preprocessed/sort/CCA066*.sorted.bam ./preprocessed/merge/TBTS5_WT_input.bam
mv ./preprocessed/sort/CCA067*.sorted.bam ./preprocessed/merge/TBTS5_MUT_K27ac.bam
mv ./preprocessed/sort/CCA068*.sorted.bam ./preprocessed/merge/TBTS5_MUT_input.bam

qsub samtools merge ./preprocessed/merge/H69_ASXL2.bam ./preprocessed/sort/CCA069*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_EZH2.bam ./preprocessed/sort/CCA070*.sorted.bam

# batch 4
qsub samtools merge ./preprocessed/merge/A035T_K27ac.bam ./preprocessed/sort/CCA071*.sorted.bam
qsub samtools merge ./preprocessed/merge/A035T_input.bam ./preprocessed/sort/CCA072*.sorted.bam
qsub samtools merge ./preprocessed/merge/C080T_K27ac.bam ./preprocessed/sort/CCA073*.sorted.bam
qsub samtools merge ./preprocessed/merge/C080T_input.bam ./preprocessed/sort/CCA074*.sorted.bam
qsub samtools merge ./preprocessed/merge/C096T_K27ac.bam ./preprocessed/sort/CCA075*.sorted.bam
qsub samtools merge ./preprocessed/merge/C096T_input.bam ./preprocessed/sort/CCA076*.sorted.bam
qsub samtools merge ./preprocessed/merge/W39T_K27ac.bam ./preprocessed/sort/CCA077*.sorted.bam
qsub samtools merge ./preprocessed/merge/W39T_input.bam ./preprocessed/sort/CCA078*.sorted.bam
qsub samtools merge ./preprocessed/merge/Y65T_K27ac.bam ./preprocessed/sort/CCA079*.sorted.bam
qsub samtools merge ./preprocessed/merge/Y65T_input.bam ./preprocessed/sort/CCA080*.sorted.bam
qsub samtools merge ./preprocessed/merge/1202T_K27ac.bam ./preprocessed/sort/CCA081*.sorted.bam
qsub samtools merge ./preprocessed/merge/1202T_input.bam ./preprocessed/sort/CCA082*.sorted.bam
qsub samtools merge ./preprocessed/merge/Y140T_K27ac.bam ./preprocessed/sort/CCA083*.sorted.bam
qsub samtools merge ./preprocessed/merge/Y140T_input.bam ./preprocessed/sort/CCA084*.sorted.bam
qsub samtools merge ./preprocessed/merge/33587479T_K27ac.bam ./preprocessed/sort/CCA085*.sorted.bam
qsub samtools merge ./preprocessed/merge/33587479T_input.bam ./preprocessed/sort/CCA086*.sorted.bam
qsub samtools merge ./preprocessed/merge/2000861T_K27ac.bam ./preprocessed/sort/CCA087*.sorted.bam
qsub samtools merge ./preprocessed/merge/2000861T_input.bam ./preprocessed/sort/CCA088*.sorted.bam
qsub samtools merge ./preprocessed/merge/17231203T_K27ac.bam ./preprocessed/sort/CCA089*.sorted.bam
qsub samtools merge ./preprocessed/merge/17231203T_input.bam ./preprocessed/sort/CCA090*.sorted.bam
qsub samtools merge ./preprocessed/merge/72600277T_K27ac.bam ./preprocessed/sort/CCA091*.sorted.bam
qsub samtools merge ./preprocessed/merge/72600277T_input.bam ./preprocessed/sort/CCA092*.sorted.bam
qsub samtools merge ./preprocessed/merge/Z3722T_K27ac.bam ./preprocessed/sort/CCA093*.sorted.bam
qsub samtools merge ./preprocessed/merge/Z3722T_input.bam ./preprocessed/sort/CCA094*.sorted.bam
qsub samtools merge ./preprocessed/merge/260418N_K27ac.bam ./preprocessed/sort/CCA095*.sorted.bam
qsub samtools merge ./preprocessed/merge/260418N_input.bam ./preprocessed/sort/CCA096*.sorted.bam
qsub samtools merge ./preprocessed/merge/807N_K27ac.bam ./preprocessed/sort/CCA097*.sorted.bam
qsub samtools merge ./preprocessed/merge/807N_input.bam ./preprocessed/sort/CCA098*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_K27ac.bam ./preprocessed/sort/CCA099*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_input.bam ./preprocessed/sort/CCA100*.sorted.bam
mv ./preprocessed/sort/CCA101*.sorted.bam ./preprocessed/merge/KKU100_K27ac.bam
mv ./preprocessed/sort/CCA102*.sorted.bam ./preprocessed/merge/KKU100_input.bam
mv ./preprocessed/sort/CCA103*.sorted.bam ./preprocessed/merge/M213_K27ac.bam
mv ./preprocessed/sort/CCA104*.sorted.bam ./preprocessed/merge/M213_input.bam

# batch 5
qsub samtools merge ./preprocessed/merge/H69_R1_K27ac.bam ./preprocessed/sort/CCA105*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_R1_input.bam ./preprocessed/sort/CCA106*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_KO_R1_K27ac.bam ./preprocessed/sort/CCA107*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_KO_R1_input.bam ./preprocessed/sort/CCA108*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_R2_K27ac.bam ./preprocessed/sort/CCA109*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_R2_input.bam ./preprocessed/sort/CCA110*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_KO_R2_K27ac.bam ./preprocessed/sort/CCA111*.sorted.bam
qsub samtools merge ./preprocessed/merge/H69_KO_R2_input.bam ./preprocessed/sort/CCA112*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_WT_R1_input.bam ./preprocessed/sort/CCA113*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_R1_K27ac.bam ./preprocessed/sort/CCA114*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_R1_input.bam ./preprocessed/sort/CCA115*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_KO_R1_K27ac.bam ./preprocessed/sort/CCA116*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_KO_R1_input.bam ./preprocessed/sort/CCA117*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_R2_K27ac.bam ./preprocessed/sort/CCA118*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_R2_input.bam ./preprocessed/sort/CCA119*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_KO_R2_K27ac.bam ./preprocessed/sort/CCA120*.sorted.bam
qsub samtools merge ./preprocessed/merge/Egi1_KO_R2_input.bam ./preprocessed/sort/CCA121*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_MUT_R1_K27ac.bam ./preprocessed/sort/CCA122*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_MUT_R1_input.bam ./preprocessed/sort/CCA123*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_R1_K27ac.bam ./preprocessed/sort/CCA124*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_R1_input.bam ./preprocessed/sort/CCA125*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_KO_R1_K27ac.bam ./preprocessed/sort/CCA126*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_KO_R1_input.bam ./preprocessed/sort/CCA127*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_R2_K27ac.bam ./preprocessed/sort/CCA128*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_R2_input.bam ./preprocessed/sort/CCA129*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_KO_R2_K27ac.bam ./preprocessed/sort/CCA130*.sorted.bam
qsub samtools merge ./preprocessed/merge/HUCCT_KO_R2_input.bam ./preprocessed/sort/CCA131*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_EV_R1_K27ac.bam ./preprocessed/sort/CCA132*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_EV_R1_input.bam ./preprocessed/sort/CCA133*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_WT_R1_K27ac.bam ./preprocessed/sort/CCA134*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_EV_R1_K27ac.bam ./preprocessed/sort/CCA135*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_EV_R1_input.bam ./preprocessed/sort/CCA136*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_WT_R1_K27ac.bam ./preprocessed/sort/CCA137*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_WT_R1_input.bam ./preprocessed/sort/CCA138*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_MUT_R1_K27ac.bam ./preprocessed/sort/CCA139*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_MUT_R1_input.bam ./preprocessed/sort/CCA140*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_EV_R2_K27ac.bam ./preprocessed/sort/CCA141*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_EV_R2_input.bam ./preprocessed/sort/CCA142*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_WT_R2_K27ac.bam ./preprocessed/sort/CCA143*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_WT_R2_input.bam ./preprocessed/sort/CCA144*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_MUT_R2_K27ac.bam ./preprocessed/sort/CCA145*.sorted.bam
qsub samtools merge ./preprocessed/merge/TFK1_MUT_R2_input.bam ./preprocessed/sort/CCA146*.sorted.bam


# batch 6
mv ./preprocessed/sort/CCA147*.sorted.bam ./preprocessed/merge/260418N_input.bam
mv ./preprocessed/sort/CCA148*.sorted.bam ./preprocessed/merge/260418N_K27ac.bam
mv ./preprocessed/sort/CCA149*.sorted.bam ./preprocessed/merge/3011118N_K27ac.bam
mv ./preprocessed/sort/CCA150*.sorted.bam ./preprocessed/merge/3011118N_input.bam
mv ./preprocessed/sort/CCA151*.sorted.bam ./preprocessed/merge/4081118N_K27ac.bam
mv ./preprocessed/sort/CCA152*.sorted.bam ./preprocessed/merge/4081118N_input.bam

qsub samtools merge ./preprocessed/merge/KKU100_R1_K27ac.bam ./preprocessed/sort/CCA153*.sorted.bam
qsub samtools merge ./preprocessed/merge/KKU100_R1_input.bam ./preprocessed/sort/CCA154*.sorted.bam
qsub samtools merge ./preprocessed/merge/MEC_K27ac.bam ./preprocessed/sort/CCA155*.sorted.bam
qsub samtools merge ./preprocessed/merge/MEC_input.bam ./preprocessed/sort/CCA156*.sorted.bam
qsub samtools merge ./preprocessed/merge/M156_K27ac.bam ./preprocessed/sort/CCA157*.sorted.bam
qsub samtools merge ./preprocessed/merge/M156_input.bam ./preprocessed/sort/CCA158*.sorted.bam
qsub samtools merge ./preprocessed/merge/M214_K27ac.bam ./preprocessed/sort/CCA159*.sorted.bam
qsub samtools merge ./preprocessed/merge/M214_input.bam ./preprocessed/sort/CCA160*.sorted.bam
qsub samtools merge ./preprocessed/merge/SNU1196_K27ac.bam ./preprocessed/sort/CCA161*.sorted.bam
qsub samtools merge ./preprocessed/merge/SNU1196_input.bam ./preprocessed/sort/CCA162*.sorted.bam
qsub samtools merge ./preprocessed/merge/SNU245_K27ac.bam ./preprocessed/sort/CCA163*.sorted.bam
qsub samtools merge ./preprocessed/merge/SNU245_input.bam ./preprocessed/sort/CCA164*.sorted.bam
qsub samtools merge ./preprocessed/merge/Z5606_K27ac.bam ./preprocessed/sort/CCA165*.sorted.bam
qsub samtools merge ./preprocessed/merge/Z5606_input.bam ./preprocessed/sort/CCA166*.sorted.bam
qsub samtools merge ./preprocessed/merge/M055_K27ac.bam ./preprocessed/sort/CCA167*.sorted.bam
qsub samtools merge ./preprocessed/merge/M055_input.bam ./preprocessed/sort/CCA168*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_K27ac.bam ./preprocessed/sort/CCA169*.sorted.bam
qsub samtools merge ./preprocessed/merge/TBTS5_input.bam ./preprocessed/sort/CCA170*.sorted.bam

# batch 7
mv ./preprocessed/sort/CHP19001_S*.sorted.bam ./preprocessed/merge/C150T_K27ac.bam
mv ./preprocessed/sort/CHP19002_S*.sorted.bam ./preprocessed/merge/2000123T_K27ac.bam
mv ./preprocessed/sort/CHP19003_S*.sorted.bam ./preprocessed/merge/31215286T_K27ac.bam
mv ./preprocessed/sort/CHP19004_S*.sorted.bam ./preprocessed/merge/45782570T_K27ac.bam
mv ./preprocessed/sort/CHP19005_S*.sorted.bam ./preprocessed/merge/B083T_K27ac.bam
mv ./preprocessed/sort/CHP19006_S*.sorted.bam ./preprocessed/merge/Y74T_K27ac.bam
mv ./preprocessed/sort/CHP19007_S*.sorted.bam ./preprocessed/merge/Z2403T_K27ac.bam
mv ./preprocessed/sort/CHP19008_S*.sorted.bam ./preprocessed/merge/23474504T_K27ac.bam
mv ./preprocessed/sort/CHP19009_S*.sorted.bam ./preprocessed/merge/1202T_K27ac.bam
mv ./preprocessed/sort/CHP19010_S*.sorted.bam ./preprocessed/merge/3001T_K27ac.bam
mv ./preprocessed/sort/CHP19011_S*.sorted.bam ./preprocessed/merge/824T_K27ac.bam
mv ./preprocessed/sort/CHP19012_S*.sorted.bam ./preprocessed/merge/86014838T_K27ac.bam
mv ./preprocessed/sort/CHP19013_S*.sorted.bam ./preprocessed/merge/90866096T_K27ac.bam
mv ./preprocessed/sort/CHP19014_S*.sorted.bam ./preprocessed/merge/A035T_K27ac.bam
mv ./preprocessed/sort/CHP19015_S*.sorted.bam ./preprocessed/merge/Y140T_K27ac.bam
mv ./preprocessed/sort/CHP19016_S*.sorted.bam ./preprocessed/merge/Y65T_K27ac.bam
mv ./preprocessed/sort/CHP19017_S*.sorted.bam ./preprocessed/merge/Z3722T_K27ac.bam
mv ./preprocessed/sort/CHP19018_S*.sorted.bam ./preprocessed/merge/Z639T_K27ac.bam
mv ./preprocessed/sort/CHP19019_S*.sorted.bam ./preprocessed/merge/A003T_K27ac.bam
mv ./preprocessed/sort/CHP19020_S*.sorted.bam ./preprocessed/merge/A042T_K27ac.bam
mv ./preprocessed/sort/CHP19021_S*.sorted.bam ./preprocessed/merge/A074T_K27ac.bam
mv ./preprocessed/sort/CHP19022_S*.sorted.bam ./preprocessed/merge/A096T_K27ac.bam
mv ./preprocessed/sort/CHP19023_S*.sorted.bam ./preprocessed/merge/A142T_K27ac.bam
mv ./preprocessed/sort/CHP19025_S*.sorted.bam ./preprocessed/merge/B011T_K27ac.bam
mv ./preprocessed/sort/CHP19026_S*.sorted.bam ./preprocessed/merge/B085T_K27ac.bam
mv ./preprocessed/sort/CHP19027_S*.sorted.bam ./preprocessed/merge/C008T_K27ac.bam
mv ./preprocessed/sort/CHP19028_S*.sorted.bam ./preprocessed/merge/C078T_K27ac.bam
mv ./preprocessed/sort/CHP19029_S*.sorted.bam ./preprocessed/merge/C105T_K27ac.bam
mv ./preprocessed/sort/CHP19030_S*.sorted.bam ./preprocessed/merge/C144T_K27ac.bam
mv ./preprocessed/sort/CHP19031_S*.sorted.bam ./preprocessed/merge/A157T_K27ac.bam
mv ./preprocessed/sort/CHP19032_S*.sorted.bam ./preprocessed/merge/C167T_K27ac.bam
mv ./preprocessed/sort/CHP19033_S*.sorted.bam ./preprocessed/merge/C176T_K27ac.bam
mv ./preprocessed/sort/CHP19034_S*.sorted.bam ./preprocessed/merge/D017T_K27ac.bam
mv ./preprocessed/sort/CHP19035_S*.sorted.bam ./preprocessed/merge/Y002T_K27ac.bam
mv ./preprocessed/sort/CHP19036_S*.sorted.bam ./preprocessed/merge/Y091T_K27ac.bam
mv ./preprocessed/sort/CHP19037_S*.sorted.bam ./preprocessed/merge/19286759T_K27ac.bam
mv ./preprocessed/sort/CHP19038_S*.sorted.bam ./preprocessed/merge/93199056T_K27ac.bam
mv ./preprocessed/sort/CHP19039_S*.sorted.bam ./preprocessed/merge/97203223T_K27ac.bam
mv ./preprocessed/sort/CHP19040_S*.sorted.bam ./preprocessed/merge/61271150T_K27ac.bam
mv ./preprocessed/sort/CHP19041_S*.sorted.bam ./preprocessed/merge/Z2778T_K27ac.bam
mv ./preprocessed/sort/CHP19042_S*.sorted.bam ./preprocessed/merge/Z321T_K27ac.bam
mv ./preprocessed/sort/CHP19043_S*.sorted.bam ./preprocessed/merge/77071507T_K27ac.bam
mv ./preprocessed/sort/CHP19044_S*.sorted.bam ./preprocessed/merge/Z508T_K27ac.bam
mv ./preprocessed/sort/CHP19045_S*.sorted.bam ./preprocessed/merge/26814273T_K27ac.bam
mv ./preprocessed/sort/CHP19046_S*.sorted.bam ./preprocessed/merge/A153T_K27ac.bam
mv ./preprocessed/sort/CHP19047_S*.sorted.bam ./preprocessed/merge/W40T_K27ac.bam
mv ./preprocessed/sort/CHP19048_S*.sorted.bam ./preprocessed/merge/A059T_K27ac.bam
mv ./preprocessed/sort/CHP19049_S*.sorted.bam ./preprocessed/merge/R149T_K27ac.bam
mv ./preprocessed/sort/CHP19050_S*.sorted.bam ./preprocessed/merge/30131247T_K27ac.bam
mv ./preprocessed/sort/CHP19052_S*.sorted.bam ./preprocessed/merge/TAIWAN_33T_K27ac.bam
mv ./preprocessed/sort/CHP19053_S*.sorted.bam ./preprocessed/merge/TAIWAN_38T_K27ac.bam
mv ./preprocessed/sort/CHP19054_S*.sorted.bam ./preprocessed/merge/TAIWAN_39T_K27ac.bam
mv ./preprocessed/sort/CHP19055_S*.sorted.bam ./preprocessed/merge/TAIWAN_43T_K27ac.bam
mv ./preprocessed/sort/CHP19056_S*.sorted.bam ./preprocessed/merge/TAIWAN_46T_K27ac.bam
mv ./preprocessed/sort/CHP19057_S*.sorted.bam ./preprocessed/merge/TAIWAN_684T_K27ac.bam
mv ./preprocessed/sort/CHP19058_S*.sorted.bam ./preprocessed/merge/17231203T_K27ac.bam
mv ./preprocessed/sort/CHP19059_S*.sorted.bam ./preprocessed/merge/A169T_K27ac.bam
mv ./preprocessed/sort/CHP19060_S*.sorted.bam ./preprocessed/merge/Z12243N_K27ac.bam
mv ./preprocessed/sort/CHP19061_S*.sorted.bam ./preprocessed/merge/Z12244N_K27ac.bam
mv ./preprocessed/sort/CHP19062_S*.sorted.bam ./preprocessed/merge/Z12249N_K27ac.bam
mv ./preprocessed/sort/CHP19063_S*.sorted.bam ./preprocessed/merge/Z12267N_K27ac.bam
mv ./preprocessed/sort/CHP19064_S*.sorted.bam ./preprocessed/merge/2143T_K27ac.bam


mv ./preprocessed/sort/CHP19001_IN_S*.sorted.bam ./preprocessed/merge/C150T_input.bam
mv ./preprocessed/sort/CHP19002_IN_S*.sorted.bam ./preprocessed/merge/2000123T_input.bam
mv ./preprocessed/sort/CHP19003_IN_S*.sorted.bam ./preprocessed/merge/31215286T_input.bam
mv ./preprocessed/sort/CHP19004_IN_S*.sorted.bam ./preprocessed/merge/45782570T_input.bam
mv ./preprocessed/sort/CHP19005_IN_S*.sorted.bam ./preprocessed/merge/B083T_input.bam
mv ./preprocessed/sort/CHP19006_IN_S*.sorted.bam ./preprocessed/merge/Y74T_input.bam
mv ./preprocessed/sort/CHP19007_IN_S*.sorted.bam ./preprocessed/merge/Z2403T_input.bam
mv ./preprocessed/sort/CHP19008_IN_S*.sorted.bam ./preprocessed/merge/23474504T_input.bam
mv ./preprocessed/sort/CHP19009_IN_S*.sorted.bam ./preprocessed/merge/1202T_input.bam
mv ./preprocessed/sort/CHP19010_IN_S*.sorted.bam ./preprocessed/merge/3001T_input.bam
mv ./preprocessed/sort/CHP19011_IN_S*.sorted.bam ./preprocessed/merge/824T_input.bam
mv ./preprocessed/sort/CHP19012_IN_S*.sorted.bam ./preprocessed/merge/86014838T_input.bam
mv ./preprocessed/sort/CHP19013_IN_S*.sorted.bam ./preprocessed/merge/90866096T_input.bam
mv ./preprocessed/sort/CHP19014_IN_S*.sorted.bam ./preprocessed/merge/A035T_input.bam
mv ./preprocessed/sort/CHP19015_IN_S*.sorted.bam ./preprocessed/merge/Y140T_input.bam
mv ./preprocessed/sort/CHP19016_IN_S*.sorted.bam ./preprocessed/merge/Y65T_input.bam
mv ./preprocessed/sort/CHP19017_IN_S*.sorted.bam ./preprocessed/merge/Z3722T_input.bam
mv ./preprocessed/sort/CHP19018_IN_S*.sorted.bam ./preprocessed/merge/Z639T_input.bam
mv ./preprocessed/sort/CHP19019_IN_S*.sorted.bam ./preprocessed/merge/A003T_input.bam
mv ./preprocessed/sort/CHP19020_IN_S*.sorted.bam ./preprocessed/merge/A042T_input.bam
mv ./preprocessed/sort/CHP19021_IN_S*.sorted.bam ./preprocessed/merge/A074T_input.bam
mv ./preprocessed/sort/CHP19022_IN_S*.sorted.bam ./preprocessed/merge/A096T_input.bam
mv ./preprocessed/sort/CHP19023_IN_S*.sorted.bam ./preprocessed/merge/A142T_input.bam
mv ./preprocessed/sort/CHP19025_IN_S*.sorted.bam ./preprocessed/merge/B011T_input.bam
mv ./preprocessed/sort/CHP19026_IN_S*.sorted.bam ./preprocessed/merge/B085T_input.bam
mv ./preprocessed/sort/CHP19027_IN_S*.sorted.bam ./preprocessed/merge/C008T_input.bam
mv ./preprocessed/sort/CHP19028_IN_S*.sorted.bam ./preprocessed/merge/C078T_input.bam
mv ./preprocessed/sort/CHP19029_IN_S*.sorted.bam ./preprocessed/merge/C105T_input.bam
mv ./preprocessed/sort/CHP19030_IN_S*.sorted.bam ./preprocessed/merge/C144T_input.bam
mv ./preprocessed/sort/CHP19031_IN_S*.sorted.bam ./preprocessed/merge/A157T_input.bam
mv ./preprocessed/sort/CHP19032_IN_S*.sorted.bam ./preprocessed/merge/C167T_input.bam
mv ./preprocessed/sort/CHP19033_IN_S*.sorted.bam ./preprocessed/merge/C176T_input.bam
mv ./preprocessed/sort/CHP19034_IN_S*.sorted.bam ./preprocessed/merge/D017T_input.bam
mv ./preprocessed/sort/CHP19035_IN_S*.sorted.bam ./preprocessed/merge/Y002T_input.bam
mv ./preprocessed/sort/CHP19036_IN_S*.sorted.bam ./preprocessed/merge/Y091T_input.bam
mv ./preprocessed/sort/CHP19037_IN_S*.sorted.bam ./preprocessed/merge/19286759T_input.bam
mv ./preprocessed/sort/CHP19038_IN_S*.sorted.bam ./preprocessed/merge/93199056T_input.bam
mv ./preprocessed/sort/CHP19039_IN_S*.sorted.bam ./preprocessed/merge/97203223T_input.bam
mv ./preprocessed/sort/CHP19040_IN_S*.sorted.bam ./preprocessed/merge/61271150T_input.bam
mv ./preprocessed/sort/CHP19041_IN_S*.sorted.bam ./preprocessed/merge/Z2778T_input.bam
mv ./preprocessed/sort/CHP19042_IN_S*.sorted.bam ./preprocessed/merge/Z321T_input.bam
mv ./preprocessed/sort/CHP19043_IN_S*.sorted.bam ./preprocessed/merge/77071507T_input.bam
mv ./preprocessed/sort/CHP19044_IN_S*.sorted.bam ./preprocessed/merge/Z508T_input.bam
mv ./preprocessed/sort/CHP19045_IN_S*.sorted.bam ./preprocessed/merge/26814273T_input.bam
mv ./preprocessed/sort/CHP19046_IN_S*.sorted.bam ./preprocessed/merge/A153T_input.bam
mv ./preprocessed/sort/CHP19047_IN_S*.sorted.bam ./preprocessed/merge/W40T_input.bam
mv ./preprocessed/sort/CHP19048_IN_S*.sorted.bam ./preprocessed/merge/A059T_input.bam
mv ./preprocessed/sort/CHP19049_IN_S*.sorted.bam ./preprocessed/merge/R149T_input.bam
mv ./preprocessed/sort/CHP19050_IN_S*.sorted.bam ./preprocessed/merge/30131247T_input.bam
mv ./preprocessed/sort/CHP19052_IN_S*.sorted.bam ./preprocessed/merge/TAIWAN_33T_input.bam
mv ./preprocessed/sort/CHP19053_IN_S*.sorted.bam ./preprocessed/merge/TAIWAN_38T_input.bam
mv ./preprocessed/sort/CHP19054_IN_S*.sorted.bam ./preprocessed/merge/TAIWAN_39T_input.bam
mv ./preprocessed/sort/CHP19055_IN_S*.sorted.bam ./preprocessed/merge/TAIWAN_43T_input.bam
mv ./preprocessed/sort/CHP19056_IN_S*.sorted.bam ./preprocessed/merge/TAIWAN_46T_input.bam
mv ./preprocessed/sort/CHP19057_IN_S*.sorted.bam ./preprocessed/merge/TAIWAN_684T_input.bam
mv ./preprocessed/sort/CHP19058_IN_S*.sorted.bam ./preprocessed/merge/17231203T_input.bam
mv ./preprocessed/sort/CHP19059_IN_S*.sorted.bam ./preprocessed/merge/A169T_input.bam
mv ./preprocessed/sort/CHP19060_IN_S*.sorted.bam ./preprocessed/merge/Z12243N_input.bam
mv ./preprocessed/sort/CHP19061_IN_S*.sorted.bam ./preprocessed/merge/Z12244N_input.bam
mv ./preprocessed/sort/CHP19062_IN_S*.sorted.bam ./preprocessed/merge/Z12249N_input.bam
mv ./preprocessed/sort/CHP19063_IN_S*.sorted.bam ./preprocessed/merge/Z12267N_input.bam
mv ./preprocessed/sort/CHP19064_IN_S*.sorted.bam ./preprocessed/merge/2143T_input.bam


mv ./preprocessed/sort/A100T_H3_*.sorted.bam ./preprocessed/merge/A100T_K27ac.bam
mv ./preprocessed/sort/A100T_I_*.sorted.bam ./preprocessed/merge/A100T_input.bam
mv ./preprocessed/sort/H3_14T_*.sorted.bam ./preprocessed/merge/TAIWAN_14T_K27ac.bam
mv ./preprocessed/sort/I_14T_*.sorted.bam ./preprocessed/merge/TAIWAN_14T_input.bam
mv ./preprocessed/sort/H3_21980572T_*.sorted.bam ./preprocessed/merge/29150572T_K27ac.bam
mv ./preprocessed/sort/I_21980572T_*.sorted.bam ./preprocessed/merge/29150572T_input.bam
mv ./preprocessed/sort/K_21914178T_*.sorted.bam ./preprocessed/merge/21914187T_K27ac.bam
mv ./preprocessed/sort/I_21914178T_*.sorted.bam ./preprocessed/merge/21914187T_input.bam
mv ./preprocessed/sort/Z3585T_H3_*.sorted.bam ./preprocessed/merge/Z3585T_K27ac.bam
mv ./preprocessed/sort/Z3585T_I_*.sorted.bam ./preprocessed/merge/Z3585T_input.bam
mv ./preprocessed/sort/H3_2KU452_*.sorted.bam ./preprocessed/merge/2KU452_K27ac.bam
mv ./preprocessed/sort/I_2KU452_*.sorted.bam ./preprocessed/merge/2KU452_input.bam
mv ./preprocessed/sort/K_TFK1_*.sorted.bam ./preprocessed/merge/TFK1_K27ac.bam
mv ./preprocessed/sort/I_TFK1_*.sorted.bam ./preprocessed/merge/TFK1_input.bam

# batch 9
mv ./preprocessed/sort/CCA171*.sorted.bam ./preprocessed/merge/HUCCT_LSD1.bam
mv ./preprocessed/sort/CCA172*.sorted.bam ./preprocessed/merge/HUCCT_KO_LSD1.bam
mv ./preprocessed/sort/CCA173*.sorted.bam ./preprocessed/merge/HUCCT_SP_LSD1.bam
mv ./preprocessed/sort/CCA174*.sorted.bam ./preprocessed/merge/HUCCT_KO_SP_LSD1.bam
mv ./preprocessed/sort/CCA175*.sorted.bam ./preprocessed/merge/HUCCT_BAP1.bam
mv ./preprocessed/sort/CCA176*.sorted.bam ./preprocessed/merge/HUCCT_KO_BAP1.bam
mv ./preprocessed/sort/CCA177*.sorted.bam ./preprocessed/merge/HUCCT_SP_BAP1.bam
mv ./preprocessed/sort/CCA178*.sorted.bam ./preprocessed/merge/HUCCT_KO_SP_BAP1.bam
mv ./preprocessed/sort/CCA179*.sorted.bam ./preprocessed/merge/HUCCT_PARP1.bam
mv ./preprocessed/sort/CCA180*.sorted.bam ./preprocessed/merge/HUCCT_KO_PARP1.bam
mv ./preprocessed/sort/CCA181*.sorted.bam ./preprocessed/merge/HUCCT_SP_PARP1.bam
mv ./preprocessed/sort/CCA182*.sorted.bam ./preprocessed/merge/HUCCT_KO_SP_PARP1.bam
mv ./preprocessed/sort/CCA183*.sorted.bam ./preprocessed/merge/HUCCT_H2AK119UB.bam
mv ./preprocessed/sort/CCA184*.sorted.bam ./preprocessed/merge/HUCCT_KO_H2AK119UB.bam
mv ./preprocessed/sort/CCA185*.sorted.bam ./preprocessed/merge/HUCCT_SP_H2AK119UB.bam
mv ./preprocessed/sort/CCA186*.sorted.bam ./preprocessed/merge/HUCCT_KO_SP_H2AK119UB.bam
mv ./preprocessed/sort/CCA187*.sorted.bam ./preprocessed/merge/HUCCT_H3K9ME2.bam
mv ./preprocessed/sort/CCA188*.sorted.bam ./preprocessed/merge/HUCCT_KO_H3K9ME2.bam
mv ./preprocessed/sort/CCA189*.sorted.bam ./preprocessed/merge/HUCCT_SP_H3K9ME2.bam
mv ./preprocessed/sort/CCA190*.sorted.bam ./preprocessed/merge/HUCCT_KO_SP_H3K9ME2.bam
mv ./preprocessed/sort/CCA191*.sorted.bam ./preprocessed/merge/HUCCT_input.bam
mv ./preprocessed/sort/CCA192*.sorted.bam ./preprocessed/merge/HUCCT_KO_input.bam
mv ./preprocessed/sort/CCA193*.sorted.bam ./preprocessed/merge/HUCCT_SP_input.bam
mv ./preprocessed/sort/CCA194*.sorted.bam ./preprocessed/merge/HUCCT_KO_SP_input.bam




# count number of reads
for i in ./preprocessed/merge/*.bam
do
  qsub -o $i.numreads samtools view -c $i
done
#rm preprocessed/numreads_2_align.txt
for i in preprocessed/merge/*.numreads
do
  echo $i >> preprocessed/numreads_2_align.txt
  cat $i >> preprocessed/numreads_2_align.txt
done




# filter by mapping quality
mkdir preprocessed/mapq
for i in ./preprocessed/merge/*.bam
do
  basefilename=$(basename $i .bam)
  qsub -o ./preprocessed/mapq/$basefilename.mapq.bam samtools view -b -q 10 $i
done

# count number of reads
for i in ./preprocessed/mapq/*.bam
do
  qsub -o $i.numreads samtools view -c $i
done
#rm preprocessed/numreads_3_mapq.txt
for i in preprocessed/mapq/*.numreads
do
  echo $i >> preprocessed/numreads_3_mapq.txt
  cat $i >> preprocessed/numreads_3_mapq.txt
done




# remove dups
mkdir preprocessed/rmdup
for i in ./preprocessed/mapq/*.bam
do
  basefilename=$(basename $i .bam)
  qsub ~/samtools-0.1.19/samtools rmdup -s $i ./preprocessed/rmdup/$basefilename.rmdup.bam
done



# count number of reads
for i in ./preprocessed/rmdup/*.bam
do
  qsub -o $i.numreads samtools view -c $i
done
#rm preprocessed/numreads_4_rmdup.txt
for i in preprocessed/rmdup/*.numreads
do
  echo $i >> preprocessed/numreads_4_rmdup.txt
  cat $i >> preprocessed/numreads_4_rmdup.txt
done





# sort, generate index
mkdir final_bam_batch9
cd final_bam_batch9
for i in ../preprocessed/rmdup/*.bam
do
  basefilename=$(basename $i .mapq.rmdup.bam)
  qsub samtools sort -o $basefilename.bam $i
done
for i in *.bam
do
  qsub samtools index $i
done
for i in *.bam
do
  qsub -o $i.numreads samtools view -c $i
done






# check chip-seq quality (NSC, RSC) with phantompeakqualtools https://github.com/kundajelab/phantompeakqualtools
cd ~/CCA/CCA_ChipSeq/processed
mkdir preprocessed/phantompeakqual
cd ~/CCA/CCA_ChipSeq/processed/final_bam_batch5
for i in *.bam
do
  basefilename=$(basename $i)
  qsub /gpfs/apps/R/3.1.0/bin/Rscript ~/phantompeakqualtools/phantompeakqualtools-master/run_spp.R -c=$i -savp -odir=../preprocessed/phantompeakqual -out=../preprocessed/phantompeakqual/$basefilename.phantom.out
done



################### PEAK CALLING ################
Encode pipeline says H3K27ac is narrow peak, but H3K4me1 is broad peak:
https://www.encodeproject.org/chip-seq/histone/
broad peak marks: H3F3A, H3K27me3, H3K36me3, H3K4me1, H3K79me2, H3K79me3, H3K9me1, H3K9me2, H4K20me1
narrow peak marks: H2AFZ, H3ac, H3K27ac, H3K4me2, H3K4me3, H3K9ac
Dont estimate extsize in MACS2, just set extsize to 150. This makes downstream analysis easier so we dont need to set different extsize for each sample!

# K27ac: narrowpeak, ext150, P=.001
# Justfication for using liberal cutoff (P instead of Q): this is just to identify regions of interest for downstream analyses.
# Downstream analyses will quantify these regions and analyze based on that (eg. diffbind).
# Also tried: P=.01 (encode, but they further filter with IDR): many peaks!
mkdir ~/CCA/CCA_ChipSeq/processed/macs_peaks_batch8
cd ~/CCA/CCA_ChipSeq/processed/macs_peaks_batch8
for i in ../final_bam_batch8/*_K27ac.bam
do
  basefilename=$(basename $i _K27ac.bam)
  dirname=$(dirname $i)
  qsub macs2 callpeak -t $i -c ${dirname}/${basefilename}_input.bam -f BAM -n ${basefilename}_K27ac --nomodel --extsize 150 -g hs -p 0.001
done


# encode
for i in ../bams/*K27ac.bam
do
  basefilename=$(basename $i -H3K27ac.bam)
  dirname=$(dirname $i)
  qsub macs2 callpeak -t $i -c ${dirname}/${basefilename}-Input.bam -f BAM -n ${basefilename}_K27ac --nomodel --extsize 150 -g hs -p 0.001
done





# filter the peaks: keep only chr 1-22, sort
for i in *.narrowPeak
do
   basefilename=$(basename $i .narrowPeak)
   grep -E "^chr[0-9]+\s" $i > $basefilename.filt2.narrowPeak
   sort -k1,1V -k2,2n $basefilename.filt2.narrowPeak > $basefilename.filt.narrowPeak
done
rm *.filt2.narrowPeak







#################### Filter peaks against Encode blacklist ###############################
# Download blacklist from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/
# Remove peaks contained completely in blacklist.
# Previously: subtract blacklist region from each peak region, if remaining peak region is > 100bp, keep it.
for i in *.filt.narrowPeak
do
  basefilename=$(basename $i .filt.narrowPeak)
  bedtools intersect -a $i -b ~/CCA/genome/encode_blacklist/wgEncodeHg19ConsensusSignalArtifactRegions.bed -v -f 1.0 > $basefilename.filt3.narrowPeak
  #bedtools subtract -a $i -b ~/genome/encode_blacklist/wgEncodeHg19ConsensusSignalArtifactRegions.bed > $basefilename.filt2.narrowPeak
  #cat $basefilename.filt2.narrowPeak | awk '{if ($3-$2>=100) {print $0}}' > $basefilename.filt3.narrowPeak
  #rm $basefilename.filt2.narrowPeak
done
for i in *.filt3.narrowPeak
do
  basefilename=$(basename $i .filt3.narrowPeak)
  sort -k1,1V -k2,2n $i > $i.sort
  mv $i.sort $basefilename.filt.narrowPeak
  rm $i
done






################# Filter peaks to remove those near promoters. Enhancers = H3K27ac peaks - promoter region (TSS+/-2.5kb) #################
# Generate TSS+/-2.5kb annotations from UCSC hg19 refseq (from ROSE annotations)
# get latest refseq curated list from UCSC table browser: genes_refseq_curated.txt
# remove non-coding RNAs
cd ~/genome
grep "\\sNM_" genes_refseq_curated.txt > genes_refseq_curated_onlyNM.txt
cat genes_refseq_curated_onlyNM.txt | awk '{if ($4=="+") {print $3 "\t" $5 "\t" $13} else if ($4=="-") {print $3 "\t" $6 "\t" $13} }' > genes_refseq_curated_onlyNM.tss.txt
grep -E "chr[0-9]+\s|chr[X|Y]\s" genes_refseq_curated_onlyNM.tss.txt > genes_refseq_curated_onlyNM.tss.filt.txt
sort -k1,1V -k2,2n -k3,3d -u genes_refseq_curated_onlyNM.tss.filt.txt > genes_refseq_curated_onlyNM.tss.filt.sort.txt
cat genes_refseq_curated_onlyNM.tss.filt.sort.txt | awk '{ print $1 "\t" $2-2500 "\t" $2+2501 "\t" $3 }' > genes_refseq_curated_onlyNM.tss2.5kb.txt
bedtools merge -i genes_refseq_curated_onlyNM.tss2.5kb.txt > genes_refseq_curated_onlyNM.tss2.5kb.txt.merged
sort -k1,1V -k2,2n -k3,3n genes_refseq_curated_onlyNM.tss2.5kb.txt.merged > genes_refseq_curated_onlyNM.tss2.5kb.txt.merged.sort
mv genes_refseq_curated_onlyNM.tss2.5kb.txt.merged.sort genes_refseq_curated_onlyNM.tss2.5kb.merged.bed

rm genes_refseq_curated_onlyNM.txt
rm genes_refseq_curated_onlyNM.tss.txt
rm genes_refseq_curated_onlyNM.tss.filt.txt
rm genes_refseq_curated_onlyNM.tss.filt.sort.filt
rm genes_refseq_curated_onlyNM.tss2.5kb.txt
rm genes_refseq_curated_onlyNM.tss2.5kb.txt.merged


# Remove regions in promoter: subtract promoter region from each peak region, if remaining peak region is > 150bp, keep it.
# Previously, remove peaks contained completely in promoter (TSS +- 2.5kb). But after DiffBind recenters, it may keep peaks inside promoters!
for i in *.filt.narrowPeak
do
  basefilename=$(basename $i .filt.narrowPeak)
  #bedtools intersect -a $i -b ~/genome/genes_refseq_curated_onlyNM.tss2.5kb.merged.bed -v -f 1.0 > $basefilename.filt.outsideProm.narrowPeak
  #bedtools intersect -a $i -b ~/genome/genes_refseq_curated_onlyNM.tss2.5kb.merged.bed -wa -f 1.0 > $basefilename.filt.insideProm.narrowPeak
  bedtools subtract -a $i -b ~/CCA/genome/genes_refseq_curated_onlyNM.tss2.5kb.merged.bed > $basefilename.filt.outsideProm2.narrowPeak
  cat $basefilename.filt.outsideProm2.narrowPeak | awk '{if ($3-$2>=150) {print $0}}' > $basefilename.filt.outsideProm.narrowPeak
  rm $basefilename.filt.outsideProm2.narrowPeak
done
for i in *.filt.narrowPeak
do
  basefilename=$(basename $i .filt.narrowPeak)
  bedtools intersect -a $i -b ~/CCA/genome/genes_refseq_curated_onlyNM.tss2.5kb.merged.bed > $basefilename.filt.insideProm2.narrowPeak
  cat $basefilename.filt.insideProm2.narrowPeak | awk '{if ($3-$2>=150) {print $0}}' > $basefilename.filt.insideProm.narrowPeak
  rm $basefilename.filt.insideProm2.narrowPeak
done
for i in *.filt.outsideProm.narrowPeak
do
  sort -k1,1V -k2,2n $i > $i.sort
  mv $i.sort $i
done
for i in *.filt.insideProm.narrowPeak
do
  sort -k1,1V -k2,2n $i > $i.sort
  mv $i.sort $i
done











##########################################################################################################################################
### Generate signal tracks for visualization in UCSC gwenome browser https://github.com/taoliu/MACS/wiki/Build-Signal-Track

# call macs2 pileup
cd signal_tracks_batch8
for i in ../final_bam_batch8/*.bam
do
  qsub macs2 pileup -i $i -f BAM -o $i.bdg --extsize 150
done
mv ../final_bam_batch8/*.bdg .



# scale input control by chip libsize / input libsize
for i in *_K27ac.bam.bdg
do
  basefilename=$(basename $i _K27ac.bam.bdg)
  qsub /gpfs/apps/R/3.1.0/bin/Rscript --vanilla ../scale_input_signals_chiplibsize.R ${basefilename}_input.bam.bdg ../final_bam_batch8/${basefilename}_K27ac.bam.numreads ../final_bam_batch8/${basefilename}_input.bam.numreads ${basefilename}_input.bam.scaled.bdg
done

# generate input-subtracted pileups
for i in *_K27ac.bam.bdg
do
  basefilename=$(basename $i _K27ac.bam.bdg)
  qsub macs2 bdgcmp -t $i -c ${basefilename}_input.bam.scaled.bdg -o ${basefilename}_K27ac_subtract.bdg -m subtract
done
rm *_input.bam.bdg
rm *_input.bam.scaled.bdg
rm *_K27ac.bam.bdg



## Normalization method 2: scale by normalization factor (eg DESEQ2 or EDGER factors). Lib size: get average CHIP lib size
# NOTE: some sample names in Diffbind deseq2 normfacs output don't match the bdg files! eg. 260418N_merged=260418N. Modify the normfacs file.
for i in *_K27ac_subtract.bdg
do
  basefilename=$(basename $i _subtract.bdg)
  samplename=$(basename $i _K27ac_subtract.bdg)
  #qsub /gpfs/apps/R/3.1.0/bin/Rscript --vanilla ../scale_signals_normfac.R $i ../20190128_cca_K27ac_Deseq2NormFacs.txt $samplename 20565565 ${basefilename}_subtract_Deseq2NormFac.bdg
  qsub /gpfs/apps/R/3.1.0/bin/Rscript --vanilla ../scale_signals_normfac.R $i ../20190130_cca_cells_tisspeaks_K27ac_Deseq2NormFacs.txt $samplename 16409875 ${basefilename}_subtract_Deseq2NormFac.bdg
done
# filter: keep only chr1-22
for i in *_subtract_Deseq2NormFac.bdg
do
do
  basefilename=$(basename $i .bdg)
  qsub -o $basefilename.filt2.bdg grep -E "\"^chr[0-9]+\\s\"" $i
done
# filter: remove lines with subtracted signal <= 0
for i in *_subtract_Deseq2NormFac.filt2.bdg
do
  basefilename=$(basename $i .filt2.bdg)
  qsub -o $basefilename.filt.bdg awk "'{if (\$4 > 0) {print \$0}}'" $i
done
rm *_subtract_Deseq2NormFac.filt2.bdg
# slop
for i in *_subtract_Deseq2NormFac.filt.bdg
do
  qsub -o ${i}.slop bedtools slop -i $i -g ~/genome/human.hg19.chrlens -b 0
done
# clip
for i in *_subtract_Deseq2NormFac.filt.bdg.slop
do
  qsub ~/UCSC_tools/bedClip $i ~/genome/human.hg19.chrlens $i.clip
done
# sort
for i in *_subtract_Deseq2NormFac.filt.bdg.slop.clip
do
  qsub -o $i.sort LC_COLLATE=C sort -k1,1 -k2,2n $i
done
# bedGraphToBigWig
for i in *_subtract_Deseq2NormFac.filt.bdg.slop.clip.sort
do
  basefilename=$(basename $i .bdg.slop.clip.sort)
  qsub ~/UCSC_tools/bedGraphToBigWig $i ~/genome/human.hg19.chrlens $basefilename.bw
done
rm *_subtract_Deseq2NormFac.bdg
rm *_subtract_Deseq2NormFac.filt.bdg
rm *_subtract_Deseq2NormFac.filt.bdg.slop
rm *_subtract_Deseq2NormFac.filt.bdg.slop.clip









# host on cyverse http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Hosting
# 1. go to https://de.cyverse.org/de/
# 2. upload file
# 3. click on file, send to genome browser (this enables sharing)
# For track hub, see https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html
# - Edit files in \ucsc genome browser\myHub
# For custom track (not recommended):
# - Open UCSC genome browser, custom track, paste link. Note that for bedgraph file, must specify track type=wiggle_0









#######################################################################
############## Find motifs ############################################
#######################################################################


# generate from R analysis, copy the bedfiles to /home/gmsv1178/CCA/CCA_ChipSeq/processed/homer

#### Run Homer motif finding
# switch to conda python 3.7 environment, for Homer / deseq2
source activate py37
cd /home/gmsv1178/CCA/CCA_ChipSeq/processed/homer
for i in *.bed
do
  samplename=$(basename $i)
  echo $samplename
  qsub -pe smp 4 findMotifsGenome.pl $i ../../../genome/GRCh37.p13.genome.chr_only.fa ${samplename}.homer -size given -p 4
done
