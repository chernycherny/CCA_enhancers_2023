

library(DiffBind)

library(GenomicRanges)






###########################################################################################################################
############################# Enhancer --> target gene mapping using RNA-seq ##############################################
###########################################################################################################################



######## Calculate K27ac signals in promoters. This is to validate mapping by RNA-seq ######
# Do on cluster: get list of genes and their promoter regions (TSS +/- 1kb): 25850 ptomoters, 19225 unique genes. Transfer bed file to local
cd ~/genome
grep "\\sNM_" genes_refseq_curated.txt > genes_refseq_curated_onlyNM.txt
cat genes_refseq_curated_onlyNM.txt | awk '{if ($4=="+") {print $3 "\t" $5 "\t" $13} else if ($4=="-") {print $3 "\t" $6 "\t" $13} }' > genes_refseq_curated_onlyNM.tss.txt
grep -E "chr[0-9]+\s|chr[X|Y]\s" genes_refseq_curated_onlyNM.tss.txt > genes_refseq_curated_onlyNM.tss.filt.txt
sort -k1,1V -k2,2n -k3,3d -u genes_refseq_curated_onlyNM.tss.filt.txt > genes_refseq_curated_onlyNM.tss.filt.sort.txt
cat genes_refseq_curated_onlyNM.tss.filt.sort.txt | awk '{ print $1 "\t" $2-1000 "\t" $2+1001 "\t" $3 }' > genes_refseq_curated_onlyNM.tss1kb.bed
rm genes_refseq_curated_onlyNM.txt
rm genes_refseq_curated_onlyNM.tss.txt
rm genes_refseq_curated_onlyNM.tss.filt.txt
rm genes_refseq_curated_onlyNM.tss.filt.sort.txt

# read and format the bedfile
promoter_regions = read.table("target_gene_mapping/genes_refseq_curated_onlyNM.tss1kb.bed",stringsAsFactors = F, header=F, sep="\t")# 25850
colnames(promoter_regions) = c("chr","start","end","gene")
promoter_regions$start = promoter_regions$start+1
promoter_regions$promoter_id = paste0(promoter_regions$chr, ":", promoter_regions$start, "-", promoter_regions$end, "_", promoter_regions$gene)
promoter_regions$coord = paste0(promoter_regions$chr, ":", promoter_regions$start, "-", promoter_regions$end)
promoter_regions_gr = as(promoter_regions$coord, "GRanges")

### Use diffbind to calculate K27ac signals at promoters
# diffbind will merge overlapping regions!!! separate out the promoters into non-overlapping sets to calculate the signals
GetNonOverlappingRegions = function(regions){
  regions_overlap = data.frame(findOverlaps(regions, regions, ignore.strand=T))
  regions_overlap = regions_overlap[regions_overlap$queryHits != regions_overlap$subjectHits,]
  if (nrow(regions_overlap)==0) {
    return (list(regions))
  }
  else {
    reg1_idx = unique(regions_overlap[regions_overlap$queryHits < regions_overlap$subjectHits, "subjectHits"]) # 5212
    reg2_idx = 1:length(regions) # 25850
    reg2_idx = reg2_idx[!reg2_idx %in% reg1_idx] # 20638
    reg1 = regions[reg1_idx]
    reg2 = regions[reg2_idx]
    reg1res = GetNonOverlappingRegions(reg1)
    reg2res = GetNonOverlappingRegions(reg2)
    return (c(reg1res, reg2res))
  }
}
promoter_regions_gr_split = GetNonOverlappingRegions(promoter_regions_gr) # 11 non-overlapping sets
length(setdiff(promoter_regions_gr, do.call(c, promoter_regions_gr_split)))==0
length(promoter_regions_gr) == length(do.call(c, promoter_regions_gr_split))

# diffbind has bug where if there's only one chromosome in a set, it doesn't return the chromosome string properly!
# so add a dummy region to sets with only one chromosome
promoter_regions_gr_split_dummy = lapply(promoter_regions_gr_split, function(x){
  if (length(unique(runValue(seqnames(x)))) > 1) { x  }
  else {
    if (as.character(unique(runValue(seqnames(x))))=="chrX") {
      c(as("chrY:1-2", "GRanges"), x)
    }
    else { c(as("chrX:1-2", "GRanges"), x) }
  }
})

# calculate signals in each of the non-overlapping sets
load("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/Rdata_cca_k27ac_0_prelim_final.Rdata")
K27ac_dba_promoter = K27ac_dba
rm(K27ac_dba, K27ac_enhsup_dba, K27ac_origPeaks, superenhancers)
K27ac_dba_promoter_split = lapply(promoter_regions_gr_split_dummy, function(x){dba.count(K27ac_dba_promoter, peaks=x, fragmentSize = 150)})
K27ac_promoter_split = lapply(K27ac_dba_promoter_split, function(x){
  dba.peakset(x, bRetrieve = T, DataType = DBA_DATA_FRAME)
})
all.equal(sapply(K27ac_promoter_split, nrow), sapply(promoter_regions_gr_split_dummy, length)) # 25852 rows total
apply(sapply(K27ac_promoter_split, colnames), 1, function(x){length(unique(x))==1}) # should be all true

# combine the signals from the non-overlapping sets, remove the dummy regions
K27ac_promoter = do.call(rbind, K27ac_promoter_split) # 25852
K27ac_promoter = K27ac_promoter[!((K27ac_promoter$CHR=="chrX" | K27ac_promoter$CHR=="chrY") & K27ac_promoter$START==1 & K27ac_promoter$END==2), ] # 25850
nrow(K27ac_promoter) == nrow(promoter_regions)

# casting to dataframe may change some of the sample names, so rename them. First check that sample ordering is as expected
as.character(dba.show(K27ac_dba_promoter)$ID)
colnames(K27ac_promoter)
# rename samples, columns
colnames(K27ac_promoter)[4:ncol(K27ac_promoter)] = as.character(dba.show(K27ac_dba_promoter)$ID)
colnames(K27ac_promoter)[colnames(K27ac_promoter)=="CHR"] = "chr"
colnames(K27ac_promoter)[colnames(K27ac_promoter)=="START"] = "start"
colnames(K27ac_promoter)[colnames(K27ac_promoter)=="END"] = "end"
rownames(K27ac_promoter) = NULL

# add gene column
K27ac_promoter$coord = paste0(K27ac_promoter$chr, ":", K27ac_promoter$start, "-", K27ac_promoter$end)
K27ac_promoter = K27ac_promoter[!duplicated(K27ac_promoter$coord),] # 25764
K27ac_promoter = merge(promoter_regions[, c("gene", "coord", "promoter_id")], K27ac_promoter, by="coord", sort=F) # 25850

# remove sex chr
K27ac_promoter = K27ac_promoter[!K27ac_promoter$chr %in% c("chrX", "chrY"),] # 24688

save(K27ac_promoter, K27ac_dba_promoter_split, file="C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/Rdata_cca_k27ac_1_targetmapping_promotersignal.Rdata")








############# Calculate correlation between each enhancer and expression of genes in its TAD  ###########


rm(list=ls())
library(DiffBind)

# load enhancer signals, delete unneeded stuff to save memory
load("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/Rdata_cca_k27ac_0_prelim_final.Rdata")
rm(K27ac_enhsup_dba, K27ac_origPeaks, superenhancers)

### enhancer signals using TMM-minus (full lib size). This accounts for lib size, but not enhancer size
K27ac_dba = dba.count(K27ac_dba, peaks=NULL, score=DBA_SCORE_TMM_MINUS_FULL)
K27ac_enhancer_signalsTMM = dba.peakset(K27ac_dba, bRetrieve = T, DataType = DBA_DATA_FRAME) # 474060 enhancers
rownames(K27ac_enhancer_signalsTMM) = paste0(K27ac_enhancer_signalsTMM$CHR, ":", K27ac_enhancer_signalsTMM$START, "-", K27ac_enhancer_signalsTMM$END)
K27ac_enhancer_signalsTMM = K27ac_enhancer_signalsTMM[, 4:ncol(K27ac_enhancer_signalsTMM)]
# casting to dataframe may change some of the sample names, so rename them. First check that sample ordering is as expected
as.character(dba.show(K27ac_dba)$ID)
colnames(K27ac_enhancer_signalsTMM)
# rename samples
colnames(K27ac_enhancer_signalsTMM) = as.character(dba.show(K27ac_dba)$ID)
# rank transform
K27ac_enhancer_signalsTMM = data.frame(apply(K27ac_enhancer_signalsTMM, 2, function(x){rank(x, ties.method = "average")}), check.names = F)
K27ac_enhancer_signalsTMM$enhancer = rownames(K27ac_enhancer_signalsTMM)



### enhancer signals using RPKM. This accounts for lib size and enhancer size, but doesn't minus control reads
K27ac_dba = dba.count(K27ac_dba, peaks=NULL, score=DBA_SCORE_RPKM)
K27ac_enhancer_signalsRPKM = dba.peakset(K27ac_dba, bRetrieve = T, DataType = DBA_DATA_FRAME) # 474060 enhancers
rownames(K27ac_enhancer_signalsRPKM) = paste0(K27ac_enhancer_signalsRPKM$CHR, ":", K27ac_enhancer_signalsRPKM$START, "-", K27ac_enhancer_signalsRPKM$END)
K27ac_enhancer_signalsRPKM = K27ac_enhancer_signalsRPKM[, 4:ncol(K27ac_enhancer_signalsRPKM)]
# rename samples
colnames(K27ac_enhancer_signalsRPKM) = as.character(dba.show(K27ac_dba)$ID)
# rank transform
K27ac_enhancer_signalsRPKM = data.frame(apply(K27ac_enhancer_signalsRPKM, 2, function(x){rank(x, ties.method = "average")}), check.names = F)
K27ac_enhancer_signalsRPKM$enhancer = rownames(K27ac_enhancer_signalsRPKM)


### enhancer signals using counts-minus per log width (ie kind of like density). This accounts for enhancer size, but not lib size, but libsize doesn't matter because we rank-normalize later.
# We don't do counts-minus per kb, because longer enhancers will have lower density! Thus take log width instead.
# After rank normalization, this is almost equivalent to using TMM / logwidth, EXCEPT that TMM has pseudocount of 1, but CPLW does not.
K27ac_dba = dba.count(K27ac_dba, peaks=NULL, score=DBA_SCORE_READS_MINUS)
K27ac_enhancer_signalsCPLW = dba.peakset(K27ac_dba, bRetrieve = T, DataType = DBA_DATA_FRAME) # 474060 enhancers
rownames(K27ac_enhancer_signalsCPLW) = paste0(K27ac_enhancer_signalsCPLW$CHR, ":", K27ac_enhancer_signalsCPLW$START, "-", K27ac_enhancer_signalsCPLW$END)
# divide by log enhancer widths
K27ac_enhancer_signalsCPLW = t(apply(K27ac_enhancer_signalsCPLW, 1, function(x){
  enh_width = as.numeric(x["END"]) - as.numeric(x["START"]) + 1
  CPLW = as.numeric(x[4:ncol(K27ac_enhancer_signalsCPLW)])
  CPLW = sapply(CPLW, function(y){max(y,0)}) # Don't use pseudocount 1! Because that will set smaller enhancers to higher activity, among enhancers with no activity!
  CPLW = CPLW / log2(enh_width)
  CPLW
}))
K27ac_enhancer_signalsCPLW = data.frame(K27ac_enhancer_signalsCPLW)
# rename samples
colnames(K27ac_enhancer_signalsCPLW) = as.character(dba.show(K27ac_dba)$ID)
# rank transform
K27ac_enhancer_signalsCPLW = data.frame(apply(K27ac_enhancer_signalsCPLW, 2, function(x){rank(x, ties.method = "average")}), check.names = F)
K27ac_enhancer_signalsCPLW$enhancer = rownames(K27ac_enhancer_signalsCPLW)



# # ### Don't use TADs, because after merging TADs from different cells we get very large TADs --> 50% of enhancer targets > 1MB away!
# # read enhancer tads, for enhancers without TADs, set TAD to +/-1Mb
# enhancer_tads = read.table("./target_gene_mapping/enhancers_tads.bed", sep="\t", stringsAsFactors = F, header=F) # 80426 enhancers x 4
# colnames(enhancer_tads) = c("chr", "start", "end", "enhancer")
# enhancer_no_tad = K27ac_enhancer_signals[!K27ac_enhancer_signals$enhancer %in% enhancer_tads$enhancer, "enhancer"] #1132
# chr_lens = read.csv("./target_gene_mapping/human.hg19.chrlens", sep="\t", header=F, row.names = 1)
# enhancer_no_tad_enh = t(sapply(enhancer_no_tad, function(x){
#   xsplit = strsplit(x, split=":|-")
#   chr = xsplit[[1]][1]
#   start = as.numeric(xsplit[[1]][2])
#   end = as.numeric(xsplit[[1]][3])
#   tadstart = max(0, start - 1000000)
#   tadend = min(chr_lens[chr,1], end+1000000)
#   c(chr, tadstart, tadend, x)
# }))
# enhancer_no_tad_enh = data.frame(enhancer_no_tad_enh, stringsAsFactors = F)
# colnames(enhancer_no_tad_enh) = c("chr", "start", "end", "enhancer")
# rownames(enhancer_no_tad_enh) = NULL
# enhancer_no_tad_enh$start = as.numeric(enhancer_no_tad_enh$start)
# enhancer_no_tad_enh$end = as.numeric(enhancer_no_tad_enh$end) # 1132 rows
# enhancer_tads = rbind(enhancer_tads, enhancer_no_tad_enh) # 81558 rows
# length(intersect(enhancer_tads$enhancer, K27ac_enhancer_signals$enhancer)) == length(enhancer_tads$enhancer)
# rm(chr_lens)
# rm(enhancer_no_tad_enh)
# rm(enhancer_no_tad)


# just set each enhancer's range to +/- 1Mb
enhancer_no_tad = K27ac_enhancer_signalsCPLW$enhancer # 476325
chr_lens = read.csv("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/target_gene_mapping/human.hg19.chrlens", sep="\t", header=F, row.names = 1)
enhancer_no_tad_enh = t(sapply(enhancer_no_tad, function(x){
  xsplit = strsplit(x, split=":|-")
  chr = xsplit[[1]][1]
  start = as.numeric(xsplit[[1]][2])
  end = as.numeric(xsplit[[1]][3])
  tadstart = max(0, start - 1000000)
  tadend = min(chr_lens[chr,1], end+1000000)
  c(chr, tadstart, tadend, x)
}))
enhancer_no_tad_enh = data.frame(enhancer_no_tad_enh, stringsAsFactors = F)
colnames(enhancer_no_tad_enh) = c("chr", "start", "end", "enhancer")
rownames(enhancer_no_tad_enh) = NULL
enhancer_no_tad_enh$start = as.numeric(enhancer_no_tad_enh$start)
enhancer_no_tad_enh$end = as.numeric(enhancer_no_tad_enh$end)
enhancer_tads = enhancer_no_tad_enh # 476325 rows
length(intersect(enhancer_tads$enhancer, K27ac_enhancer_signalsTMM$enhancer)) == length(enhancer_tads$enhancer)
length(intersect(enhancer_tads$enhancer, K27ac_enhancer_signalsRPKM$enhancer)) == length(enhancer_tads$enhancer)
length(intersect(enhancer_tads$enhancer, K27ac_enhancer_signalsCPLW$enhancer)) == length(enhancer_tads$enhancer)
rm(chr_lens)
rm(enhancer_no_tad_enh)
rm(enhancer_no_tad)



# sort enhancer signals to have same row ordering as enhancer tads
K27ac_enhancer_signalsTMM = merge(enhancer_tads, K27ac_enhancer_signalsTMM, by="enhancer", sort=F)
K27ac_enhancer_signalsRPKM = merge(enhancer_tads, K27ac_enhancer_signalsRPKM, by="enhancer", sort=F)
K27ac_enhancer_signalsCPLW = merge(enhancer_tads, K27ac_enhancer_signalsCPLW, by="enhancer", sort=F)
all.equal(K27ac_enhancer_signalsTMM$enhancer, enhancer_tads$enhancer)
all.equal(K27ac_enhancer_signalsRPKM$enhancer, enhancer_tads$enhancer)
all.equal(K27ac_enhancer_signalsCPLW$enhancer, enhancer_tads$enhancer)
K27ac_enhancer_signalsTMM = K27ac_enhancer_signalsTMM[, !colnames(K27ac_enhancer_signalsTMM) %in% c("chr","start","end")]
K27ac_enhancer_signalsRPKM = K27ac_enhancer_signalsRPKM[, !colnames(K27ac_enhancer_signalsRPKM) %in% c("chr","start","end")]
K27ac_enhancer_signalsCPLW = K27ac_enhancer_signalsCPLW[, !colnames(K27ac_enhancer_signalsCPLW) %in% c("chr","start","end")]





###### read RNAseq tpm ######
# Ensure samples have same IDs as chipseq data!
rnatpm = read.table("expression/counts_matrix_tissue_CCA_HCC/tpm_matrix-genes.txt", stringsAsFactors = F, header=T, check.names = F, sep="\t") # 24622 rows, all unique gene names.
colnames(rnatpm)[colnames(rnatpm)=="gene_name"] = "gene"
rnatpm = rnatpm[, colnames(rnatpm)!="gene_id"]

# keep only samples with both enhancer and RNAseq data
K27ac_enhancer_signalsTMM_allsamps = K27ac_enhancer_signalsTMM # keep this, to calculate corr with promoter K27ac later on
K27ac_enhancer_signalsRPKM_allsamps = K27ac_enhancer_signalsRPKM # keep this, to calculate corr with promoter K27ac later on
K27ac_enhancer_signalsCPLW_allsamps = K27ac_enhancer_signalsCPLW # keep this, to calculate corr with promoter K27ac later on
samps_keep = intersect(colnames(K27ac_enhancer_signalsCPLW_allsamps)[!colnames(K27ac_enhancer_signalsCPLW_allsamps) %in% c("enhancer")], colnames(rnatpm)[!colnames(rnatpm) %in% c("gene")]) # 71: 67 samples in excel (K27ac non-Oct2018 + RNAseq) + 3 HCC + 807N
rnatpm = rnatpm[, match(c("gene", samps_keep), colnames(rnatpm))]
K27ac_enhancer_signalsTMM = K27ac_enhancer_signalsTMM[, match(c("enhancer", samps_keep), colnames(K27ac_enhancer_signalsTMM))]
K27ac_enhancer_signalsRPKM = K27ac_enhancer_signalsRPKM[, match(c("enhancer", samps_keep), colnames(K27ac_enhancer_signalsRPKM))]
K27ac_enhancer_signalsCPLW = K27ac_enhancer_signalsCPLW[, match(c("enhancer", samps_keep), colnames(K27ac_enhancer_signalsCPLW))]



######## Get RNAseq signal for each promoters ######
# get genes and promoters
promoter_regions = read.table("target_gene_mapping/genes_refseq_curated_onlyNM.tss.filt.sort.txt",stringsAsFactors = F, header=F, sep="\t")# 25850, 19225 unique genes.
colnames(promoter_regions) = c("chr","start","gene")
promoter_regions$end = promoter_regions$start
# # promoter_regions$start = promoter_regions$start+1
promoter_regions$promoter_id = paste0(promoter_regions$chr, ":", promoter_regions$start, "-", promoter_regions$end, "_", promoter_regions$gene)
promoter_regions$coord = paste0(promoter_regions$chr, ":", promoter_regions$start, "-", promoter_regions$end)
promoter_regions = promoter_regions[, c("coord", "gene", "promoter_id", "chr", "start", "end")]

# create promoter rnaseq data
promoter_rnatpm = merge(promoter_regions, rnatpm, by="gene", sort=F) # 25820 rows, 19195 unique genes
rm(promoter_regions)
rm(rnatpm)


# remove sex chr, genes expressed in <5% of samples, rank transform
promoter_rnatpm = promoter_rnatpm[!promoter_rnatpm$chr %in% c("chrX", "chrY"),] # 24658, 18331 genes
promoter_rnatpm_expr = apply(promoter_rnatpm[,!colnames(promoter_rnatpm) %in% c("coord","gene","promoter_id","chr","start","end")], 1, function(x){sum(x>.5)}) # .5 cutoff
promoter_rnatpm_expr = promoter_rnatpm_expr / sum(!colnames(promoter_rnatpm) %in% c("coord","gene","promoter_id","chr","start","end")) > .05 # 22156
promoter_rnatpm = promoter_rnatpm[promoter_rnatpm_expr,] # 22156, 16078 genes
promoter_rnatpm[,!colnames(promoter_rnatpm) %in% c("coord","gene","promoter_id","chr","start","end")] = apply(promoter_rnatpm[,!colnames(promoter_rnatpm) %in% c("coord","gene","promoter_id","chr","start","end")], 2, function(x){rank(x, ties.method = "average")})
rm(promoter_rnatpm_expr)




# check that enhancer signals have same sample ordering as promoter signals
all.equal(colnames(K27ac_enhancer_signalsTMM)[colnames(K27ac_enhancer_signalsTMM)!="enhancer"], colnames(promoter_rnatpm)[!colnames(promoter_rnatpm) %in% c("coord", "gene", "chr", "start", "end", "promoter_id")])
all.equal(colnames(K27ac_enhancer_signalsRPKM)[colnames(K27ac_enhancer_signalsRPKM)!="enhancer"], colnames(promoter_rnatpm)[!colnames(promoter_rnatpm) %in% c("coord", "gene", "chr", "start", "end", "promoter_id")])
all.equal(colnames(K27ac_enhancer_signalsCPLW)[colnames(K27ac_enhancer_signalsCPLW)!="enhancer"], colnames(promoter_rnatpm)[!colnames(promoter_rnatpm) %in% c("coord", "gene", "chr", "start", "end", "promoter_id")])

# check that enhancer tads have same row ordering as enhancer signals
all.equal(enhancer_tads$enhancer, K27ac_enhancer_signalsTMM$enhancer)
all.equal(enhancer_tads$enhancer, K27ac_enhancer_signalsRPKM$enhancer)
all.equal(enhancer_tads$enhancer, K27ac_enhancer_signalsCPLW$enhancer)


# get promoters for each enhancer
library(GenomicRanges)
enhancer_tads_gr = GRanges(seqnames=Rle(enhancer_tads$chr), ranges=IRanges(start=enhancer_tads$start+1, end=enhancer_tads$end), strand=Rle("+"), enhancer=enhancer_tads$enhancer)
promoters_gr = GRanges(seqnames=Rle(promoter_rnatpm$chr), ranges=IRanges(start=promoter_rnatpm$start, end=promoter_rnatpm$end), strand=Rle("+"), gene=promoter_rnatpm$gene)
enhancer_tads_promoters = data.frame(findOverlaps(enhancer_tads_gr, promoters_gr, ignore.strand=T)) # 10616644 overlaps
colnames(enhancer_tads_promoters) = c("enhancer_idx", "promoter_idx")
rm(promoters_gr)
rm(enhancer_tads_gr)

# save promoter and enhancer data to calculate correlation on cluster
promoter_data = as.matrix(promoter_rnatpm[, !colnames(promoter_rnatpm) %in% c("coord", "gene", "chr", "start", "end", "promoter_id")]) # 22156 promoters x 71
enhancer_data = as.matrix(K27ac_enhancer_signalsTMM[, colnames(K27ac_enhancer_signalsTMM)!="enhancer"]) #  474060 enhancers x 71
all.equal(colnames(promoter_data), colnames(enhancer_data))
save(promoter_data, enhancer_data, enhancer_tads_promoters, file="target_gene_mapping/enhancerTMM_promoter_input.RData")

enhancer_data = as.matrix(K27ac_enhancer_signalsRPKM[, colnames(K27ac_enhancer_signalsRPKM)!="enhancer"]) #  474060 enhancers x 71
all.equal(colnames(promoter_data), colnames(enhancer_data))
save(promoter_data, enhancer_data, enhancer_tads_promoters, file="target_gene_mapping/enhancerRPKM_promoter_input.RData")

enhancer_data = as.matrix(K27ac_enhancer_signalsCPLW[, colnames(K27ac_enhancer_signalsCPLW)!="enhancer"]) #  474060 enhancers x 71
all.equal(colnames(promoter_data), colnames(enhancer_data))
save(promoter_data, enhancer_data, enhancer_tads_promoters, file="target_gene_mapping/enhancerCPLW_promoter_input.RData")


# too big! run on cluster
# zzz = t(apply(enhancer_tads_promoters, 1, function(x){
#   enhancer_signal = enhancer_data[as.numeric(x["enhancer_idx"]), ]
#   promoter_signal = promoter_data[as.numeric(x["promoter_idx"]), ]
#   corresult = cor.test(enhancer_signal, promoter_signal, method="spearman")
#   c(corresult$estimate, corresult$p.value)
# }))
# colnames(zzz) = c("rho", "p")
# enhancer_tads_promoters = cbind(enhancer_tads_promoters, zzz)




# upload to server CCA_ChipSeq/processed/target_gene_mapping/calc_coef/
# run on server:
mkdir cor_output_TMM
for i in {1..10616644..1000000}
do
  end_idx=$((i + 1000000 - 1))
  echo $i $end_idx
  qsub /gpfs/apps/R/3.1.0/bin/Rscript --vanilla calc_coef.R enhancerTMM_promoter_input.RData cor_output_TMM/enhancerTMM_promoter_cor_${i}.txt $i $end_idx
done
mkdir cor_output_RPKM
for i in {1..10616644..1000000}
do
  end_idx=$((i + 1000000 - 1))
  echo $i $end_idx
  qsub /gpfs/apps/R/3.1.0/bin/Rscript --vanilla calc_coef.R enhancerRPKM_promoter_input.RData cor_output_RPKM/enhancerRPKM_promoter_cor_${i}.txt $i $end_idx
done
mkdir cor_output_CPLW
for i in {1..10616644..1000000}
do
  end_idx=$((i + 1000000 - 1))
  echo $i $end_idx
  qsub /gpfs/apps/R/3.1.0/bin/Rscript --vanilla calc_coef.R enhancerCPLW_promoter_input.RData cor_output_CPLW/enhancerCPLW_promoter_cor_${i}.txt $i $end_idx
done

# combine all outputs
rm enhancerTMM_promoter_cor.txt
for i in {1..10616644..1000000}
do
  cat cor_output_TMM/enhancerTMM_promoter_cor_${i}.txt >> enhancerTMM_promoter_cor.txt
done
rm enhancerRPKM_promoter_cor.txt
for i in {1..10616644..1000000}
do
  cat cor_output_RPKM/enhancerRPKM_promoter_cor_${i}.txt >> enhancerRPKM_promoter_cor.txt
done
rm enhancerCPLW_promoter_cor.txt
for i in {1..10616644..1000000}
do
  cat cor_output_CPLW/enhancerCPLW_promoter_cor_${i}.txt >> enhancerCPLW_promoter_cor.txt
done



# Download enhancer_promoter_cor.txt from server, load calculated correlations
corr_data_TMM = read.table("target_gene_mapping/enhancerTMM_promoter_cor.txt",  sep="\t", stringsAsFactors = F, header=F)  #  10616644 rows
corr_data_RPKM = read.table("target_gene_mapping/enhancerRPKM_promoter_cor.txt",  sep="\t", stringsAsFactors = F, header=F)  #  10616644 rows
corr_data_CPLW = read.table("target_gene_mapping/enhancerCPLW_promoter_cor.txt",  sep="\t", stringsAsFactors = F, header=F)  #  10616644 rows
colnames(corr_data_TMM) = c("enhancer_idx", "promoter_idx", "rho", "p")
colnames(corr_data_RPKM) = c("enhancer_idx", "promoter_idx", "rho", "p")
colnames(corr_data_CPLW) = c("enhancer_idx", "promoter_idx", "rho", "p")


# test a few to check results
corr_data_TMM[1,]
cor.test(as.numeric(K27ac_enhancer_signalsTMM[1,2:ncol(K27ac_enhancer_signalsTMM)]), as.numeric(promoter_rnatpm[1, 7:ncol(promoter_rnatpm)]), method="spearman")
corr_data_TMM[5000000,]
cor.test(as.numeric(K27ac_enhancer_signalsTMM[199444,2:ncol(K27ac_enhancer_signalsTMM)]), as.numeric(promoter_rnatpm[18738, 7:ncol(promoter_rnatpm)]), method="spearman")
corr_data_TMM[10616644,]
cor.test(as.numeric(K27ac_enhancer_signalsTMM[474060,2:ncol(K27ac_enhancer_signalsTMM)]), as.numeric(promoter_rnatpm[11143, 7:ncol(promoter_rnatpm)]), method="spearman")

corr_data_RPKM[1,]
cor.test(as.numeric(K27ac_enhancer_signalsRPKM[1,2:ncol(K27ac_enhancer_signalsRPKM)]), as.numeric(promoter_rnatpm[1, 7:ncol(promoter_rnatpm)]), method="spearman")
corr_data_RPKM[5000000,]
cor.test(as.numeric(K27ac_enhancer_signalsRPKM[199444,2:ncol(K27ac_enhancer_signalsRPKM)]), as.numeric(promoter_rnatpm[18738, 7:ncol(promoter_rnatpm)]), method="spearman")
corr_data_RPKM[10616644,]
cor.test(as.numeric(K27ac_enhancer_signalsRPKM[474060,2:ncol(K27ac_enhancer_signalsRPKM)]), as.numeric(promoter_rnatpm[11143, 7:ncol(promoter_rnatpm)]), method="spearman")

corr_data_CPLW[1,]
cor.test(as.numeric(K27ac_enhancer_signalsCPLW[1,2:ncol(K27ac_enhancer_signalsCPLW)]), as.numeric(promoter_rnatpm[1, 7:ncol(promoter_rnatpm)]), method="spearman")
corr_data_CPLW[5000000,]
cor.test(as.numeric(K27ac_enhancer_signalsCPLW[199444,2:ncol(K27ac_enhancer_signalsCPLW)]), as.numeric(promoter_rnatpm[18738, 7:ncol(promoter_rnatpm)]), method="spearman")
corr_data_CPLW[10616644,]
cor.test(as.numeric(K27ac_enhancer_signalsCPLW[474060,2:ncol(K27ac_enhancer_signalsCPLW)]), as.numeric(promoter_rnatpm[11143, 7:ncol(promoter_rnatpm)]), method="spearman")

plot(corr_data_RPKM[1:10000,"rho"], corr_data_CPLW[1:10000,"rho"])
plot(corr_data_TMM[1:10000,"rho"], corr_data_CPLW[1:10000,"rho"])



##### Get target mapping from correlation data

## calculate enhancer targets from corr_data
# corr_data indices refer to enhancer_tads and promoter_rnatpm
# K27ac_enh and K27ac_prom are only used to calculate enhancer-promoter K27ac correlation
get_target_map = function(corr_data, enhancer_tads, promoter_rnatpm, K27ac_enh, K27ac_prom, criteria="normal", maxdist=100000000, rho_thr=.5) {

  corr_data$enhancer = enhancer_tads[corr_data$enhancer_idx, "enhancer"]
  corr_data$promoter = promoter_rnatpm[corr_data$promoter_idx, "promoter_id"]

  # # adjust p-values per enhancer
  # corr_data_split = split(corr_data, corr_data$enhancer_idx)
  # corr_data_split_padj = lapply(corr_data_split, function(x){
  #   cbind(as.matrix(x), p.adjust(x$p, method="BH"))
  # })
  # corr_data = do.call(rbind, corr_data_split_padj)
  # rm(corr_data_split)
  # rm(corr_data_split_padj)
  # corr_data = data.frame(corr_data) # 3106605 x 5
  # colnames(corr_data)[5] = c("padj")
  # rownames(corr_data) = NULL

  # keep only rows with p < .05 & rho >= rho_thr
  corr_data = corr_data[corr_data$p<.05 & corr_data$rho>=rho_thr,]
  rownames(corr_data) = NULL

  # for each enhancer, best_target is the target with highest rho, where p<.05 and rho>=rho_thr.
  # Also keep additional targets whose rho are within 0.05 of best_target, and p<.05 and rho>=rho_thr.
  if (criteria=="normal") {
    corr_data_split = split(corr_data, corr_data$enhancer_idx)
    zzz = lapply(corr_data_split, function(x){
      x$is_target = 0
      max_rho_idx = which.max(x$rho)
      if (x$p[max_rho_idx] < .05 & x$rho[max_rho_idx] >= rho_thr) {
        # set best target
        x$is_target[max_rho_idx] = 1
        # set additional targets
        x[x$p < .05 & x$rho >= rho_thr & x$rho >= x$rho[max_rho_idx] - 0.05, "is_target"] = 1
      }
      x
    })
    corr_data = do.call(rbind, zzz)
    rm(zzz)
    rm(corr_data_split)

    # More liberal criteria: for each enhancer, targets are those with p<.05 and rho >= rho_thr
  } else if (criteria=="liberal") {
    corr_data_split = split(corr_data, corr_data$enhancer_idx)
    zzz = lapply(corr_data_split, function(x){
      x$is_target = 0
      x[x$rho >= rho_thr & x$p < .05, "is_target"] = 1
      x
    })
    corr_data = do.call(rbind, zzz)
    rm(zzz)
    rm(corr_data_split)
  }
  rownames(corr_data) = NULL
  target_map = corr_data[corr_data$is_target==1,]


  # calculate distances. If a TSS is inside enhancer, set dist to 0
  target_map$dist =  apply(target_map, 1, function(x){
    tss = as.numeric(strsplit(as.character(x["promoter"]), split="[:_-]", perl = T)[[1]][2])
    enh_start = as.numeric(strsplit(as.character(x["enhancer"]), split="[:_-]", perl = T)[[1]][2])
    enh_end = as.numeric(strsplit(as.character(x["enhancer"]), split="[:_-]", perl = T)[[1]][3])
    if (enh_start <= tss & tss <= enh_end) { 0
    } else {    min(abs(enh_start-tss), abs(enh_end-tss)) }
  })
  target_map = target_map[target_map$dist <= maxdist,]

  # extract genes
  target_map$gene = sapply(strsplit(target_map$promoter, split="_"), function(x){x[2]})
  target_map = target_map[, c("enhancer", "gene", "promoter", "rho", "p", "dist")]
  rownames(target_map) = NULL

  # calculate enhancer-promoter K27ac correlation
  if (!is.null(K27ac_prom)) {
    # keep only relevant data in enhancer and promoter signals
    rownames(K27ac_enh) = K27ac_enh$enhancer
    K27ac_enh = K27ac_enh[, !colnames(K27ac_enh) %in% "enhancer"]
    rownames(K27ac_prom) = K27ac_prom$promoter_id
    K27ac_prom = K27ac_prom[, !colnames(K27ac_prom) %in% c("coord","gene","promoter_id","chr","start","end")]

    # check that enhancer signals have same sample ordering as promoter signals
    if (!all.equal(colnames(K27ac_enh), colnames(K27ac_prom))) { return(NULL)}

    # calculate correlation between enhancer and promoter K27ac
    zz = apply(target_map, 1, function(x){
      enh_sig = as.numeric(K27ac_enh[as.character(x["enhancer"]), ])
      prom_sig = as.numeric(K27ac_prom[as.character(x["promoter"]), ])
      corcoef = cor.test(enh_sig, prom_sig, method="spearman")
      data.frame(promK27ac_rho=corcoef$estimate, promK27ac_p=corcoef$p.value)
    })
    zz = do.call(rbind, zz)
    target_map = cbind(target_map, zz)
    rownames(target_map) = NULL
  } else {
    target_map = cbind(target_map, promK27ac_rho=as.numeric(NA), promK27ac_p=as.numeric(NA))
  }

  ### append enhancers without targets
  enh_no_targets = enhancer_tads$enhancer[!enhancer_tads$enhancer %in% target_map$enhancer]
  enh_no_targets = data.frame(enhancer=enh_no_targets, gene=as.character(NA), promoter=as.character(NA), rho=as.numeric(NA), p=as.numeric(NA), dist=as.numeric(NA), promK27ac_rho=as.numeric(NA), promK27ac_p=as.numeric(NA), stringsAsFactors = F)
  target_map = rbind(target_map, enh_no_targets) # 181078 rows, 140385 enhancers
  return(target_map)
}



# load promoter K27ac data for validation (check if enhancer and promoter K27ac are correlated), rank transform
# convert the promoter K27ac promoter ids from promoter regions to TSS
load("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/Rdata_cca_k27ac_1_targetmapping_promotersignal.RData")
rm(K27ac_dba_promoter_split)
K27ac_promoter[,!colnames(K27ac_promoter) %in% c("coord","gene","promoter_id","chr","start","end")] = apply(K27ac_promoter[,!colnames(K27ac_promoter) %in% c("coord","gene","promoter_id","chr","start","end")], 2, function(x){rank(x, ties.method = "average")})
K27ac_promoter$promoter_id = sapply(strsplit(K27ac_promoter$coord, split="[:-]"), function(x){
  chr = as.character(x[1])
  pos = mean(as.numeric(c(x[2], x[3]))) - 1
  paste0(chr, ":", pos, "-", pos)
})
K27ac_promoter$promoter_id = paste0(K27ac_promoter$promoter_id, "_", K27ac_promoter$gene)

### get target mapping: rho.5
target_map_TMM_rho.5 = get_target_map(corr_data=corr_data_TMM, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsTMM_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.5) # 17218/474060 enhancers mapped to  4748 genes
target_map_RPKM_rho.5 = get_target_map(corr_data=corr_data_RPKM, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsRPKM_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.5) # 36524/474060 enhancers mapped to 6990 genes
target_map_CPLW_rho.5 = get_target_map(corr_data=corr_data_CPLW, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsCPLW_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.5) #  18069/474060 enhancers mapped to 4779 genes
target_map_CPLW_rho.5_liberal = get_target_map(corr_data=corr_data_CPLW, enhancer_tads, promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsCPLW_allsamps, K27ac_prom=K27ac_promoter, criteria="liberal", maxdist=1000000, rho_thr=.5) #  18069/474060 enhancers mapped to 5019 genes

### get target mapping: rho.4
target_map_TMM_rho.4 = get_target_map(corr_data=corr_data_TMM, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsTMM_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.4) # 67653/474060 enhancers mapped to 10870 genes
target_map_RPKM_rho.4 = get_target_map(corr_data=corr_data_RPKM, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsRPKM_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.4) # 123074/474060 enhancers mapped to 12999 genes
target_map_CPLW_rho.4 = get_target_map(corr_data=corr_data_CPLW, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsCPLW_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.4) # 68258/474060 enhancers mapped to 10992 genes
target_map_CPLW_rho.4_liberal = get_target_map(corr_data=corr_data_CPLW, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsCPLW_allsamps, K27ac_prom=K27ac_promoter, criteria="liberal", maxdist=1000000, rho_thr=.4) # 68258/474060 enhancers mapped to 11613 genes
write.table(target_map_TMM_rho.4, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_TMM_rho.4.txt", sep="\t", quote = F, row.names=F, col.names=T)
write.table(target_map_RPKM_rho.4, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_RPKM_rho.4.txt", sep="\t", quote = F, row.names=F, col.names=T)
write.table(target_map_CPLW_rho.4, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_CPLW_rho.4.txt", sep="\t", quote = F, row.names=F, col.names=T)
write.table(target_map_CPLW_rho.4_liberal, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_CPLW_rho.4_liberal.txt", sep="\t", quote = F, row.names=F, col.names=T)


### get target mapping: rho.3
target_map_TMM_rho.3 = get_target_map(corr_data=corr_data_TMM, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsTMM_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.3) # 191912/474060 enhancers mapped to 15068 genes
target_map_RPKM_rho.3 = get_target_map(corr_data=corr_data_RPKM, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsRPKM_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.3) # 273022/474060 enhancers mapped to 15627 genes
target_map_CPLW_rho.3 = get_target_map(corr_data=corr_data_CPLW, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsCPLW_allsamps, K27ac_prom=K27ac_promoter, maxdist=1000000, rho_thr=.3) # 189514/474060 enhancers mapped to 15235 genes
target_map_CPLW_rho.3_liberal = get_target_map(corr_data=corr_data_CPLW, enhancer_tads=enhancer_tads, promoter_rnatpm=promoter_rnatpm, K27ac_enh=K27ac_enhancer_signalsCPLW_allsamps, K27ac_prom=K27ac_promoter, criteria="liberal", maxdist=1000000, rho_thr=.3) # 189514/474060 enhancers mapped to 15602 genes
write.table(target_map_TMM_rho.3, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_TMM_rho.3.txt", sep="\t", quote = F, row.names=F, col.names=T)
write.table(target_map_RPKM_rho.3, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_RPKM_rho.3.txt", sep="\t", quote = F, row.names=F, col.names=T)
write.table(target_map_CPLW_rho.3, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_CPLW_rho.3.txt", sep="\t", quote = F, row.names=F, col.names=T)
write.table(target_map_CPLW_rho.3_liberal, file="target_gene_mapping/20190628_enhancer_target_map_noTAD_CPLW_rho.3_liberal.txt", sep="\t", quote = F, row.names=F, col.names=T)
