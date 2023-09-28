
########### Get Deseq2 normalization factors. If samps==NULL, use all samples in my_dba
get_deseq2_normfacs = function(my_dba, samps=NULL, debugPlot=F){

  # get counts directly. Assume that control reads are scaled!
  my_dba2 = dba.count(my_dba, peaks=NULL, score=DBA_SCORE_READS_MINUS)
  zz = dba.peakset(my_dba2, bRetrieve = T, DataType = DBA_DATA_FRAME)
  zz = zz[, 4:ncol(zz)]
  zz = apply(zz, 2, function(x){sapply(x, function(y){max(y,1)})})
  colnames(zz) = as.character(dba.show(my_dba)$ID)
  if (!is.null(samps)) {
    zz = zz[, colnames(zz) %in% samps]
  }
  apply(zz, 2, summary)

  ### calculate pseudosample: geometric mean of each enhancer
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  zz_gm_mean = apply(zz, 1, gm_mean)

  ### for each sample, calculate ratio of counts to pseudosample. for each sample, get the median of ratios as the normalization factor
  zz_ratio = apply(zz, 2, function(x){x/zz_gm_mean})
  apply(zz_ratio, 2, summary)

  if (debugPlot==T) {
    plot_sample = "21914187"
    plot_ratios = sort(zz_ratio[, plot_sample])
    plot_ratios = plot_ratios[seq(from=1, to=length(plot_ratios), by=10)]
    plot(plot_ratios, 1:length(plot_ratios), pch=20, cex=.5, xlim=c(0,50))
    abline(h=length(plot_ratios)/2)
    abline(v=median(plot_ratios))

    plot_sample = "29150572"
    plot_ratios = sort(zz_ratio[, plot_sample])
    plot_ratios = plot_ratios[seq(from=1, to=length(plot_ratios), by=10)]
    points(plot_ratios, 1:length(plot_ratios), pch=20, cex=.5, col="red")
    abline(v=median(plot_ratios), col="red")
    abline(h=length(plot_ratios)/2)

    plot_sample = "A100"
    plot_ratios = sort(zz_ratio[, plot_sample])
    plot_ratios = plot_ratios[seq(from=1, to=length(plot_ratios), by=10)]
    points(plot_ratios, 1:length(plot_ratios), pch=20, cex=.5, col="green")
    abline(v=median(plot_ratios), col="green")
    abline(h=length(plot_ratios)/2)

    plot_sample = "807N"
    plot_ratios = sort(zz_ratio[, plot_sample])
    plot_ratios = plot_ratios[seq(from=1, to=length(plot_ratios), by=10)]
    points(plot_ratios, 1:length(plot_ratios), pch=20, cex=.5, col="blue")
    abline(v=median(plot_ratios), col="blue")
    abline(h=length(plot_ratios)/2)

  }

  ### get norm facs
  normfacs = apply(zz_ratio, 2, median)
  # names(normfacs) = as.character(dba.show(my_dba2)$ID)
  return(normfacs)
}


######## Create list of col_data used to run differential expression analysis, corresponding to diffbind contrasts
# Input: list of contrasts for differential binding
# Output: list of col_data for differential expression analysis
create_diffexpr_contrasts = function(dba_contrasts, expr_samples) {
  col_data_list = lapply(dba_contrasts, function(x){
    grp1_samps = names(x$group1)[x$group1]
    grp2_samps = names(x$group2)[x$group2]
    col_data = data.frame(sample = c(grp1_samps,grp2_samps), stringsAsFactors=F)
    col_data$contrast = factor(c(rep(x$name1, length(grp1_samps)), rep(x$name2, length(grp2_samps))), levels=c(x$name2, x$name1))
    col_data = col_data[col_data$sample %in% expr_samples,]
  })
  return (col_data_list)
}


######### Run differential expression analyses of a lis of col_Data
run_diffexpr_contrasts = function(my_col_data_list, my_expr_counts_data) {
  diffexpr_list = lapply(1:length(my_col_data_list), function(x){
    if (length(unique(my_col_data_list[[x]]$contrast))<=1) {
      list(diffexpr=NULL, sizeFactors=NULL)
    }
    else {
      dds = DESeqDataSetFromMatrix(countData = my_expr_counts_data[, my_col_data_list[[x]]$sample], colData = my_col_data_list[[x]], design = ~ contrast)
      dds = DESeq(dds)
      diffexpr = results(dds, name=paste0("contrast_", levels(my_col_data_list[[x]]$contrast)[2], "_vs_", levels(my_col_data_list[[x]]$contrast)[1]))
      diffexpr = diffexpr[order(diffexpr$stat, decreasing = T),]
      list(diffexpr=diffexpr, sizeFactors=dds$sizeFactor)
    }
  })
}



##### Compare differential binding results from 4 possible methods:
# 1. Deseq2, effective library size (ie. normalize by its median ratio method)
# 2. Deseq2, full library size (ie. normalize by library size only)
# 3. EdgeR, effective library size (ie. normalize by TMM effective library size)
# 4. EdgeR, full library size (ie. normalize by TMM full library size)
# filterEnh: only consider enhancers in this list. This is used when analyzing superenhancers diffbind, where I gave it enhancers+superenhancers so that normalization factors won't be based on superenhancers.
# Note that when merging enhancers with target map, enhancers not in corr_data will be removed anyway.
# When filterEnh is not NULL, adjusted P-values will be recalculated.
compare_diff_enhancers_methods = function(my_H3K27ac_dba_effLib, my_H3K27ac_dba_fullLib, contrast, corr_data, diffexprlist, my_expr_tpm_data, filterEnh=NULL) {


  # keep only unique enhancer-gene pairs in target map
  corr_data = corr_data[!is.na(corr_data$gene),c("enhancer", "gene")]
  corr_data = corr_data[!duplicated(corr_data),]
  # keep only non-NA diffexpr foldchange
  diffexpr = diffexprlist[[contrast]]
  diffexpr = diffexpr$diffexpr
  diffexpr = data.frame(gene=rownames(diffexpr), expr_fold=diffexpr$log2FoldChange, expr_stat=diffexpr$stat, expr_padj=diffexpr$padj, stringsAsFactors = F)
  diffexpr = diffexpr[!is.na(diffexpr$expr_fold),]


  ########### print some configuration values
  #  Deseq2 eff lib (normalize by median ratio): normalized score ~= counts / fac
  print(paste0("Deseq2, effective library size (ie. normalize by median ratio): score ~= counts / fac"), quote = F)
  print(paste0("FullLib = ", my_H3K27ac_dba_effLib$contrasts[[contrast]]$DESeq2$bFullLibrarySize), quote = F)
  print(paste0("SubControl = ", my_H3K27ac_dba_effLib$contrasts[[contrast]]$DESeq2$bSubControl), quote = F)
  print(paste0("samples = ", paste0(c(names(my_H3K27ac_dba_effLib$contrasts[[contrast]]$group1)[my_H3K27ac_dba_effLib$contrasts[[contrast]]$group1],
                                      names(my_H3K27ac_dba_effLib$contrasts[[contrast]]$group2)[my_H3K27ac_dba_effLib$contrasts[[contrast]]$group2]),collapse=" ")), quote = F)
  print(paste0("Facs = ", paste0(my_H3K27ac_dba_effLib$contrasts[[contrast]]$DESeq2$facs, collapse=" ")), quote = F)

  #  Deseq2 full lib (normalize by library size): normalized score ~= counts / fac
  print(paste0(""), quote = F)
  print(paste0(""), quote = F)
  print(paste0("Deseq2, full library size (ie. normalize by library size): score ~= counts / fac"), quote = F)
  print(paste0("FullLib = ", my_H3K27ac_dba_fullLib$contrasts[[contrast]]$DESeq2$bFullLibrarySize), quote = F)
  print(paste0("SubControl = ", my_H3K27ac_dba_fullLib$contrasts[[contrast]]$DESeq2$bSubControl), quote = F)
  print(paste0("samples = ", paste0(c(names(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group1)[my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group1],
                                      names(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group2)[my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group2]),collapse=" ")), quote = F)
  print(paste0("Facs = ", paste0(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$DESeq2$facs, collapse=" ")), quote = F)

  # EdgeR eff lib (normalize by TMM effective lib size): normalized score ~= counts / (lib.size * norm.factor) * (some large constant factor)
  print(paste0(""), quote = F)
  print(paste0(""), quote = F)
  print(paste0("EdgeR, effective library size (ie. normalize by TMM eff lib size): score ~= counts / (lib.size * norm.factor) * C"), quote = F)
  print(paste0("FullLib = ", my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$bFullLibrarySize), quote = F)
  print(paste0("SubControl = ", my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$bSubControl), quote = F)
  print(paste0("samples = ", paste0(rownames(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples), collapse=" ")), quote = F)
  print(paste0("Facs = ", paste0(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$facs, collapse=" ")), quote = F)
  print(paste0("lib.size = ", paste0(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples$lib.size, collapse=" ")), quote = F)
  print(paste0("norm.factors = ", paste0(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples$norm.factors, collapse=" ")), quote = F)
  print(paste0("lib.size * norm.factors = ", paste0(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples$lib.size * my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples$norm.factors, collapse=" ")), quote = F)

  # EdgeR full lib (normalize by TMM full lib size): normalized score ~= counts / (lib.size * norm.factor) * (some large constant factor)
  print(paste0(""), quote = F)
  print(paste0(""), quote = F)
  print(paste0("EdgeR, full library size (ie. normalize by TMM full lib size): score ~= counts / (lib.size * norm.factor) * C"), quote = F)
  print(paste0("FullLib = ", my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$bFullLibrarySize), quote = F)
  print(paste0("SubControl = ", my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$bSubControl), quote = F)
  print(paste0("samples = ", paste0(rownames(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples), collapse=" ")), quote = F)
  print(paste0("lib.size = ", paste0(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples$lib.size, collapse=" ")), quote = F)
  print(paste0("norm.factors = ", paste0(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples$norm.factors, collapse=" ")), quote = F)
  print(paste0("lib.size * norm.factors = ", paste0(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples$lib.size * my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples$norm.factors, collapse=" ")), quote = F)

  ################# Plot normalization factors
  layout(mat=matrix(c(1:4), 2, 2, byrow = T))
  # deseq
  barplot(my_H3K27ac_dba_effLib$contrasts[[contrast]]$DESeq2$facs, main="Deseq2, Lib=Eff", names.arg=rownames(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples), las=2, cex.names=.6, ylab="fac")
  # deseq
  barplot(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$DESeq2$facs, main="Deseq2, Lib=Full", names.arg=rownames(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples), las=2, cex.names=.6, ylab="fac")
  # edgeR
  barplot(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples$lib.size * my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples$norm.factors, main="edgeR, Lib=Eff", names.arg=rownames(my_H3K27ac_dba_effLib$contrasts[[contrast]]$edgeR$samples), las=2, cex.names=.6, ylab="libsize * normfac")
  # edgeR
  barplot(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples$lib.size * my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples$norm.factors, main="edgeR, Lib=Full", names.arg=rownames(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$edgeR$samples), las=2, cex.names=.6, ylab="libsize * normfac")


  ##################### MA plot
  layout(mat=matrix(c(1:8), 2, 4, byrow = T))
  # full library size = F
  dba.plotMA(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, bNormalized=F, factor="DESEQ2 norm=F lib=Eff", cex.main=.75)
  dba.plotMA(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, bNormalized=T, factor="DESEQ2 norm=T lib=Eff", cex.main=.75)
  dba.plotMA(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_EDGER, bNormalized=F, factor="EDGER norm=F lib=Eff", cex.main=.75)
  dba.plotMA(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_EDGER, bNormalized=T, factor="EDGER norm=T lib=Eff", cex.main=.75)
  # full library size = T
  dba.plotMA(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, bNormalized=F, factor="DESEQ2 norm=F lib=Full", cex.main=.75)
  dba.plotMA(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, bNormalized=T, factor="DESEQ2 norm=T lib=Full", cex.main=.75)
  dba.plotMA(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_EDGER, bNormalized=F, factor="EDGER norm=F lib=Full", cex.main=.75)
  dba.plotMA(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_EDGER, bNormalized=T, factor="EDGER norm=T lib=Full", cex.main=.75)



  ############## Plot enhancer foldchange vs expression foldchange
  # Using simplistic calculation of expression foldchange
  grp1_samps = names(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group1)[my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group1]
  grp1_samps = intersect(grp1_samps, colnames(my_expr_tpm_data))
  grp2_samps = names(my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group2)[my_H3K27ac_dba_fullLib$contrasts[[contrast]]$group2]
  grp2_samps = intersect(grp2_samps, colnames(my_expr_tpm_data))
  simple_expr_change = my_expr_tpm_data[rownames(my_expr_tpm_data) %in% diffexpr$gene, colnames(my_expr_tpm_data) %in% c(grp1_samps, grp2_samps)]
  simple_expr_change = data.frame(gene=rownames(simple_expr_change),
                                  expr_fold_simple=apply(simple_expr_change, 1, function(x){log2(median(sapply(x[grp1_samps], function(y){max(.5,y)})) / median(sapply(x[grp2_samps], function(y){max(.5,y)})))}), stringsAsFactors = F)


  layout(mat=matrix(c(1:8), 2, 4, byrow = T))
  ## Deseq2, effective library size
  diff_enhancers_deseq_effLib = data.frame(dba.report(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, th=1), stringsAsFactors = F)
  diff_enhancers_deseq_effLib$enhancer = paste0(diff_enhancers_deseq_effLib$seqnames, ":", diff_enhancers_deseq_effLib$start, "-", diff_enhancers_deseq_effLib$end)
  if (!is.null(filterEnh)) {
    diff_enhancers_deseq_effLib = diff_enhancers_deseq_effLib[diff_enhancers_deseq_effLib$enhancer %in% filterEnh,]
    diff_enhancers_deseq_effLib$FDR = p.adjust(diff_enhancers_deseq_effLib$p.value, method="BH")
  }
  diff_enhancers_deseq_effLib = diff_enhancers_deseq_effLib[, c("enhancer", "Fold", "FDR")]
  diff_enhancers_deseq_effLib = merge(diff_enhancers_deseq_effLib, corr_data, by="enhancer", all.x=T, sort=F)
  diff_enhancers_deseq_effLib = merge(diff_enhancers_deseq_effLib, diffexpr, by="gene", sort=F)
  diff_enhancers_deseq_effLib = merge(diff_enhancers_deseq_effLib, simple_expr_change, by="gene", sort=F)
  # Simplistic expression foldchange
  smoothScatter(x=diff_enhancers_deseq_effLib$Fold, y=diff_enhancers_deseq_effLib$expr_fold_simple, xlab="Enhancer log2FC", ylab="Expression log2FC (Simplistic)", main="Enh: Deseq2 lib=Eff. Exp: Simple")
  corcoef = cor.test(diff_enhancers_deseq_effLib$Fold, diff_enhancers_deseq_effLib$expr_fold_simple, method="spearman")
  abline(lm(diff_enhancers_deseq_effLib$expr_fold_simple ~ diff_enhancers_deseq_effLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))
  # Deseq2 expression foldchange
  smoothScatter(x=diff_enhancers_deseq_effLib$Fold, y=diff_enhancers_deseq_effLib$expr_fold, xlab="Enhancer log2FC", ylab="Expression log2FC (Deseq2)", main="Enh: Deseq2 lib=Eff. Exp: Deseq2")
  corcoef = cor.test(diff_enhancers_deseq_effLib$Fold, diff_enhancers_deseq_effLib$expr_fold, method="spearman")
  abline(lm(diff_enhancers_deseq_effLib$expr_fold ~ diff_enhancers_deseq_effLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))

  # EdgeR, effective library size
  diff_enhancers_edger_effLib = data.frame(dba.report(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_EDGER, th=1), stringsAsFactors = F)
  diff_enhancers_edger_effLib$enhancer = paste0(diff_enhancers_edger_effLib$seqnames, ":", diff_enhancers_edger_effLib$start, "-", diff_enhancers_edger_effLib$end)
  if (!is.null(filterEnh)) {
    diff_enhancers_edger_effLib = diff_enhancers_edger_effLib[diff_enhancers_edger_effLib$enhancer %in% filterEnh,]
    diff_enhancers_edger_effLib$FDR = p.adjust(diff_enhancers_edger_effLib$p.value, method="BH")
  }
  diff_enhancers_edger_effLib = diff_enhancers_edger_effLib[, c("enhancer", "Fold", "FDR")]
  diff_enhancers_edger_effLib = merge(diff_enhancers_edger_effLib, corr_data, by="enhancer", all.x=T, sort=F)
  diff_enhancers_edger_effLib = merge(diff_enhancers_edger_effLib, diffexpr, by="gene", sort=F)
  diff_enhancers_edger_effLib = merge(diff_enhancers_edger_effLib, simple_expr_change, by="gene", sort=F)
  # Simplistic expression foldchange
  smoothScatter(x=diff_enhancers_edger_effLib$Fold, y=diff_enhancers_edger_effLib$expr_fold_simple, xlab="Enhancer log2FC", ylab="Expression log2FC (Simplistic)", main="Enh: EdgeR lib=Eff. Exp: Simple")
  corcoef = cor.test(diff_enhancers_edger_effLib$Fold, diff_enhancers_edger_effLib$expr_fold_simple, method="spearman")
  abline(lm(diff_enhancers_edger_effLib$expr_fold_simple ~ diff_enhancers_edger_effLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))
  # Deseq2 expression foldchange
  smoothScatter(x=diff_enhancers_edger_effLib$Fold, y=diff_enhancers_edger_effLib$expr_fold, xlab="Enhancer log2FC", ylab="Expression log2FC (Deseq2)", main="Enh: EdgeR lib=Eff. Exp: Deseq2")
  corcoef = cor.test(diff_enhancers_edger_effLib$Fold, diff_enhancers_edger_effLib$expr_fold, method="spearman")
  abline(lm(diff_enhancers_edger_effLib$expr_fold ~ diff_enhancers_edger_effLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))

  # Deseq2, full library size
  diff_enhancers_deseq_fullLib = data.frame(dba.report(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, th=1), stringsAsFactors = F)
  diff_enhancers_deseq_fullLib$enhancer = paste0(diff_enhancers_deseq_fullLib$seqnames, ":", diff_enhancers_deseq_fullLib$start, "-", diff_enhancers_deseq_fullLib$end)
  if (!is.null(filterEnh)) {
    diff_enhancers_deseq_fullLib = diff_enhancers_deseq_fullLib[diff_enhancers_deseq_fullLib$enhancer %in% filterEnh,]
    diff_enhancers_deseq_fullLib$FDR = p.adjust(diff_enhancers_deseq_fullLib$p.value, method="BH")
  }
  diff_enhancers_deseq_fullLib = diff_enhancers_deseq_fullLib[, c("enhancer", "Fold", "FDR")]
  diff_enhancers_deseq_fullLib = merge(diff_enhancers_deseq_fullLib, corr_data, by="enhancer", all.x=T, sort=F)
  diff_enhancers_deseq_fullLib = merge(diff_enhancers_deseq_fullLib, diffexpr, by="gene", sort=F)
  diff_enhancers_deseq_fullLib = merge(diff_enhancers_deseq_fullLib, simple_expr_change, by="gene", sort=F)
  # Simplistic expression foldchange
  smoothScatter(x=diff_enhancers_deseq_fullLib$Fold, y=diff_enhancers_deseq_fullLib$expr_fold_simple, xlab="Enhancer log2FC", ylab="Expression log2FC (Simplistic)", main="Enh: Deseq2 lib=Full. Exp: Simple")
  corcoef = cor.test(diff_enhancers_deseq_fullLib$Fold, diff_enhancers_deseq_fullLib$expr_fold_simple, method="spearman")
  abline(lm(diff_enhancers_deseq_fullLib$expr_fold_simple ~ diff_enhancers_deseq_fullLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))
  # Deseq2 expression foldchange
  smoothScatter(x=diff_enhancers_deseq_fullLib$Fold, y=diff_enhancers_deseq_fullLib$expr_fold, xlab="Enhancer log2FC", ylab="Expression log2FC (Deseq2)", main="Enh: DESEQ2 lib=Full. Exp: Deseq2")
  corcoef = cor.test(diff_enhancers_deseq_fullLib$Fold, diff_enhancers_deseq_fullLib$expr_fold, method="spearman")
  abline(lm(diff_enhancers_deseq_fullLib$expr_fold ~ diff_enhancers_deseq_fullLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))

  # EdgeR, full library size
  diff_enhancers_edger_fullLib = data.frame(dba.report(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_EDGER, th=1), stringsAsFactors = F)
  diff_enhancers_edger_fullLib$enhancer = paste0(diff_enhancers_edger_fullLib$seqnames, ":", diff_enhancers_edger_fullLib$start, "-", diff_enhancers_edger_fullLib$end)
  if (!is.null(filterEnh)) {
    diff_enhancers_edger_fullLib = diff_enhancers_edger_fullLib[diff_enhancers_edger_fullLib$enhancer %in% filterEnh,]
    diff_enhancers_edger_fullLib$FDR = p.adjust(diff_enhancers_edger_fullLib$p.value, method="BH")
  }
  diff_enhancers_edger_fullLib = diff_enhancers_edger_fullLib[, c("enhancer", "Fold", "FDR")]
  diff_enhancers_edger_fullLib = merge(diff_enhancers_edger_fullLib, corr_data, by="enhancer", all.x=T, sort=F)
  diff_enhancers_edger_fullLib = merge(diff_enhancers_edger_fullLib, diffexpr, by="gene", sort=F)
  diff_enhancers_edger_fullLib = merge(diff_enhancers_edger_fullLib, simple_expr_change, by="gene", sort=F)
  # Simplistic expression foldchange
  smoothScatter(x=diff_enhancers_edger_fullLib$Fold, y=diff_enhancers_edger_fullLib$expr_fold_simple, xlab="Enhancer log2FC", ylab="Expression log2FC (Simplistic)", main="Enh: EdgeR lib=Full. Exp: Simple")
  corcoef = cor.test(diff_enhancers_edger_fullLib$Fold, diff_enhancers_edger_fullLib$expr_fold_simple, method="spearman")
  abline(lm(diff_enhancers_edger_fullLib$expr_fold_simple ~ diff_enhancers_edger_fullLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))
  # Deseq2 expression foldchange
  smoothScatter(x=diff_enhancers_edger_fullLib$Fold, y=diff_enhancers_edger_fullLib$expr_fold, xlab="Enhancer log2FC", ylab="Expression log2FC (Deseq2)", main="Enh: EdgeR lib=Full. Exp: Deseq2")
  corcoef = cor.test(diff_enhancers_edger_fullLib$Fold, diff_enhancers_edger_fullLib$expr_fold, method="spearman")
  abline(lm(diff_enhancers_edger_fullLib$expr_fold ~ diff_enhancers_edger_fullLib$Fold))
  abline(h=0,v=0, col="blue")
  text(x=grconvertX(0.5, from="npc"), y=grconvertY(.97, from="npc"), labels=paste0("Spearman rho=",signif(corcoef$estimate,digits=4), ", P=",signif(corcoef$p.value, digits=4)))



  ######################### Compare the results for the 4 different methods
  # Deseq2, effective library size
  if (!is.null(filterEnh)) {
    diff_enhancers_deseq_effLib = data.frame(dba.report(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, th=1), stringsAsFactors = F)
    diff_enhancers_deseq_effLib$enhancer = paste0(diff_enhancers_deseq_effLib$seqnames, ":", diff_enhancers_deseq_effLib$start, "-", diff_enhancers_deseq_effLib$end)
    diff_enhancers_deseq_effLib = diff_enhancers_deseq_effLib[diff_enhancers_deseq_effLib$enhancer %in% filterEnh,]
    diff_enhancers_deseq_effLib$FDR = p.adjust(diff_enhancers_deseq_effLib$p.value, method="BH")
    diff_enhancers_deseq_effLib = diff_enhancers_deseq_effLib[diff_enhancers_deseq_effLib$FDR<.05,]
  }
  else {
    diff_enhancers_deseq_effLib = data.frame(dba.report(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, th=0.05), stringsAsFactors = F)
    diff_enhancers_deseq_effLib$enhancer = paste0(diff_enhancers_deseq_effLib$seqnames, ":", diff_enhancers_deseq_effLib$start, "-", diff_enhancers_deseq_effLib$end)
  }

  # Deseq2, full library size (this is DiffBind's default method)
  if (!is.null(filterEnh)) {
    diff_enhancers_deseq_fullLib = data.frame(dba.report(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, th=1), stringsAsFactors = F)
    diff_enhancers_deseq_fullLib$enhancer = paste0(diff_enhancers_deseq_fullLib$seqnames, ":", diff_enhancers_deseq_fullLib$start, "-", diff_enhancers_deseq_fullLib$end)
    diff_enhancers_deseq_fullLib = diff_enhancers_deseq_fullLib[diff_enhancers_deseq_fullLib$enhancer %in% filterEnh,]
    diff_enhancers_deseq_fullLib$FDR = p.adjust(diff_enhancers_deseq_fullLib$p.value, method="BH")
    diff_enhancers_deseq_fullLib = diff_enhancers_deseq_fullLib[diff_enhancers_deseq_fullLib$FDR<.05,]
  }
  else {
    diff_enhancers_deseq_fullLib = data.frame(dba.report(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, th=0.05), stringsAsFactors = F)
    diff_enhancers_deseq_fullLib$enhancer = paste0(diff_enhancers_deseq_fullLib$seqnames, ":", diff_enhancers_deseq_fullLib$start, "-", diff_enhancers_deseq_fullLib$end)
  }

  # EdgeR, effective library size
  if (!is.null(filterEnh)) {
    diff_enhancers_edger_effLib = data.frame(dba.report(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, th=1), stringsAsFactors = F)
    diff_enhancers_edger_effLib$enhancer = paste0(diff_enhancers_edger_effLib$seqnames, ":", diff_enhancers_edger_effLib$start, "-", diff_enhancers_edger_effLib$end)
    diff_enhancers_edger_effLib = diff_enhancers_edger_effLib[diff_enhancers_edger_effLib$enhancer %in% filterEnh,]
    diff_enhancers_edger_effLib$FDR = p.adjust(diff_enhancers_edger_effLib$p.value, method="BH")
    diff_enhancers_edger_effLib = diff_enhancers_edger_effLib[diff_enhancers_edger_effLib$FDR<.05,]
  }
  else {
    diff_enhancers_edger_effLib = data.frame(dba.report(my_H3K27ac_dba_effLib, contrast=contrast, method=DBA_DESEQ2, th=0.05), stringsAsFactors = F)
    diff_enhancers_edger_effLib$enhancer = paste0(diff_enhancers_edger_effLib$seqnames, ":", diff_enhancers_edger_effLib$start, "-", diff_enhancers_edger_effLib$end)
  }

  # EdgeR, full library size
  if (!is.null(filterEnh)) {
    diff_enhancers_edger_fullLib = data.frame(dba.report(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, th=1), stringsAsFactors = F)
    diff_enhancers_edger_fullLib$enhancer = paste0(diff_enhancers_edger_fullLib$seqnames, ":", diff_enhancers_edger_fullLib$start, "-", diff_enhancers_edger_fullLib$end)
    diff_enhancers_edger_fullLib = diff_enhancers_edger_fullLib[diff_enhancers_edger_fullLib$enhancer %in% filterEnh,]
    diff_enhancers_edger_fullLib$FDR = p.adjust(diff_enhancers_edger_fullLib$p.value, method="BH")
    diff_enhancers_edger_fullLib = diff_enhancers_edger_fullLib[diff_enhancers_edger_fullLib$FDR<.05,]
  }
  else {
    diff_enhancers_edger_fullLib = data.frame(dba.report(my_H3K27ac_dba_fullLib, contrast=contrast, method=DBA_DESEQ2, th=0.05), stringsAsFactors = F)
    diff_enhancers_edger_fullLib$enhancer = paste0(diff_enhancers_edger_fullLib$seqnames, ":", diff_enhancers_edger_fullLib$start, "-", diff_enhancers_edger_fullLib$end)
  }


  ### plot venn diagram
  layout(mat=matrix(c(1,2), 2, 1))
  venn(data=list("Deseq2 effLib" = diff_enhancers_deseq_effLib[diff_enhancers_deseq_effLib$Fold>0, "enhancer"],
                 "Deseq2 fullLib" = diff_enhancers_deseq_fullLib[diff_enhancers_deseq_fullLib$Fold>0, "enhancer"],
                 "EdgeR effLib" = diff_enhancers_edger_effLib[diff_enhancers_edger_effLib$Fold>0, "enhancer"],
                 "EdgeR fullLib" = diff_enhancers_edger_fullLib[diff_enhancers_edger_fullLib$Fold>0, "enhancer"]))
  venn(data=list("Deseq2 effLib" = diff_enhancers_deseq_effLib[diff_enhancers_deseq_effLib$Fold<0, "enhancer"],
                 "Deseq2 fullLib" = diff_enhancers_deseq_fullLib[diff_enhancers_deseq_fullLib$Fold<0, "enhancer"],
                 "EdgeR effLib" = diff_enhancers_edger_effLib[diff_enhancers_edger_effLib$Fold<0, "enhancer"],
                 "EdgeR fullLib" = diff_enhancers_edger_fullLib[diff_enhancers_edger_fullLib$Fold<0, "enhancer"]))

}




analyze_geneset_enrichment = function(genes_interest, genes_univ, pathways_list) {
  ### fisher enrichment on genes with enhancer gain
  zz = lapply(pathways_list, function(x){
    genes_univ = unique(na.omit(genes_univ))
    pathway_genes = intersect(x, genes_univ)
    non_pathway_genes = genes_univ[!genes_univ %in% pathway_genes]
    genes_interest = intersect(genes_interest, genes_univ)
    non_genes_interest = genes_univ[!genes_univ %in% genes_interest]
    fish = fisher.test(matrix(c(length(intersect(pathway_genes, genes_interest)), length(intersect(non_pathway_genes, genes_interest)),
                                length(intersect(pathway_genes, non_genes_interest)), length(intersect(non_pathway_genes, non_genes_interest))), ncol=2), alternative="greater")
    data.frame(est=fish$estimate, p=fish$p.value, num_int=length(intersect(pathway_genes, genes_interest)), num_gs=length(pathway_genes), num_interest=length(genes_interest),
               int_genes = paste0(sort(intersect(pathway_genes, genes_interest)), collapse=","))
  })
  zz = do.call(rbind,zz)
  zz$padj = p.adjust(zz$p, method="BH")
  zz$pathway = rownames(zz)
  zz = zz[, c("pathway", "p", "padj", "est", "num_int", "num_gs", "num_interest", "int_genes")]
  zz = zz[order(zz$est<1, zz$p),]
  rownames(zz) = NULL
  return(zz)
}





get_gsea_genelist_filterExpr = function(K27ac_dba, contrast, method, corr_data, diffexprlist, filterEnh=NULL, filterExpr) {

  # keep only unique enhancer-gene pairs in target map
  corr_data = corr_data[!is.na(corr_data$gene),c("enhancer", "gene")]
  corr_data = corr_data[!duplicated(corr_data),]
  # keep only non-NA diffexpr foldchange
  if (filterExpr) {
    diffexpr = diffexprlist[[contrast]]
    diffexpr = diffexpr$diffexpr
    diffexpr = data.frame(gene=rownames(diffexpr), expr_fold=diffexpr$log2FoldChange, expr_stat=diffexpr$stat, expr_p=diffexpr$pvalue, expr_padj=diffexpr$padj, stringsAsFactors = F)
    diffexpr = diffexpr[!is.na(diffexpr$expr_fold),]
  }

  # Get full enhancer diffbind
  full_diffbind = data.frame(dba.report(K27ac_dba, contrast=contrast, method=method, th=1), stringsAsFactors = F)
  full_diffbind$enhancer = paste0(full_diffbind$seqnames, ":", full_diffbind$start, "-", full_diffbind$end)
  # if filterEnh not null, then filter enhancers and recalculate FDRs for enhancers
  if (!is.null(filterEnh)) {
    full_diffbind = full_diffbind[full_diffbind$enhancer %in% filterEnh,]
    full_diffbind$FDR = p.adjust(full_diffbind$p.value, method="BH")
  }
  full_diffbind = data.frame(enhancer=full_diffbind$enhancer, grp1_conc=full_diffbind[,7], grp2_conc=full_diffbind[,8], Fold=full_diffbind$Fold, FDR=full_diffbind$FDR, stringsAsFactors = F)
  full_diffbind = full_diffbind[full_diffbind$enhancer %in% corr_data$enhancer, ]
  full_diffbind = merge(full_diffbind, corr_data, by="enhancer", sort=F)
  if (filterExpr) {
    full_diffbind = merge(full_diffbind, diffexpr, by="gene", all.x=T, sort=F)
  }


  # Without filtering by gene expression
  if (!filterExpr) {
    full_diffbind_splitgene = split(full_diffbind, full_diffbind$gene)
    full_diffbind_collapsegene = sapply(full_diffbind_splitgene, function(x){
      log2(sum(2**as.numeric(x$grp1_conc))) - log2(sum(2**as.numeric(x$grp2_conc)))
    })
    return (full_diffbind_collapsegene)
  }

  # Filtered / validated by gene expression:
  if (filterExpr) {
    # # Medium strict: If enhancer FDR < .05, then expression should be significant in same direction (otherwise set group conc to 0). But if enhancer FDR >=.05, then expression should just be in same direction.
    # # Don't do this, because we don't want the stat to be dependent on P values which is affected by sample size!
    # full_diffbind_filtExp = full_diffbind[!is.na(full_diffbind$expr_fold) & !is.na(full_diffbind$expr_p),]
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp2_conc"] = 0
    # # More strict: If expr P >.05 or Fold is opposite from enhancer, then set enhancer group concentrations to 0
    # # Don't do this, because we don't want the stat to be dependent on P values which is affected by sample size!
    # full_diffbind_filtExp = full_diffbind[!is.na(full_diffbind$expr_fold) & !is.na(full_diffbind$expr_p),]
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp2_conc"] = 0
    # More lenient: If enhancer Fold and expr Fold are opposite, then set enhancer group concentrations to 0
    full_diffbind_filtExp = full_diffbind[!is.na(full_diffbind$expr_fold),]
    full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp1_conc"] = 0
    full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp2_conc"] = 0
    full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp1_conc"] = 0
    full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp2_conc"] = 0
    full_diffbind_filtExp_splitgene = split(full_diffbind_filtExp, full_diffbind_filtExp$gene)
    full_diffbind_filtExp_collapsegene = sapply(full_diffbind_filtExp_splitgene, function(x){
      log2(sum(2**as.numeric(x$grp1_conc))) - log2(sum(2**as.numeric(x$grp2_conc)))
    })
    return(full_diffbind_filtExp_collapsegene)
  }
}







# filterCell: 0 = don't filter by cell line, 1 = filter by cell line enhancer existence, 2 = filter by existence & same direction
get_gsea_genelist_multicomp_filterExpr_Cell = function(K27ac_dba1, K27ac_dba2, contrasts, group_interest, method, diffexprlist1, diffexprlist2, corr_data, filterEnh=NULL, filterExpr,
                                                       filterCell, K27ac_dba.cell, contrast.cell, group_interest.cell, enh_present_thr) {
  print(paste0("num contrasts: ", length(contrasts)))

  # keep only unique enhancer-gene pairs in target map
  corr_data = corr_data[!is.na(corr_data$gene),c("enhancer", "gene")]
  corr_data = corr_data[!duplicated(corr_data),]

  ### Get expression changes
  # keep only non-NA diffexpr foldchange
  if (filterExpr) {
    # create list of contrasts' DFs
    diffexprs_df = list()
    for (cc in 1:length(contrasts)) {
      zz_diffexpr = diffexprlist1[[contrasts[cc]]]
      zz_diffexpr = zz_diffexpr$diffexpr
      zz_diffexpr = data.frame(gene=rownames(zz_diffexpr), expr_fold=zz_diffexpr$log2FoldChange, expr_stat=zz_diffexpr$stat, expr_p=zz_diffexpr$pvalue, expr_padj=zz_diffexpr$padj, stringsAsFactors = F)
      zz_diffexpr = zz_diffexpr[!is.na(zz_diffexpr$expr_fold),]
      if (group_interest[1]==2) {
        zz_diffexpr$expr_fold1 = - zz_diffexpr$expr_fold1
        zz_diffexpr$expr_stat1 = - zz_diffexpr$expr_stat1
      }
      diffexprs_df[[cc]] = zz_diffexpr
    }
    # merge the contrats' DFs into single DF
    diffexpr_merge = diffexprs_df[[1]]
    if (length(diffexprs_df) > 1) {
      for (ii in 2:length(diffexprs_df)) {
        diffexpr_merge = merge(diffexpr_merge, diffexprs_df[[ii]], by="gene", suffixes=c("", paste0(".",ii)))
      }
    }
    colnames(diffexpr_merge)[2:5] = c("expr_fold.1", "expr_stat.1", "expr_p.1", "expr_padj.1")
    # If any of fold changes are in opposite directions, then set to 0. Otherwise take the "weaker" values
    fold_idx = grep('^expr_fold', colnames(diffexpr_merge))
    zzz = t(apply(diffexpr_merge, 1, function(x){
      opposite_fc = FALSE
      for (ii in 2:length(fold_idx)) {
        if (as.numeric(x[fold_idx[1]]) * as.numeric(x[fold_idx[ii]]) <= 0) {
          opposite_fc = TRUE
        }
      }
      if (opposite_fc==TRUE) {
        c(0, 0, 1, 1)
      } else {
        # get contrast with min FC
        minfc_name = names(x[fold_idx])[which.min(abs(as.numeric(x[fold_idx])))]
        minfc_idx = substring(minfc_name, 11)
        c(as.numeric(x[paste0("expr_fold.", minfc_idx)]), as.numeric(x[paste0("expr_stat.", minfc_idx)]), as.numeric(x[paste0("expr_p.", minfc_idx)]), as.numeric(x[paste0("expr_padj.", minfc_idx)]))
      }
    }))
    zzz = data.frame(zzz)
    colnames(zzz) = c("expr_fold", "expr_stat", "expr_p", "expr_padj")
    diffexpr_merge = cbind(diffexpr_merge, zzz)
    diffexpr = diffexpr_merge[, c("gene", "expr_fold", "expr_stat", "expr_p", "expr_padj")]
  }

  #### Get full enhancer diffbind
  # create list of contrasts' DFs
  full_diffbinds_df = list()
  for (cc in 1:length(contrasts)) {
    zz_full_diffbind = data.frame(dba.report(K27ac_dba1, contrast=contrasts[cc], method=method, th=1, bFlip=group_interest[1]==2), stringsAsFactors = F)
    zz_full_diffbind$enhancer = paste0(zz_full_diffbind$seqnames, ":", zz_full_diffbind$start, "-", zz_full_diffbind$end)
    # if filterEnh not null, then filter enhancers and recalculate FDRs for enhancers
    if (!is.null(filterEnh)) {
      zz_full_diffbind = zz_full_diffbind[zz_full_diffbind$enhancer %in% filterEnh,]
      zz_full_diffbind$FDR = p.adjust(zz_full_diffbind$p.value, method="BH")
    }
    zz_full_diffbind = data.frame(enhancer=zz_full_diffbind$enhancer, grp1_conc=zz_full_diffbind[,7], grp2_conc=zz_full_diffbind[,8], Fold=zz_full_diffbind$Fold, FDR=zz_full_diffbind$FDR, stringsAsFactors = F)
    full_diffbinds_df[[cc]] = zz_full_diffbind
  }
  # merge the contrats' DFs into single DF
  full_diffbind_merge = full_diffbinds_df[[1]]
  if (length(full_diffbinds_df) > 1) {
    for (ii in 2:length(full_diffbinds_df)) {
      full_diffbind_merge = merge(full_diffbind_merge, full_diffbinds_df[[ii]], by="enhancer", suffixes=c("", paste0(".",ii)))
    }
  }
  colnames(full_diffbind_merge)[2:5] = c("grp1_conc.1", "grp2_conc.1", "Fold.1", "FDR.1")

  # If any of fold changes are in opposite directions, then set to 0. Otherwise take the "weaker" values
  fold_idx = grep('^Fold.', colnames(full_diffbind_merge))
  zzz = t(apply(full_diffbind_merge, 1, function(x){
    opposite_fc = FALSE
    for (ii in 2:length(fold_idx)) {
      if (as.numeric(x[fold_idx[1]]) * as.numeric(x[fold_idx[ii]]) <= 0) {
        opposite_fc = TRUE
      }
      # print(paste0(ii, (as.numeric(x[fold_idx[1]]) * as.numeric(x[fold_idx[ii]]) <= 0), opposite_fc))
    }
    if (opposite_fc==TRUE) {
      c(0, 0, 0, 1)
    } else {
      # get contrast with min FC
      minfc_name = names(x[fold_idx])[which.min(abs(as.numeric(x[fold_idx])))]
      minfc_idx = substring(minfc_name, 6)
      c(as.numeric(x[paste0("grp1_conc.", minfc_idx)]), as.numeric(x[paste0("grp2_conc.", minfc_idx)]), as.numeric(x[paste0("Fold.", minfc_idx)]), as.numeric(x[paste0("FDR.", minfc_idx)]))
      # c(as.numeric(x[paste0("grp1_conc.", minfc_idx)]), as.numeric(x[paste0("grp2_conc.", minfc_idx)]), as.numeric(x["Fold1"]), max(as.numeric(x["FDR1"]),as.numeric(x["FDR2"])))
    }
  }))
  zzz = data.frame(zzz)
  colnames(zzz) = c("grp1_conc", "grp2_conc", "Fold", "FDR")
  full_diffbind_merge = cbind(full_diffbind_merge, zzz)
  full_diffbind = full_diffbind_merge[, c("enhancer", "grp1_conc", "grp2_conc", "Fold", "FDR")]

  ## add genes and expression
  full_diffbind = full_diffbind[full_diffbind$enhancer %in% corr_data$enhancer, ]
  full_diffbind = merge(full_diffbind, corr_data, by="enhancer", sort=F)
  if (filterExpr) {
    full_diffbind = merge(full_diffbind, diffexpr, by="gene", sort=F)
  }


  # filter enhancers by whether they exist in cell line.
  # By existence in cell: If an enhancer is gained (Fold>0), then it should exist in cell-line group1. Lost = should exist in cell-line group2. Otherwise, set grp_conc and fold to 0
  # By enhancer changes: If an enhancer is gained (Fold>0), then it should also be gained in cell-line. Otherwise, set grp_conc and fold to 0
  if (filterCell==1 | filterCell==2) {
    full_diffbind.cell = data.frame(dba.report(K27ac_dba.cell, contrast=contrast.cell, method=method, th=1, bFlip=group_interest.cell==2), stringsAsFactors = F)
    full_diffbind.cell$enhancer = paste0(full_diffbind.cell$seqnames, ":", full_diffbind.cell$start, "-", full_diffbind.cell$end)
    # if filterEnh not null, then filter enhancers and recalculate FDRs for enhancers
    if (!is.null(filterEnh)) {
      full_diffbind.cell = full_diffbind.cell[full_diffbind.cell$enhancer %in% filterEnh,]
      full_diffbind.cell$FDR = p.adjust(full_diffbind.cell$p.value, method="BH")
    }
    full_diffbind.cell$ConcPLW_1 = (2**full_diffbind.cell[,7]) / log2(full_diffbind.cell$width)
    full_diffbind.cell$ConcPLW_2 = (2**full_diffbind.cell[,8]) / log2(full_diffbind.cell$width)
    ConcPLW_dist = c(full_diffbind.cell$ConcPLW_1, full_diffbind.cell$ConcPLW_2)
    full_diffbind.cell = full_diffbind.cell[, c("enhancer", "Fold", "p.value", "FDR", "ConcPLW_1", "ConcPLW_2")]
    colnames(full_diffbind.cell) = c("enhancer", "cell_Fold", "cell_p.value", "cell_FDR", "cell_ConcPLW_1", "cell_ConcPLW_2")
    full_diffbind_filtCell = merge(full_diffbind, full_diffbind.cell, by="enhancer", sort=F)
    ConcPLW_dist_thr = quantile(ConcPLW_dist, enh_present_thr)

    # filter by existence in cell-line
    full_diffbind_filtCell[full_diffbind_filtCell$Fold>0 & full_diffbind_filtCell$grp1_conc <= ConcPLW_dist_thr, "grp1_conc"] = 0
    full_diffbind_filtCell[full_diffbind_filtCell$Fold>0 & full_diffbind_filtCell$grp1_conc <= ConcPLW_dist_thr, "grp2_conc"] = 0
    full_diffbind_filtCell[full_diffbind_filtCell$Fold>0 & full_diffbind_filtCell$grp1_conc <= ConcPLW_dist_thr, "Fold"] = 0
    full_diffbind_filtCell[full_diffbind_filtCell$Fold<0 & full_diffbind_filtCell$grp2_conc <= ConcPLW_dist_thr, "grp1_conc"] = 0
    full_diffbind_filtCell[full_diffbind_filtCell$Fold<0 & full_diffbind_filtCell$grp2_conc <= ConcPLW_dist_thr, "grp2_conc"] = 0
    full_diffbind_filtCell[full_diffbind_filtCell$Fold<0 & full_diffbind_filtCell$grp2_conc <= ConcPLW_dist_thr, "Fold"] = 0

    # filter by enhancer change in cell-line
    if (filterCell==2) {
      full_diffbind_filtCell[full_diffbind_filtCell$Fold>0 & full_diffbind_filtCell$cell_Fold < 0, "grp1_conc"] = 0
      full_diffbind_filtCell[full_diffbind_filtCell$Fold>0 & full_diffbind_filtCell$cell_Fold < 0, "grp2_conc"] = 0
      full_diffbind_filtCell[full_diffbind_filtCell$Fold>0 & full_diffbind_filtCell$cell_Fold < 0, "Fold"] = 0
      full_diffbind_filtCell[full_diffbind_filtCell$Fold<0 & full_diffbind_filtCell$cell_Fold > 0, "grp1_conc"] = 0
      full_diffbind_filtCell[full_diffbind_filtCell$Fold<0 & full_diffbind_filtCell$cell_Fold > 0, "grp2_conc"] = 0
      full_diffbind_filtCell[full_diffbind_filtCell$Fold<0 & full_diffbind_filtCell$cell_Fold > 0, "Fold"] = 0
    }

    full_diffbind = full_diffbind_filtCell
  }


  # Without filtering by gene expression
  if (!filterExpr) {
    full_diffbind_splitgene = split(full_diffbind, full_diffbind$gene)
    full_diffbind_collapsegene = sapply(full_diffbind_splitgene, function(x){
      log2(sum(2**as.numeric(x$grp1_conc))) - log2(sum(2**as.numeric(x$grp2_conc)))
    })
    return (full_diffbind_collapsegene)
  }

  # Filtered / validated by gene expression:
  if (filterExpr) {
    # # Medium strict: If enhancer FDR < .05, then expression should be significant in same direction (otherwise set group conc to 0). But if enhancer FDR >=.05, then expression should just be in same direction.
    # # Don't do this, because we don't want the stat to be dependent on P values which is affected by sample size!
    # full_diffbind_filtExp = full_diffbind[!is.na(full_diffbind$expr_fold) & !is.na(full_diffbind$expr_p),]
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR<.05 & full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$FDR>=.05 & full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp2_conc"] = 0
    # # More strict: If expr P >.05 or Fold is opposite from enhancer, then set enhancer group concentrations to 0
    # # Don't do this, because we don't want the stat to be dependent on P values which is affected by sample size!
    # full_diffbind_filtExp = full_diffbind[!is.na(full_diffbind$expr_fold) & !is.na(full_diffbind$expr_p),]
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold < 0), "grp2_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp1_conc"] = 0
    # full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & (full_diffbind_filtExp$expr_p>=.05 | full_diffbind_filtExp$expr_fold > 0), "grp2_conc"] = 0
    # More lenient: If enhancer Fold and expr Fold are opposite, then set enhancer group concentrations to 0
    full_diffbind_filtExp = full_diffbind[!is.na(full_diffbind$expr_fold),]
    full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp1_conc"] = 0
    full_diffbind_filtExp[full_diffbind_filtExp$Fold > 0 & full_diffbind_filtExp$expr_fold < 0, "grp2_conc"] = 0
    full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp1_conc"] = 0
    full_diffbind_filtExp[full_diffbind_filtExp$Fold < 0 & full_diffbind_filtExp$expr_fold > 0, "grp2_conc"] = 0
    full_diffbind_filtExp_splitgene = split(full_diffbind_filtExp, full_diffbind_filtExp$gene)
    full_diffbind_filtExp_collapsegene = sapply(full_diffbind_filtExp_splitgene, function(x){
      log2(sum(2**as.numeric(x$grp1_conc))) - log2(sum(2**as.numeric(x$grp2_conc)))
    })
    return(full_diffbind_filtExp_collapsegene)
  }
}


















########### Validate differential enhancers from multiple comparisons, with a validation set
# Eg. Fluke-specific enhancers: Gain = gain in Fluke vs N, gain in Fluke vs Non-Fluke. Lost: likewise
# Inputs: constrasts.x = vector of contrasts. group_interest.x = indices of the group of interest in each of the contrasts. Likewise for contrasts.y & group_interest.y.
#
# Strategy:
#  1. Get gain in first comparison, then among these check gain in second comparison (shrunken hypothesis space)
#  2. Tissue: expression change in same direction for both comparisons (FDR<.05)
#  3. Cell line: enhancer present in the higher condition for both comparisons  (2**conc / log2(enh_width) is not in bottom enh_present_thr)
#  4. Cell line: enhancer change in same direction for both comparisons (Fold>.2)
#  5. Cell line: expression change in same direction for both comparisons (Fold>.1)
# filterEnh: only consider enhancers in this list. This is used when analyzing superenhancers diffbind, where I gave it enhancers+superenhancers so that normalization factors won't be based on superenhancers.
# Note that when merging enhancers with target map, enhancers not in corr_data will be removed anyway.
# When filterEnh is not NULL, adjusted P-values will be recalculated.
# test: validate_diff_enhancers_multicomp_exp_cellline(K27ac_dba.x.1=K27ac_dba_deseq_contrasts1, K27ac_dba.x.2=K27ac_dba_deseq_contrasts1, contrasts.x=c(1, 5), group_interest.x=c(1,1), method.x=DBA_DESEQ2, diffexprlist.x.1=diffexpr_list_contrasts1, diffexprlist.x.2=diffexpr_list_contrasts1,
#                                                      K27ac_dba.y.1=K27ac_cells_tisspeaks_dba, K27ac_dba.y.2=K27ac_cells_tisspeaks_dba, contrasts.y=c(2,4), group_interest.y=c(1,1), method.y=DBA_DESEQ2, diffexprlist.y.1=diffexpr_cells_list, diffexprlist.y.2=diffexpr_cells_list,
#                                                      corr_data=target_map, enh_present_thr=.5, filterEnh=NULL, pathways_list=msigdb_H_list, exp_p_or_q="p")
validate_diff_enhancers_multicomp_exp_cellline = function(K27ac_dba.x.1, K27ac_dba.x.2, contrasts.x, group_interest.x, method.x, diffexprlist.x.1, diffexprlist.x.2,
                                                          K27ac_dba.y.1, K27ac_dba.y.2, contrasts.y, group_interest.y, method.y, diffexprlist.y.1, diffexprlist.y.2,
                                                          corr_data, enh_present_thr=.5, filterEnh=NULL, pathways_list, enh_p_or_q="q", exp_p_or_q="p") {


  # keep only unique enhancer-gene pairs in target map
  corr_data = corr_data[!is.na(corr_data$gene),c("enhancer", "gene")]
  corr_data = corr_data[!duplicated(corr_data),]

  ### get all diff expr, keep only non-NA diffexpr foldchange
  # tissue enhancers, contrast 1
  diffexpr.x.1 = diffexprlist.x.1[[contrasts.x[1]]]
  diffexpr.x.1 = diffexpr.x.1$diffexpr
  diffexpr.x.1 = data.frame(gene=rownames(diffexpr.x.1), expr_fold.x.1=diffexpr.x.1$log2FoldChange, expr_stat.x.1=diffexpr.x.1$stat, expr_p.x.1=diffexpr.x.1$pvalue, expr_padj.x.1=diffexpr.x.1$padj, stringsAsFactors = F)
  diffexpr.x.1 = diffexpr.x.1[!is.na(diffexpr.x.1$expr_fold.x.1),]
  if (group_interest.x[1]==2) {
    diffexpr.x.1$expr_fold.x.1 = - diffexpr.x.1$expr_fold.x.1
    diffexpr.x.1$expr_stat.x.1 = - diffexpr.x.1$expr_stat.x.1
  }
  # tissue enhancers, contrast 2
  diffexpr.x.2 = diffexprlist.x.2[[contrasts.x[2]]]
  diffexpr.x.2 = diffexpr.x.2$diffexpr
  diffexpr.x.2 = data.frame(gene=rownames(diffexpr.x.2), expr_fold.x.2=diffexpr.x.2$log2FoldChange, expr_stat.x.2=diffexpr.x.2$stat, expr_p.x.2=diffexpr.x.2$pvalue, expr_padj.x.2=diffexpr.x.2$padj, stringsAsFactors = F)
  diffexpr.x.2 = diffexpr.x.2[!is.na(diffexpr.x.2$expr_fold.x.2),]
  if (group_interest.x[2]==2) {
    diffexpr.x.2$expr_fold.x.2 = - diffexpr.x.2$expr_fold.x.2
    diffexpr.x.2$expr_stat.x.2 = - diffexpr.x.2$expr_stat.x.2
  }
  # cell-line enhancers, contrast 1
  diffexpr.y.1 = diffexprlist.y.1[[contrasts.y[1]]]
  diffexpr.y.1 = diffexpr.y.1$diffexpr
  diffexpr.y.1 = data.frame(gene=rownames(diffexpr.y.1), expr_fold.y.1=diffexpr.y.1$log2FoldChange, expr_stat.y.1=diffexpr.y.1$stat, expr_p.y.1=diffexpr.y.1$pvalue, expr_padj.y.1=diffexpr.y.1$padj, stringsAsFactors = F)
  diffexpr.y.1 = diffexpr.y.1[!is.na(diffexpr.y.1$expr_fold.y.1),]
  if (group_interest.y[1]==2) {
    diffexpr.y.1$expr_fold.y.1 = - diffexpr.y.1$expr_fold.y.1
    diffexpr.y.1$expr_stat.y.1 = - diffexpr.y.1$expr_stat.y.1
  }
  # cell-line enhancers, contrast 2
  diffexpr.y.2 = diffexprlist.y.2[[contrasts.y[2]]]
  diffexpr.y.2 = diffexpr.y.2$diffexpr
  diffexpr.y.2 = data.frame(gene=rownames(diffexpr.y.2), expr_fold.y.2=diffexpr.y.2$log2FoldChange, expr_stat.y.2=diffexpr.y.2$stat, expr_p.y.2=diffexpr.y.2$pvalue, expr_padj.y.2=diffexpr.y.2$padj, stringsAsFactors = F)
  diffexpr.y.2 = diffexpr.y.2[!is.na(diffexpr.y.2$expr_fold.y.2),]
  if (group_interest.y[2]==2) {
    diffexpr.y.2$expr_fold.y.2 = - diffexpr.y.2$expr_fold.y.2
    diffexpr.y.2$expr_stat.y.2 = - diffexpr.y.2$expr_stat.y.2
  }

  ### Create table of ALL enhancers, with enhancer and expression FC in both contrasts, in tissue and cell line. No need to recalculate FDRs (except for when filterEnh not null
  # tissue enhancers, contrast 1
  diff_enhancers.x.1 = data.frame(dba.report(K27ac_dba.x.1, contrast=contrasts.x[1], method=method.x, th=1, bFlip=group_interest.x[1]==2), stringsAsFactors = F)
  diff_enhancers.x.1$enhancer = paste0(diff_enhancers.x.1$seqnames, ":", diff_enhancers.x.1$start, "-", diff_enhancers.x.1$end)
  diff_enhancers.x.1$ConcPLW_1 = (2**diff_enhancers.x.1[,7]) / log2(diff_enhancers.x.1$width)
  diff_enhancers.x.1$ConcPLW_2 = (2**diff_enhancers.x.1[,8]) / log2(diff_enhancers.x.1$width)
  diff_enhancers.x.1 = diff_enhancers.x.1[, c("enhancer", "width", "ConcPLW_1", "ConcPLW_2", "Fold", "p.value", "FDR")]
  colnames(diff_enhancers.x.1) = c("enhancer", "width", "ConcPLW_1.x.1", "ConcPLW_2.x.1", "Fold.x.1", "p.value.x.1", "FDR.x.1")
  # tissue enhancers, contrast 2
  diff_enhancers.x.2 = data.frame(dba.report(K27ac_dba.x.2, contrast=contrasts.x[2], method=method.x, th=1, bFlip=group_interest.x[2]==2), stringsAsFactors = F)
  diff_enhancers.x.2$enhancer = paste0(diff_enhancers.x.2$seqnames, ":", diff_enhancers.x.2$start, "-", diff_enhancers.x.2$end)
  diff_enhancers.x.2$ConcPLW_1 = (2**diff_enhancers.x.2[,7]) / log2(diff_enhancers.x.2$width)
  diff_enhancers.x.2$ConcPLW_2 = (2**diff_enhancers.x.2[,8]) / log2(diff_enhancers.x.2$width)
  diff_enhancers.x.2 = diff_enhancers.x.2[, c("enhancer", "ConcPLW_1", "ConcPLW_2", "Fold", "p.value", "FDR")]
  colnames(diff_enhancers.x.2) = c("enhancer", "ConcPLW_1.x.2", "ConcPLW_2.x.2", "Fold.x.2", "p.value.x.2", "FDR.x.2")
  # cell-line enhancers, contrast 1
  diff_enhancers.y.1 = data.frame(dba.report(K27ac_dba.y.1, contrast=contrasts.y[1], method=method.y, th=1, bFlip=group_interest.y[1]==2), stringsAsFactors = F)
  diff_enhancers.y.1$enhancer = paste0(diff_enhancers.y.1$seqnames, ":", diff_enhancers.y.1$start, "-", diff_enhancers.y.1$end)
  diff_enhancers.y.1$ConcPLW_1 = (2**diff_enhancers.y.1[,7]) / log2(diff_enhancers.y.1$width)
  diff_enhancers.y.1$ConcPLW_2 = (2**diff_enhancers.y.1[,8]) / log2(diff_enhancers.y.1$width)
  diff_enhancers.y.1 = diff_enhancers.y.1[, c("enhancer", "ConcPLW_1", "ConcPLW_2", "Fold", "p.value", "FDR")]
  colnames(diff_enhancers.y.1) = c("enhancer", "ConcPLW_1.y.1", "ConcPLW_2.y.1", "Fold.y.1", "p.value.y.1", "FDR.y.1")
  # cell-line enhancers, contrast 2
  diff_enhancers.y.2 = data.frame(dba.report(K27ac_dba.y.2, contrast=contrasts.y[2], method=method.y, th=1, bFlip=group_interest.y[2]==2), stringsAsFactors = F)
  diff_enhancers.y.2$enhancer = paste0(diff_enhancers.y.2$seqnames, ":", diff_enhancers.y.2$start, "-", diff_enhancers.y.2$end)
  diff_enhancers.y.2$ConcPLW_1 = (2**diff_enhancers.y.2[,7]) / log2(diff_enhancers.y.2$width)
  diff_enhancers.y.2$ConcPLW_2 = (2**diff_enhancers.y.2[,8]) / log2(diff_enhancers.y.2$width)
  diff_enhancers.y.2 = diff_enhancers.y.2[, c("enhancer", "ConcPLW_1", "ConcPLW_2", "Fold", "p.value", "FDR")]
  colnames(diff_enhancers.y.2) = c("enhancer", "ConcPLW_1.y.2", "ConcPLW_2.y.2", "Fold.y.2", "p.value.y.2", "FDR.y.2")

  # combine
  diff_enhancers_merge = merge(diff_enhancers.x.1, diff_enhancers.x.2, by="enhancer", sort=F)
  diff_enhancers_merge = merge(diff_enhancers_merge, diff_enhancers.y.1, by="enhancer", sort=F)
  diff_enhancers_merge = merge(diff_enhancers_merge, diff_enhancers.y.2, by="enhancer", sort=F) # 474060
  # if filterEnh not null, then filter enhancers and recalculate FDRs for enhancers
  if (!is.null(filterEnh)) {
    diff_enhancers_merge = diff_enhancers_merge[diff_enhancers_merge$enhancer %in% filterEnh,]
    diff_enhancers_merge$FDR.x.1 = p.adjust(diff_enhancers_merge$p.value.x.1, method="BH")
    diff_enhancers_merge$FDR.x.2 = p.adjust(diff_enhancers_merge$p.value.x.2, method="BH")
    diff_enhancers_merge$FDR.y.1 = p.adjust(diff_enhancers_merge$p.value.y.1, method="BH")
    diff_enhancers_merge$FDR.y.2 = p.adjust(diff_enhancers_merge$p.value.y.2, method="BH")
  }
  # get distribution of cell line enhancer signals
  ConcPLW_dist_1.y.1 = diff_enhancers_merge$ConcPLW_1.y.1 # contrast 1, group 1 (use this for contrast 1 gain)
  ConcPLW_dist_2.y.1 = diff_enhancers_merge$ConcPLW_2.y.1 # contrast 1, group 2 (use this for contrast 1 loss)
  ConcPLW_dist_1.y.2 = diff_enhancers_merge$ConcPLW_1.y.2 # contrast 2, group 1 (use this for contrast 2 gain)
  ConcPLW_dist_2.y.2 = diff_enhancers_merge$ConcPLW_2.y.2 # contrast 2, group 2 (use this for contrast 2 loss)
  ConcPLW_dist.y.1 = c(diff_enhancers_merge$ConcPLW_1.y.1, diff_enhancers_merge$ConcPLW_2.y.1) # contrast 1
  ConcPLW_dist.y.2 = c(diff_enhancers_merge$ConcPLW_1.y.2, diff_enhancers_merge$ConcPLW_2.y.2) # contrast 2
  # add genes and expression
  diff_enhancers_mergegene = merge(diff_enhancers_merge, corr_data, by="enhancer", all.x=T, sort=F)
  diff_enhancers_mergegene = merge(diff_enhancers_mergegene, diffexpr.x.1, by="gene", all.x=T, sort=F)
  diff_enhancers_mergegene = merge(diff_enhancers_mergegene, diffexpr.x.2, by="gene", all.x=T, sort=F)
  diff_enhancers_mergegene = merge(diff_enhancers_mergegene, diffexpr.y.1, by="gene", all.x=T, sort=F)
  diff_enhancers_mergegene = merge(diff_enhancers_mergegene, diffexpr.y.2, by="gene", all.x=T, sort=F)


  #### Create table of only enhancers gained or lost in tissue contrast 1. Recalculate FDRs for cellline enhancers, and all expr
  #### Note that if I filter FDR<.05 on the full list, somehow its slightly different from doing dba.report(th=.05)! So call dba.report again
  ## Get enhancers gained/lost in tissue in contrast1
  if (!is.null(filterEnh)) {
    zz_diff_enhancers.x.1 = data.frame(dba.report(K27ac_dba.x.1, contrast=contrasts.x[1], method=method.x, th=1, bFlip=group_interest.x[1]==2), stringsAsFactors = F)
    zz_diff_enhancers.x.1$enhancer = paste0(zz_diff_enhancers.x.1$seqnames, ":", zz_diff_enhancers.x.1$start, "-", zz_diff_enhancers.x.1$end)
    zz_diff_enhancers.x.1 = zz_diff_enhancers.x.1[zz_diff_enhancers.x.1$enhancer %in% filterEnh,]
    zz_diff_enhancers.x.1$FDR = p.adjust(zz_diff_enhancers.x.1$p.value, method="BH")
    if (enh_p_or_q=="q")
      zz_diff_enhancers.x.1 = zz_diff_enhancers.x.1[zz_diff_enhancers.x.1$FDR<.05,]
    else if (enh_p_or_q=="p")
      zz_diff_enhancers.x.1 = zz_diff_enhancers.x.1[zz_diff_enhancers.x.1$p.value<.05,]
  }
  else {
    if (enh_p_or_q=="q")
      zz_diff_enhancers.x.1 = data.frame(dba.report(K27ac_dba.x.1, contrast=contrasts.x[1], method=method.x, th=0.05, bFlip=group_interest.x[1]==2), stringsAsFactors = F)
    else if (enh_p_or_q=="p")
      zz_diff_enhancers.x.1 = data.frame(dba.report(K27ac_dba.x.1, contrast=contrasts.x[1], method=method.x, bUsePval=TRUE, th=0.05, bFlip=group_interest.x[1]==2), stringsAsFactors = F)
    zz_diff_enhancers.x.1$enhancer = paste0(zz_diff_enhancers.x.1$seqnames, ":", zz_diff_enhancers.x.1$start, "-", zz_diff_enhancers.x.1$end)
  }
  zz_diff_enhancers.x.1_gain = unique(zz_diff_enhancers.x.1[zz_diff_enhancers.x.1$Fold>0, "enhancer"]) # 85300
  zz_diff_enhancers.x.1_lost = unique(zz_diff_enhancers.x.1[zz_diff_enhancers.x.1$Fold<0, "enhancer"]) # 36242
  # Extract the enhancers gained/lost in tissue from the full table calculated above
  diff_enhancers_merge_enhX1_gain = diff_enhancers_merge[diff_enhancers_merge$enhancer %in% zz_diff_enhancers.x.1_gain,] # 85300
  diff_enhancers_merge_enhX1_lost = diff_enhancers_merge[diff_enhancers_merge$enhancer %in% zz_diff_enhancers.x.1_lost,] # 36242
  # Recalculate FDRs for tissue enhancers contrast 2
  diff_enhancers_merge_enhX1_gain$FDR.x.2.shrunk = p.adjust(diff_enhancers_merge_enhX1_gain$p.value.x.2, method="BH") # 30034 --> 31642
  diff_enhancers_merge_enhX1_lost$FDR.x.2.shrunk = p.adjust(diff_enhancers_merge_enhX1_lost$p.value.x.2, method="BH") # 11412 --> 11632
  # Recalculate FDRs for cell-line enhancers contrasts 1 & 2
  diff_enhancers_merge_enhX1_gain$FDR.y.1.shrunk = p.adjust(diff_enhancers_merge_enhX1_gain$p.value.y.1, method="BH") # 488 --> 110 Using P: 5699
  diff_enhancers_merge_enhX1_gain$FDR.y.2.shrunk = p.adjust(diff_enhancers_merge_enhX1_gain$p.value.y.2, method="BH") # 536 -->  301 Using P: 4478
  diff_enhancers_merge_enhX1_lost$FDR.y.1.shrunk = p.adjust(diff_enhancers_merge_enhX1_lost$p.value.y.1, method="BH") # 457 --> 168 Using P: 3438
  diff_enhancers_merge_enhX1_lost$FDR.y.2.shrunk = p.adjust(diff_enhancers_merge_enhX1_lost$p.value.y.2, method="BH") #  42 --> 0 Using P: 2078
  # Reorder the columns so that recalculated FDRs are next to the original FDRs
  diff_enhancers_merge_enhX1_gain = diff_enhancers_merge_enhX1_gain[, c("enhancer", "width", "ConcPLW_1.x.1", "ConcPLW_2.x.1", "Fold.x.1", "p.value.x.1", "FDR.x.1",
                                                                        "ConcPLW_1.x.2", "ConcPLW_2.x.2", "Fold.x.2", "p.value.x.2", "FDR.x.2", "FDR.x.2.shrunk",
                                                                        "ConcPLW_1.y.1", "ConcPLW_2.y.1", "Fold.y.1", "p.value.y.1", "FDR.y.1", "FDR.y.1.shrunk",
                                                                        "ConcPLW_1.y.2", "ConcPLW_2.y.2", "Fold.y.2", "p.value.y.2", "FDR.y.2", "FDR.y.2.shrunk")]
  diff_enhancers_merge_enhX1_lost = diff_enhancers_merge_enhX1_lost[, c("enhancer", "width", "ConcPLW_1.x.1", "ConcPLW_2.x.1", "Fold.x.1", "p.value.x.1", "FDR.x.1",
                                                                        "ConcPLW_1.x.2", "ConcPLW_2.x.2", "Fold.x.2", "p.value.x.2", "FDR.x.2", "FDR.x.2.shrunk",
                                                                        "ConcPLW_1.y.1", "ConcPLW_2.y.1", "Fold.y.1", "p.value.y.1", "FDR.y.1", "FDR.y.1.shrunk",
                                                                        "ConcPLW_1.y.2", "ConcPLW_2.y.2", "Fold.y.2", "p.value.y.2", "FDR.y.2", "FDR.y.2.shrunk")]
  ## add genes and expression, recalculate FDRs
  # gained enhancers
  diff_enhancers_mergegene_enhX1_gain = merge(diff_enhancers_merge_enhX1_gain, corr_data, by="enhancer", all.x=T, sort=F) # 114194 rows, 8598 genes
  zz_diffexpr.x.1 = diffexpr.x.1[diffexpr.x.1$gene %in% na.omit(diff_enhancers_mergegene_enhX1_gain$gene),]
  zz_diffexpr.x.1$expr_padj.x.1.shrunk = p.adjust(zz_diffexpr.x.1$expr_p.x.1, method="BH")
  zz_diffexpr.x.2 = diffexpr.x.2[diffexpr.x.2$gene %in% na.omit(diff_enhancers_mergegene_enhX1_gain$gene),]
  zz_diffexpr.x.2$expr_padj.x.2.shrunk = p.adjust(zz_diffexpr.x.2$expr_p.x.2, method="BH")
  zz_diffexpr.y.1 = diffexpr.y.1[diffexpr.y.1$gene %in% na.omit(diff_enhancers_mergegene_enhX1_gain$gene),]
  zz_diffexpr.y.1$expr_padj.y.1.shrunk = p.adjust(zz_diffexpr.y.1$expr_p.y.1, method="BH")
  zz_diffexpr.y.2 = diffexpr.y.2[diffexpr.y.2$gene %in% na.omit(diff_enhancers_mergegene_enhX1_gain$gene),]
  zz_diffexpr.y.2$expr_padj.y.2.shrunk = p.adjust(zz_diffexpr.y.2$expr_p.y.2, method="BH")
  diff_enhancers_mergegene_enhX1_gain = merge(diff_enhancers_mergegene_enhX1_gain, zz_diffexpr.x.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 2105 --> 2058 Using P:
  diff_enhancers_mergegene_enhX1_gain = merge(diff_enhancers_mergegene_enhX1_gain, zz_diffexpr.x.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 1592 --> 1719 Using P:
  diff_enhancers_mergegene_enhX1_gain = merge(diff_enhancers_mergegene_enhX1_gain, zz_diffexpr.y.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 620 --> 605 Using P:
  diff_enhancers_mergegene_enhX1_gain = merge(diff_enhancers_mergegene_enhX1_gain, zz_diffexpr.y.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 226 --> 214 Using P:
  # lost enhancers
  diff_enhancers_mergegene_enhX1_lost = merge(diff_enhancers_merge_enhX1_lost, corr_data, by="enhancer", all.x=T, sort=F)
  zz_diffexpr.x.1 = diffexpr.x.1[diffexpr.x.1$gene %in% na.omit(diff_enhancers_mergegene_enhX1_lost$gene),]
  zz_diffexpr.x.1$expr_padj.x.1.shrunk = p.adjust(zz_diffexpr.x.1$expr_p.x.1, method="BH")
  zz_diffexpr.x.2 = diffexpr.x.2[diffexpr.x.2$gene %in% na.omit(diff_enhancers_mergegene_enhX1_lost$gene),]
  zz_diffexpr.x.2$expr_padj.x.2.shrunk = p.adjust(zz_diffexpr.x.2$expr_p.x.2, method="BH")
  zz_diffexpr.y.1 = diffexpr.y.1[diffexpr.y.1$gene %in% na.omit(diff_enhancers_mergegene_enhX1_lost$gene),]
  zz_diffexpr.y.1$expr_padj.y.1.shrunk = p.adjust(zz_diffexpr.y.1$expr_p.y.1, method="BH")
  zz_diffexpr.y.2 = diffexpr.y.2[diffexpr.y.2$gene %in% na.omit(diff_enhancers_mergegene_enhX1_lost$gene),]
  zz_diffexpr.y.2$expr_padj.y.2.shrunk = p.adjust(zz_diffexpr.y.2$expr_p.y.2, method="BH")
  diff_enhancers_mergegene_enhX1_lost = merge(diff_enhancers_mergegene_enhX1_lost, zz_diffexpr.x.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 1367 --> 1447 Using P:
  diff_enhancers_mergegene_enhX1_lost = merge(diff_enhancers_mergegene_enhX1_lost, zz_diffexpr.x.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes  742 --> 795 Using P:
  diff_enhancers_mergegene_enhX1_lost = merge(diff_enhancers_mergegene_enhX1_lost, zz_diffexpr.y.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 171 --> 175 Using P:
  diff_enhancers_mergegene_enhX1_lost = merge(diff_enhancers_mergegene_enhX1_lost, zz_diffexpr.y.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 56 --> 59 Using P:
  rm(zz_diffexpr.x.1, zz_diffexpr.x.2, zz_diffexpr.y.1, zz_diffexpr.y.2)


  ## Among enhancers gained/lost in tissue in contrast1, get those gained/lost in tissue in contrast2 (using shrunken hypothesis space)
  # Extract the enhancers gained/lost in tissue from the full table calculated above
  diff_enhancers_merge_enhX_gain = diff_enhancers_merge[diff_enhancers_merge$enhancer %in% zz_diff_enhancers.x.1_gain,] # 85300
  diff_enhancers_merge_enhX_lost = diff_enhancers_merge[diff_enhancers_merge$enhancer %in% zz_diff_enhancers.x.1_lost,] # 36242
  # Recalculate FDRs for tissue enhancers contrast 2, keep those FDR<.05
  diff_enhancers_merge_enhX_gain$FDR.x.2.shrunk = p.adjust(diff_enhancers_merge_enhX_gain$p.value.x.2, method="BH") # 30034 --> 31642
  diff_enhancers_merge_enhX_lost$FDR.x.2.shrunk = p.adjust(diff_enhancers_merge_enhX_lost$p.value.x.2, method="BH") # 11412 --> 11632
  if (enh_p_or_q=="q") {
    diff_enhancers_merge_enhX_gain = diff_enhancers_merge_enhX_gain[diff_enhancers_merge_enhX_gain$FDR.x.2.shrunk < .05 & diff_enhancers_merge_enhX_gain$Fold.x.2 > 0,] # 31642
    diff_enhancers_merge_enhX_lost = diff_enhancers_merge_enhX_lost[diff_enhancers_merge_enhX_lost$FDR.x.2.shrunk < .05 & diff_enhancers_merge_enhX_lost$Fold.x.2 < 0,] # 11632
  }
  else if (enh_p_or_q=="p") {
    diff_enhancers_merge_enhX_gain = diff_enhancers_merge_enhX_gain[diff_enhancers_merge_enhX_gain$p.value.x.2 < .05 & diff_enhancers_merge_enhX_gain$Fold.x.2 > 0,] # 31642
    diff_enhancers_merge_enhX_lost = diff_enhancers_merge_enhX_lost[diff_enhancers_merge_enhX_lost$p.value.x.2 < .05 & diff_enhancers_merge_enhX_lost$Fold.x.2 < 0,] # 11632
  }
  # Recalculate FDRs for cell-line enhancers contrasts 1 & 2
  diff_enhancers_merge_enhX_gain$FDR.y.1.shrunk = p.adjust(diff_enhancers_merge_enhX_gain$p.value.y.1, method="BH") # 206 --> 45
  diff_enhancers_merge_enhX_gain$FDR.y.2.shrunk = p.adjust(diff_enhancers_merge_enhX_gain$p.value.y.2, method="BH") # 153 --> 114
  diff_enhancers_merge_enhX_lost$FDR.y.1.shrunk = p.adjust(diff_enhancers_merge_enhX_lost$p.value.y.1, method="BH") # 186 --> 96
  diff_enhancers_merge_enhX_lost$FDR.y.2.shrunk = p.adjust(diff_enhancers_merge_enhX_lost$p.value.y.2, method="BH") # 10 --> 0
  # Reorder the columns so that recalculated FDRs are next to the original FDRs
  diff_enhancers_merge_enhX_gain = diff_enhancers_merge_enhX_gain[, c("enhancer", "width", "ConcPLW_1.x.1", "ConcPLW_2.x.1", "Fold.x.1", "p.value.x.1", "FDR.x.1",
                                                                      "ConcPLW_1.x.2", "ConcPLW_2.x.2", "Fold.x.2", "p.value.x.2", "FDR.x.2", "FDR.x.2.shrunk",
                                                                      "ConcPLW_1.y.1", "ConcPLW_2.y.1", "Fold.y.1", "p.value.y.1", "FDR.y.1", "FDR.y.1.shrunk",
                                                                      "ConcPLW_1.y.2", "ConcPLW_2.y.2", "Fold.y.2", "p.value.y.2", "FDR.y.2", "FDR.y.2.shrunk")]
  diff_enhancers_merge_enhX_lost = diff_enhancers_merge_enhX_lost[, c("enhancer", "width", "ConcPLW_1.x.1", "ConcPLW_2.x.1", "Fold.x.1", "p.value.x.1", "FDR.x.1",
                                                                      "ConcPLW_1.x.2", "ConcPLW_2.x.2", "Fold.x.2", "p.value.x.2", "FDR.x.2", "FDR.x.2.shrunk",
                                                                      "ConcPLW_1.y.1", "ConcPLW_2.y.1", "Fold.y.1", "p.value.y.1", "FDR.y.1", "FDR.y.1.shrunk",
                                                                      "ConcPLW_1.y.2", "ConcPLW_2.y.2", "Fold.y.2", "p.value.y.2", "FDR.y.2", "FDR.y.2.shrunk")]
  ## add genes and expression, recalculate FDRs
  # gained enhancers
  diff_enhancers_mergegene_enhX_gain = merge(diff_enhancers_merge_enhX_gain, corr_data, by="enhancer", all.x=T, sort=F) # 46713 rows, 4261 genes
  zz_diffexpr.x.1 = diffexpr.x.1[diffexpr.x.1$gene %in% na.omit(diff_enhancers_mergegene_enhX_gain$gene),]
  zz_diffexpr.x.1$expr_padj.x.1.shrunk = p.adjust(zz_diffexpr.x.1$expr_p.x.1, method="BH")
  zz_diffexpr.x.2 = diffexpr.x.2[diffexpr.x.2$gene %in% na.omit(diff_enhancers_mergegene_enhX_gain$gene),]
  zz_diffexpr.x.2$expr_padj.x.2.shrunk = p.adjust(zz_diffexpr.x.2$expr_p.x.2, method="BH")
  zz_diffexpr.y.1 = diffexpr.y.1[diffexpr.y.1$gene %in% na.omit(diff_enhancers_mergegene_enhX_gain$gene),]
  zz_diffexpr.y.1$expr_padj.y.1.shrunk = p.adjust(zz_diffexpr.y.1$expr_p.y.1, method="BH")
  zz_diffexpr.y.2 = diffexpr.y.2[diffexpr.y.2$gene %in% na.omit(diff_enhancers_mergegene_enhX_gain$gene),]
  zz_diffexpr.y.2$expr_padj.y.2.shrunk = p.adjust(zz_diffexpr.y.2$expr_p.y.2, method="BH")
  diff_enhancers_mergegene_enhX_gain = merge(diff_enhancers_mergegene_enhX_gain, zz_diffexpr.x.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 1531 --> 1598
  diff_enhancers_mergegene_enhX_gain = merge(diff_enhancers_mergegene_enhX_gain, zz_diffexpr.x.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 1525 --> 1670
  diff_enhancers_mergegene_enhX_gain = merge(diff_enhancers_mergegene_enhX_gain, zz_diffexpr.y.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 385 --> 393
  diff_enhancers_mergegene_enhX_gain = merge(diff_enhancers_mergegene_enhX_gain, zz_diffexpr.y.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 111 --> 111
  # lost enhancers
  diff_enhancers_mergegene_enhX_lost = merge(diff_enhancers_merge_enhX_lost, corr_data, by="enhancer", all.x=T, sort=F) # 12148 rows, 1117 genes
  zz_diffexpr.x.1 = diffexpr.x.1[diffexpr.x.1$gene %in% na.omit(diff_enhancers_mergegene_enhX_lost$gene),]
  zz_diffexpr.x.1$expr_padj.x.1.shrunk = p.adjust(zz_diffexpr.x.1$expr_p.x.1, method="BH")
  zz_diffexpr.x.2 = diffexpr.x.2[diffexpr.x.2$gene %in% na.omit(diff_enhancers_mergegene_enhX_lost$gene),]
  zz_diffexpr.x.2$expr_padj.x.2.shrunk = p.adjust(zz_diffexpr.x.2$expr_p.x.2, method="BH")
  zz_diffexpr.y.1 = diffexpr.y.1[diffexpr.y.1$gene %in% na.omit(diff_enhancers_mergegene_enhX_lost$gene),]
  zz_diffexpr.y.1$expr_padj.y.1.shrunk = p.adjust(zz_diffexpr.y.1$expr_p.y.1, method="BH")
  zz_diffexpr.y.2 = diffexpr.y.2[diffexpr.y.2$gene %in% na.omit(diff_enhancers_mergegene_enhX_lost$gene),]
  zz_diffexpr.y.2$expr_padj.y.2.shrunk = p.adjust(zz_diffexpr.y.2$expr_p.y.2, method="BH")
  diff_enhancers_mergegene_enhX_lost = merge(diff_enhancers_mergegene_enhX_lost, zz_diffexpr.x.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 766 --> 807
  diff_enhancers_mergegene_enhX_lost = merge(diff_enhancers_mergegene_enhX_lost, zz_diffexpr.x.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 622 --> 697
  diff_enhancers_mergegene_enhX_lost = merge(diff_enhancers_mergegene_enhX_lost, zz_diffexpr.y.1, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 111 --> 119
  diff_enhancers_mergegene_enhX_lost = merge(diff_enhancers_mergegene_enhX_lost, zz_diffexpr.y.2, by="gene", all.x=T, sort=F) # Adjust FDR before merging: genes 35 --> 38
  rm(zz_diffexpr.x.1, zz_diffexpr.x.2, zz_diffexpr.y.1, zz_diffexpr.y.2)
  rm(zz_diff_enhancers.x.1, zz_diff_enhancers.x.1_gain, zz_diff_enhancers.x.1_lost)


  ### Validate tissue enhancers using tissue expression
  if (exp_p_or_q=="q") {
    # tissue gained enhancers
    enhX_gain_expX1 = diff_enhancers_mergegene_enhX_gain[!is.na(diff_enhancers_mergegene_enhX_gain$gene) & !is.na(diff_enhancers_mergegene_enhX_gain$expr_padj.x.1.shrunk) &
                                                           diff_enhancers_mergegene_enhX_gain$expr_padj.x.1.shrunk<.05 & diff_enhancers_mergegene_enhX_gain$expr_fold.x.1>0, ]  # 1598 genes
    enhX_gain_expX = diff_enhancers_mergegene_enhX_gain[!is.na(diff_enhancers_mergegene_enhX_gain$gene) & !is.na(diff_enhancers_mergegene_enhX_gain$expr_padj.x.1.shrunk) &
                                                          diff_enhancers_mergegene_enhX_gain$expr_padj.x.1.shrunk<.05 & diff_enhancers_mergegene_enhX_gain$expr_fold.x.1>0 &
                                                          !is.na(diff_enhancers_mergegene_enhX_gain$expr_padj.x.2.shrunk) &
                                                          diff_enhancers_mergegene_enhX_gain$expr_padj.x.2.shrunk<.05 & diff_enhancers_mergegene_enhX_gain$expr_fold.x.2>0  , ]  # 814 genes
    # tissue lost enhancers
    enhX_lost_expX1 = diff_enhancers_mergegene_enhX_lost[!is.na(diff_enhancers_mergegene_enhX_lost$gene) & !is.na(diff_enhancers_mergegene_enhX_lost$expr_padj.x.1.shrunk) &
                                                           diff_enhancers_mergegene_enhX_lost$expr_padj.x.1.shrunk<.05 & diff_enhancers_mergegene_enhX_lost$expr_fold.x.1<0, ]  # 807 gemes
    enhX_lost_expX = diff_enhancers_mergegene_enhX_lost[!is.na(diff_enhancers_mergegene_enhX_lost$gene) & !is.na(diff_enhancers_mergegene_enhX_lost$expr_padj.x.1.shrunk) &
                                                          diff_enhancers_mergegene_enhX_lost$expr_padj.x.1.shrunk<.05 & diff_enhancers_mergegene_enhX_lost$expr_fold.x.1<0 &
                                                          !is.na(diff_enhancers_mergegene_enhX_lost$expr_padj.x.2.shrunk) &
                                                          diff_enhancers_mergegene_enhX_lost$expr_padj.x.2.shrunk<.05 & diff_enhancers_mergegene_enhX_lost$expr_fold.x.2<0, ]  # 484 genes
  } else if (exp_p_or_q=="p") {
    # tissue gained enhancers
    enhX_gain_expX1 = diff_enhancers_mergegene_enhX_gain[!is.na(diff_enhancers_mergegene_enhX_gain$gene) & !is.na(diff_enhancers_mergegene_enhX_gain$expr_p.x.1) &
                                                           diff_enhancers_mergegene_enhX_gain$expr_p.x.1<.05 & diff_enhancers_mergegene_enhX_gain$expr_fold.x.1>0, ]  # 1849 genes
    enhX_gain_expX = diff_enhancers_mergegene_enhX_gain[!is.na(diff_enhancers_mergegene_enhX_gain$gene) & !is.na(diff_enhancers_mergegene_enhX_gain$expr_p.x.1) &
                                                          diff_enhancers_mergegene_enhX_gain$expr_p.x.1<.05 & diff_enhancers_mergegene_enhX_gain$expr_fold.x.1>0 &
                                                          !is.na(diff_enhancers_mergegene_enhX_gain$expr_p.x.2) &
                                                          diff_enhancers_mergegene_enhX_gain$expr_p.x.2<.05 & diff_enhancers_mergegene_enhX_gain$expr_fold.x.2>0  , ]  # 1041 genes
    # tissue lost enhancers
    enhX_lost_expX1 = diff_enhancers_mergegene_enhX_lost[!is.na(diff_enhancers_mergegene_enhX_lost$gene) & !is.na(diff_enhancers_mergegene_enhX_lost$expr_p.x.1) &
                                                           diff_enhancers_mergegene_enhX_lost$expr_p.x.1<.05 & diff_enhancers_mergegene_enhX_lost$expr_fold.x.1<0, ]  # 820 gemes
    enhX_lost_expX = diff_enhancers_mergegene_enhX_lost[!is.na(diff_enhancers_mergegene_enhX_lost$gene) & !is.na(diff_enhancers_mergegene_enhX_lost$expr_p.x.1) &
                                                          diff_enhancers_mergegene_enhX_lost$expr_p.x.1<.05 & diff_enhancers_mergegene_enhX_lost$expr_fold.x.1<0 &
                                                          !is.na(diff_enhancers_mergegene_enhX_lost$expr_p.x.2) &
                                                          diff_enhancers_mergegene_enhX_lost$expr_p.x.2<.05 & diff_enhancers_mergegene_enhX_lost$expr_fold.x.2<0, ]  # 521 genes
  }


  ### Validate tissue enhancers validated by tissue expr, using cell-line enhancers existence
  # tissue gained enhancers
  enhX_gain_expX_enhY1 = enhX_gain_expX[enhX_gain_expX$ConcPLW_1.y.1 > quantile(ConcPLW_dist.y.1, enh_present_thr), ] # 1001 genes
  enhX_gain_expX_enhY = enhX_gain_expX[enhX_gain_expX$ConcPLW_1.y.1 > quantile(ConcPLW_dist.y.1, enh_present_thr) &
                                         enhX_gain_expX$ConcPLW_1.y.2 > quantile(ConcPLW_dist.y.2, enh_present_thr)  , ] # 990 genes
  # tissue lost enhancers
  enhX_lost_expX_enhY1 = enhX_lost_expX[enhX_lost_expX$ConcPLW_2.y.1 > quantile(ConcPLW_dist.y.1, enh_present_thr), ] # 370 genes
  enhX_lost_expX_enhY = enhX_lost_expX[enhX_lost_expX$ConcPLW_2.y.1 > quantile(ConcPLW_dist.y.1, enh_present_thr) &
                                         enhX_lost_expX$ConcPLW_2.y.2 > quantile(ConcPLW_dist.y.2, enh_present_thr),  ] # 286 genes

  ### print and plot info about ConcPLW distribution
  layout(mat=matrix(1:6, 2, 3, byrow = T))
  hist(ConcPLW_dist.y.1[ConcPLW_dist.y.1 < quantile(ConcPLW_dist.y.1,.99)],breaks=500, xlab="Concentration per log width in Contrast 1, in Y samples")
  abline(v=quantile(ConcPLW_dist.y.1, enh_present_thr))
  smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_1.y.1, xlab="Log enhancer width", ylab="Conc per log width, Contrast 1, Y Group 1", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_1.y.1,.99)))
  abline(h=quantile(ConcPLW_dist.y.1, enh_present_thr))
  smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_2.y.1, xlab="Log enhancer width", ylab="Conc per log width, Contrast 1, Y Group 2", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_2.y.1,.99)))
  abline(h=quantile(ConcPLW_dist.y.1, enh_present_thr))

  hist(ConcPLW_dist.y.2[ConcPLW_dist.y.2 < quantile(ConcPLW_dist.y.2,.99)],breaks=500, xlab="Concentration per log width in Contrast 2, in Y samples")
  abline(v=quantile(ConcPLW_dist.y.2, enh_present_thr))
  smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_1.y.2, xlab="Log enhancer width", ylab="Conc per log width, Contrast 2, Y Group 1", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_1.y.2,.99)))
  abline(h=quantile(ConcPLW_dist.y.2, enh_present_thr))
  smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_2.y.2, xlab="Log enhancer width", ylab="Conc per log width, Contrast 2, Y Group 2", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_2.y.2,.99)))
  abline(h=quantile(ConcPLW_dist.y.2, enh_present_thr))

  print(paste0("Concentration per log width in Contrast 1, in Y samples (n=", length(ConcPLW_dist.y.1), "):"), quote=F)
  print(summary(ConcPLW_dist.y.1))
  print(paste0("enh_present_thr = ", enh_present_thr, ", cutoff = ", quantile(ConcPLW_dist.y.1, enh_present_thr)), quote=F)
  print("", quote=F)
  print(paste0("Concentration per log width in Contrast 2, in Y samples (n=", length(ConcPLW_dist.y.2), "):"), quote=F)
  print(summary(ConcPLW_dist.y.2))
  print(paste0("enh_present_thr = ", enh_present_thr, ", cutoff = ", quantile(ConcPLW_dist.y.2, enh_present_thr)), quote=F)
  print("--------------------------", quote=F)


  # # tissue gained enhancers
  # enhX_gain_expX_enhY1 = enhX_gain_expX[enhX_gain_expX$ConcPLW_1.y.1 > quantile(ConcPLW_dist_1.y.1, enh_present_thr), ] # 1001 genes
  # enhX_gain_expX_enhY = enhX_gain_expX[enhX_gain_expX$ConcPLW_1.y.1 > quantile(ConcPLW_dist_1.y.1, enh_present_thr) &
  #                                        enhX_gain_expX$ConcPLW_1.y.2 > quantile(ConcPLW_dist_1.y.2, enh_present_thr)  , ] # 990 genes
  # # tissue lost enhancers
  # enhX_lost_expX_enhY1 = enhX_lost_expX[enhX_lost_expX$ConcPLW_2.y.1 > quantile(ConcPLW_dist_2.y.1, enh_present_thr), ] # 370 genes
  # enhX_lost_expX_enhY = enhX_lost_expX[enhX_lost_expX$ConcPLW_2.y.1 > quantile(ConcPLW_dist_2.y.1, enh_present_thr) &
  #                                        enhX_lost_expX$ConcPLW_2.y.2 > quantile(ConcPLW_dist_2.y.2, enh_present_thr),  ] # 286 genes
  #
  # ### print and plot info about ConcPLW distribution
  # layout(mat=matrix(1:8, 2, 4, byrow = T))
  # hist(ConcPLW_dist_1.y.1[ConcPLW_dist_1.y.1 < quantile(ConcPLW_dist_1.y.1,.5)],breaks=100, xlab="Concentration per log width in Contrast 1, Y Group 1")
  # abline(v=quantile(ConcPLW_dist_1.y.1, enh_present_thr))
  #
  # hist(ConcPLW_dist_2.y.1[ConcPLW_dist_2.y.1 < quantile(ConcPLW_dist_2.y.1,.5)],breaks=100, xlab="Concentration per log width in Contrast 1, Y Group 2")
  # abline(v=quantile(ConcPLW_dist_2.y.1, enh_present_thr))
  #
  # smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_1.y.1, xlab="Log enhancer width", ylab="Conc per log width, Contrast 1, Y Group 1", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_1.y.1,.9)))
  # abline(h=quantile(ConcPLW_dist_1.y.1, enh_present_thr))
  # smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_2.y.1, xlab="Log enhancer width", ylab="Conc per log width, Contrast 1, Y Group 2", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_2.y.1,.9)))
  # abline(h=quantile(ConcPLW_dist_2.y.1, enh_present_thr))
  #
  # hist(ConcPLW_dist_1.y.2[ConcPLW_dist_1.y.2 < quantile(ConcPLW_dist_1.y.2,.5)],breaks=100, xlab="Concentration per log width in Contrast 2, Y Group 1")
  # abline(v=quantile(ConcPLW_dist_1.y.2, enh_present_thr))
  #
  # hist(ConcPLW_dist_2.y.2[ConcPLW_dist_2.y.2 < quantile(ConcPLW_dist_2.y.2,.5)],breaks=100, xlab="Concentration per log width in Contrast 2, Y Group 2")
  # abline(v=quantile(ConcPLW_dist_2.y.2, enh_present_thr))
  #
  # smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_1.y.2, xlab="Log enhancer width", ylab="Conc per log width, Contrast 2, Y Group 1", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_1.y.2,.9)))
  # abline(h=quantile(ConcPLW_dist_1.y.2, enh_present_thr))
  # smoothScatter(x=log2(diff_enhancers_merge$width), y=diff_enhancers_merge$ConcPLW_2.y.2, xlab="Log enhancer width", ylab="Conc per log width, Contrast 2, Y Group 2", ylim=c(0,quantile(diff_enhancers_merge$ConcPLW_2.y.2,.9)))
  # abline(h=quantile(ConcPLW_dist_2.y.2, enh_present_thr))






  ### Validate tissue enhancers validated by tissue expr and cell-line enhancers existence, using cell-line enhancers gain/lost
  # tissue gained enhancers
  enhX_gain_expX_enhYgain1 = enhX_gain_expX_enhY[enhX_gain_expX_enhY$Fold.y.1 > .2, ] # 926 genes
  enhX_gain_expX_enhYgain = enhX_gain_expX_enhY[enhX_gain_expX_enhY$Fold.y.1 > .2 & enhX_gain_expX_enhY$Fold.y.2 > .2, ] # 744 genes
  # tissue lost enhancers
  enhX_lost_expX_enhYlost1 = enhX_lost_expX_enhY[enhX_lost_expX_enhY$Fold.y.1 < -.2, ] # 250 genes
  enhX_lost_expX_enhYlost = enhX_lost_expX_enhY[enhX_lost_expX_enhY$Fold.y.1 < -.2 & enhX_lost_expX_enhY$Fold.y.2 < -.2, ] # 222 genes


  ### Validate tissue enhancers validated by tissue expr and cell-line enhancers existence and gain/lost, using cell-line expr
  # tissue gained enhancers
  enhX_gain_expX_enhYgain_expY1 =  enhX_gain_expX_enhYgain[!is.na(enhX_gain_expX_enhYgain$expr_fold.y.1) & enhX_gain_expX_enhYgain$expr_fold.y.1 > .1,] # 432 genes
  enhX_gain_expX_enhYgain_expY =  enhX_gain_expX_enhYgain[!is.na(enhX_gain_expX_enhYgain$expr_fold.y.1) & enhX_gain_expX_enhYgain$expr_fold.y.1 > .1 &
                                                            !is.na(enhX_gain_expX_enhYgain$expr_fold.y.2) & enhX_gain_expX_enhYgain$expr_fold.y.2 > .1,] # 199 genes
  # tissue lost enhancers
  enhX_lost_expX_enhYlost_expY1 =  enhX_lost_expX_enhYlost[!is.na(enhX_lost_expX_enhYlost$expr_fold.y.1) & enhX_lost_expX_enhYlost$expr_fold.y.1 < -.1,] # 97 genes
  enhX_lost_expX_enhYlost_expY =  enhX_lost_expX_enhYlost[!is.na(enhX_lost_expX_enhYlost$expr_fold.y.1) & enhX_lost_expX_enhYlost$expr_fold.y.1 < -.1 &
                                                            !is.na(enhX_lost_expX_enhYlost$expr_fold.y.2) & enhX_lost_expX_enhYlost$expr_fold.y.2 < -.1,] # 50 genes


  #### barplot enhancers validation
  layout(mat=matrix(c(1,2), 1, 2))
  enhancer_val = matrix(c(length(unique(diff_enhancers_mergegene_enhX1_gain$enhancer)), length(unique(diff_enhancers_mergegene_enhX1_lost$enhancer)),
                          length(unique(diff_enhancers_mergegene_enhX_gain$enhancer)), length(unique(diff_enhancers_mergegene_enhX_lost$enhancer)),
                          length(unique(diff_enhancers_mergegene_enhX_gain[!is.na(diff_enhancers_mergegene_enhX_gain$gene),"enhancer"])), length(unique(diff_enhancers_mergegene_enhX_lost[!is.na(diff_enhancers_mergegene_enhX_lost$gene),"enhancer"])),
                          length(unique(enhX_gain_expX1$enhancer)), length(unique(enhX_lost_expX1$enhancer)),
                          length(unique(enhX_gain_expX$enhancer)), length(unique(enhX_lost_expX$enhancer)),
                          length(unique(enhX_gain_expX_enhY1$enhancer)), length(unique(enhX_lost_expX_enhY1$enhancer)),
                          length(unique(enhX_gain_expX_enhY$enhancer)), length(unique(enhX_lost_expX_enhY$enhancer)),
                          length(unique(enhX_gain_expX_enhYgain1$enhancer)), length(unique(enhX_lost_expX_enhYlost1$enhancer)),
                          length(unique(enhX_gain_expX_enhYgain$enhancer)), length(unique(enhX_lost_expX_enhYlost$enhancer)),
                          length(unique(enhX_gain_expX_enhYgain_expY1$enhancer)), length(unique(enhX_lost_expX_enhYlost_expY1$enhancer)),
                          length(unique(enhX_gain_expX_enhYgain_expY$enhancer)), length(unique(enhX_lost_expX_enhYlost_expY$enhancer))), 11, 2, byrow=T)
  colnames(enhancer_val) = c("Gained", "Lost")
  rownames(enhancer_val) = c("All 1", "All 2", "Mapped to gene", "Val exprX 1", "Val exprX 2", "Val enhY exist 1", "Val enhY exist 2", "Val enhY change 1", "Val enhY change 2", "Val exprY 1", "Val exprY 2")
  barplot(enhancer_val, beside=T,
          col=c("white", "beige", "grey", "cadetblue1", "cadetblue3", "aquamarine1", "aquamarine3", "blue", "darkblue", "gray16", "black"),
          names.arg=c("Gained enhancers", "Lost enhancers"), main="Enhancers validation", ylab="Enhancers", cex.names = .7, cex.main=.7)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend =c("All 1", "All 2", "Mapped to gene", "Val exprX 1", "Val exprX 2", "Val enhY exist 1", "Val enhY exist 2", "Val enhY change 1", "Val enhY change 2", "Val exprY 1", "Val exprY 2"),
         fill = c("white", "beige", "grey", "cadetblue1", "cadetblue3", "aquamarine1", "aquamarine3", "blue", "darkblue", "gray16", "black"), cex=.7)
  print("", quote=F)
  print("Enhancer validation:", quote=F)
  print(enhancer_val)
  print("--------------------------", quote=F)


  ## barplot gene level validation
  layout(mat=matrix(c(1,2), 1, 2))
  gene_val = matrix(c(length(unique(na.omit(diff_enhancers_mergegene_enhX1_gain$gene))), length(unique(na.omit(diff_enhancers_mergegene_enhX1_lost$gene))),
                      length(unique(na.omit(diff_enhancers_mergegene_enhX_gain$gene))), length(unique(na.omit(diff_enhancers_mergegene_enhX_lost$gene))),
                      length(unique(enhX_gain_expX1$gene)), length(unique(enhX_lost_expX1$gene)),
                      length(unique(enhX_gain_expX$gene)), length(unique(enhX_lost_expX$gene)),
                      length(unique(enhX_gain_expX_enhY1$gene)), length(unique(enhX_lost_expX_enhY1$gene)),
                      length(unique(enhX_gain_expX_enhY$gene)), length(unique(enhX_lost_expX_enhY$gene)),
                      length(unique(enhX_gain_expX_enhYgain1$gene)), length(unique(enhX_lost_expX_enhYlost1$gene)),
                      length(unique(enhX_gain_expX_enhYgain$gene)), length(unique(enhX_lost_expX_enhYlost$gene)),
                      length(unique(enhX_gain_expX_enhYgain_expY1$gene)), length(unique(enhX_lost_expX_enhYlost_expY1$gene)),
                      length(unique(enhX_gain_expX_enhYgain_expY$gene)), length(unique(enhX_lost_expX_enhYlost_expY$gene))), 10, 2, byrow=T)
  colnames(gene_val) = c("Gained", "Lost")
  rownames(gene_val) = c("All 1", "All 2", "Val exprX 1", "Val exprX 2", "Val enhY exist 1", "Val enhY exist 2", "Val enhY change 1", "Val enhY change 2", "Val exprY 1", "Val exprY 2")
  barplot(gene_val, beside=T,
          col=c("white", "beige", "cadetblue1", "cadetblue3", "aquamarine1", "aquamarine3", "blue", "darkblue", "gray16", "black"),
          names.arg=c("Gained enhancers", "Lost enhancers"), main="Enhancers validation (gene level)", ylab="Genes", cex.names = .7, cex.main=.7)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend =c("All 1", "All 2", "Val exprX 1", "Val exprX 2", "Val enhY exist 1", "Val enhY exist 2", "Val enhY change 1", "Val enhY change 2", "Val exprY 1", "Val exprY 2"),
         fill = c("white", "beige", "cadetblue1", "cadetblue3", "aquamarine1", "aquamarine3", "blue", "darkblue", "gray16", "black"), cex=.7)
  print("", quote=F)
  print("Enhancer validation (gene level):", quote=F)
  print(gene_val)
  print("--------------------------", quote=F)


  ### print enriched genesets at each validation step
  ## Gain
  enhX_gain_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(diff_enhancers_mergegene_enhX_gain$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_gain_genesets[enhX_gain_genesets$padj<.05 & enhX_gain_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(diff_enhancers_mergegene_enhX_gain$gene))), " genes, ", sum(enhX_gain_genesets$padj<.05 & enhX_gain_genesets$est>1),  " genesets gained in X. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_gain_expX_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_gain_expX$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_gain_expX_genesets[enhX_gain_expX_genesets$padj<.05 & enhX_gain_expX_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_gain_expX$gene))), " genes, ", sum(enhX_gain_expX_genesets$padj<.05 & enhX_gain_expX_genesets$est>1),  " genesets gained in X, Val exprX. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_gain_expX_enhY_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_gain_expX_enhY$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_gain_expX_enhY_genesets[enhX_gain_expX_enhY_genesets$padj<.05 & enhX_gain_expX_enhY_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_gain_expX_enhY$gene))), " genes, ", sum(enhX_gain_expX_enhY_genesets$padj<.05 & enhX_gain_expX_enhY_genesets$est>1),  " genesets gained in X, Val exprX enhYexist. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_gain_expX_enhYgain_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_gain_expX_enhYgain$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_gain_expX_enhYgain_genesets[enhX_gain_expX_enhYgain_genesets$padj<.05 & enhX_gain_expX_enhYgain_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_gain_expX_enhYgain$gene))), " genes, ", sum(enhX_gain_expX_enhYgain_genesets$padj<.05 & enhX_gain_expX_enhYgain_genesets$est>1),  " genesets gained in X, Val exprX enhYexist enhYgain. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_gain_expX_enhYgain_expY_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_gain_expX_enhYgain_expY$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_gain_expX_enhYgain_expY_genesets[enhX_gain_expX_enhYgain_expY_genesets$padj<.05 & enhX_gain_expX_enhYgain_expY_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_gain_expX_enhYgain_expY$gene))), " genes, ", sum(enhX_gain_expX_enhYgain_expY_genesets$padj<.05 & enhX_gain_expX_enhYgain_expY_genesets$est>1),  " genesets gained in X, Val exprX enhYexist enhYgain exprY. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  print("--------------------------", quote=F)
  ## Lost
  enhX_lost_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(diff_enhancers_mergegene_enhX_lost$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_lost_genesets[enhX_lost_genesets$padj<.05 & enhX_lost_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(diff_enhancers_mergegene_enhX_lost$gene))), " genes, ", sum(enhX_lost_genesets$padj<.05 & enhX_lost_genesets$est>1),  " genesets lost in X. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_lost_expX_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_lost_expX$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_lost_expX_genesets[enhX_lost_expX_genesets$padj<.05 & enhX_lost_expX_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_lost_expX$gene))), " genes, ", sum(enhX_lost_expX_genesets$padj<.05 & enhX_lost_expX_genesets$est>1),  " genesets lost in X, Val exprX. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_lost_expX_enhY_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_lost_expX_enhY$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_lost_expX_enhY_genesets[enhX_lost_expX_enhY_genesets$padj<.05 & enhX_lost_expX_enhY_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_lost_expX_enhY$gene))), " genes, ", sum(enhX_lost_expX_enhY_genesets$padj<.05 & enhX_lost_expX_enhY_genesets$est>1),  " genesets lost in X, Val exprX enhYexist. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_lost_expX_enhYlost_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_lost_expX_enhYlost$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_lost_expX_enhYlost_genesets[enhX_lost_expX_enhYlost_genesets$padj<.05 & enhX_lost_expX_enhYlost_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_lost_expX_enhYlost$gene))), " genes, ", sum(enhX_lost_expX_enhYlost_genesets$padj<.05 & enhX_lost_expX_enhYlost_genesets$est>1),  " genesets lost in X, Val exprX enhYexist enhYlost. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  enhX_lost_expX_enhYlost_expY_genesets = analyze_geneset_enrichment(genes_interest=unique(na.omit(enhX_lost_expX_enhYlost_expY$gene)), genes_univ=unique(corr_data$gene), pathways_list=pathways_list)
  zz = enhX_lost_expX_enhYlost_expY_genesets[enhX_lost_expX_enhYlost_expY_genesets$padj<.05 & enhX_lost_expX_enhYlost_expY_genesets$est>1,]
  zz = head(zz[order(zz$p), ], n=30)
  print("", quote=F)
  print(paste0(length(unique(na.omit(enhX_lost_expX_enhYlost_expY$gene))), " genes, ", sum(enhX_lost_expX_enhYlost_expY_genesets$padj<.05 & enhX_lost_expX_enhYlost_expY_genesets$est>1),  " genesets lost in X, Val exprX enhYexist enhYlost exprY. Top 30:"), quote=F)
  if (nrow(zz)>0) {
    tmp = apply(zz, 1, function(x){ print(paste0(as.character(x["pathway"]), ",  p=", signif(as.numeric(x["p"]),digits=3), ", padj=", signif(as.numeric(x["padj"]),digits=3), ", int=", as.numeric(x["num_int"]), "/", as.numeric(x["num_gs"]), ", genes=", substr(as.character(x["int_genes"]),1,70)), quote=F) })
  }
  print("--------------------------", quote=F)



  return(list(gain=diff_enhancers_mergegene_enhX_gain,
              gain_exp=enhX_gain_expX, gain_exp_enh=enhX_gain_expX_enhY, gain_exp_enhgain=enhX_gain_expX_enhYgain, gain_exp_enhgain_exp=enhX_gain_expX_enhYgain_expY,
              gain_exp_genesets=enhX_gain_expX_genesets, gain_exp_enh_genesets=enhX_gain_expX_enhY_genesets, gain_exp_enhgain_genesets=enhX_gain_expX_enhYgain_genesets, gain_exp_enhgain_exp_genesets=enhX_gain_expX_enhYgain_expY_genesets,
              lost=diff_enhancers_mergegene_enhX_lost,
              lost_exp=enhX_lost_expX, lost_exp_enh=enhX_lost_expX_enhY, lost_exp_enhlost=enhX_lost_expX_enhYlost, lost_exp_enhlost_exp=enhX_lost_expX_enhYlost_expY,
              lost_exp_genesets=enhX_lost_expX_genesets, lost_exp_enh_genesets=enhX_lost_expX_enhY_genesets, lost_exp_enhlost_genesets=enhX_lost_expX_enhYlost_genesets, lost_exp_enhlost_exp_genesets=enhX_lost_expX_enhYlost_expY_genesets))

}
