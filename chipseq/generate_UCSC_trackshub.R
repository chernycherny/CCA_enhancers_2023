

########################################################################################################################
######################################### Generate trackDb.txt file for UCSC genome browser tracks hub #################
########################################################################################################################



############# Tissues

# sample names should correspond to the bam files
trackDb = data.frame(sample = c("21914187T", "29150572T", "A100T", "Z3585T", 
                                 "33587479T", "20020012T", "17231203T", "72600277T", "Z3722T", "A035T", "C080T", "C096T", "W39T", "Y65T", "Y140T", "1202T", "807N", 
                                 "260418N_merged", "3011118N", "4081118N"), stringsAsFactors = F)

trackDb$path = c(rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch1/",4),
                 rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch4/", 13),
                 rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch6/", 3))

trackDb$track = paste0("X", trackDb$sample, "_K27ac_subtract_Deseq2NormFac")

trackDb$bigDataUrl = paste0(trackDb$path, trackDb$sample, "_K27ac_subtract_Deseq2NormFac.filt.bw")

trackDb$shortLabel = paste0(trackDb$sample, "_Deseq")
trackDb$longLabel = paste0(trackDb$sample, "_K27ac_Deseq")
trackDb$type = "bigWig"
trackDb$visibility = "full"
trackDb$maxHeightPixels = 32
trackDb$viewLimits = "0:2"
trackDb$viewLimitsMax = "0:2"
trackDb$gridDefault = "on"

sample_ordering = c("807N", "260418N_merged", "3011118N", "4081118N",
                    "21914187T", "29150572T", "Z3585T", "33587479T", "20020012T", "17231203T", "72600277T", "Z3722T", "1202T", 
                    "A035T", "C080T", "C096T", "W39T", "Y65T", "Y140T", "A100T")
trackDb$priority = match(trackDb$sample, sample_ordering)
trackDb = trackDb[order(trackDb$priority),]

# write file
file.remove("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/ucsc genome browser/myHub/hg19/trackDb.txt")
trackDb_list = split(trackDb, trackDb$priority)
trackDb_list = lapply(trackDb_list, function(x){
  y = t(x[, c("track", "bigDataUrl", "shortLabel", "longLabel", "type", "visibility", "maxHeightPixels", "viewLimits", "viewLimitsMax", "gridDefault", "priority")])
  y = rbind(y, "")
  write.table(y, file="C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/ucsc genome browser/myHub/hg19/trackDb.txt", append=T, quote=F, sep="\t", col.names=F)
})





############# Cell-lines

# sample names should correspond to the bam files
trackDb = data.frame(sample = c("21914187T", "29150572T", "A100T", "Z3585T", 
                                "33587479T", "20020012T", "17231203T", "72600277T", "Z3722T", "A035T", "C080T", "C096T", "W39T", "Y65T", "Y140T", "1202T", "807N", 
                                "260418N_merged", "3011118N", "4081118N"), stringsAsFactors = F)

trackDb$path = c(rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch1/",4),
                 rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch4/", 13),
                 rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch6/", 3))

trackDb$track = paste0("X", trackDb$sample, "_K27ac_subtract_Deseq2NormFac")

trackDb$bigDataUrl = paste0(trackDb$path, trackDb$sample, "_K27ac_subtract_Deseq2NormFac.filt.bw")

trackDb$shortLabel = paste0(trackDb$sample, "_Deseq")
trackDb$longLabel = paste0(trackDb$sample, "_K27ac_Deseq")
trackDb$type = "bigWig"
trackDb$visibility = "full"
trackDb$maxHeightPixels = 32
trackDb$viewLimits = "0:2"
trackDb$viewLimitsMax = "0:2"
trackDb$gridDefault = "on"

sample_ordering = c("807N", "260418N_merged", "3011118N", "4081118N",
                    "21914187T", "29150572T", "Z3585T", "33587479T", "20020012T", "17231203T", "72600277T", "Z3722T", "1202T", 
                    "A035T", "C080T", "C096T", "W39T", "Y65T", "Y140T", "A100T")
trackDb$priority = match(trackDb$sample, sample_ordering)
trackDb = trackDb[order(trackDb$priority),]

# write file
file.remove("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/ucsc genome browser/myHub/hg19/trackDb.txt")
trackDb_list = split(trackDb, trackDb$priority)
trackDb_list = lapply(trackDb_list, function(x){
  y = t(x[, c("track", "bigDataUrl", "shortLabel", "longLabel", "type", "visibility", "maxHeightPixels", "viewLimits", "viewLimitsMax", "gridDefault", "priority")])
  y = rbind(y, "")
  write.table(y, file="C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/ucsc genome browser/myHub/hg19/trackDb.txt", append=T, quote=F, sep="\t", col.names=F)
})




########## SPMR-normalized 
trackDb = data.frame(sample = c("21914187T", "29150572T", "A100T", "Z3585T", 
                                "33587479T", "20020012T", "17231203T", "72600277T", "Z3722T", "A035T", "C080T", "C096T", "W39T", "Y65T", "Y140T", "1202T", "807N", 
                                "260418N_merged", "3011118N", "4081118N"), stringsAsFactors = F)

trackDb$path = c(rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch1/",4),
                 rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch4/", 13),
                 rep("https://de.cyverse.org/anon-files/iplant/home/cherny/ucsc_data/signal_tracks_batch6/", 3))

trackDb$track = paste0("X", trackDb$sample, "_K27ac_subtract_SPMR")

trackDb$bigDataUrl = paste0(trackDb$path, trackDb$sample, "_K27ac_subtract_SPMR.filt.bw")

trackDb$shortLabel = paste0(trackDb$sample, "_SPMR")
trackDb$longLabel = paste0(trackDb$sample, "_K27ac_SPMR")
trackDb$type = "bigWig"
trackDb$visibility = "hide"
trackDb$maxHeightPixels = 32
trackDb$viewLimits = "0:2"
trackDb$viewLimitsMax = "0:2"
trackDb$gridDefault = "on"

sample_ordering = c("807N", "260418N_merged", "3011118N", "4081118N",
                    "21914187T", "29150572T", "Z3585T", "33587479T", "20020012T", "17231203T", "72600277T", "Z3722T", "1202T", 
                    "A035T", "C080T", "C096T", "W39T", "Y65T", "Y140T", "A100T")
trackDb$priority = match(trackDb$sample, sample_ordering)
trackDb = trackDb[order(trackDb$priority),]

# write file
file.remove("C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/ucsc genome browser/myHub_SPMR/hg19/trackDb.txt")
trackDb_list = split(trackDb, trackDb$priority)
trackDb_list = lapply(trackDb_list, function(x){
  y = t(x[, c("track", "bigDataUrl", "shortLabel", "longLabel", "type", "visibility", "maxHeightPixels", "viewLimits", "viewLimitsMax", "gridDefault", "priority")])
  y = rbind(y, "")
  write.table(y, file="C:/Users/cherny/Google Drive/CCA-epigenetics/CCA_enhancers_project/ucsc genome browser/myHub_SPMR/hg19/trackDb.txt", append=T, quote=F, sep="\t", col.names=F)
})




