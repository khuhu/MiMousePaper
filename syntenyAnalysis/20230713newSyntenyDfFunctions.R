### make hg19 to mm10 format for synteny 
### get mm10 and hg19 cyto Df

# sytenyPortalHg19Mm10 <- readLines("/mnt/DATA6/kevin_recovery/apps/circos/syntenyPortalhg19mm10Results.txt")
# hg19Str <- sytenyPortalHg19Mm10[grep("hg19\\.", sytenyPortalHg19Mm10)]
# hg19Str <- str_remove(unlist(lapply(str_split(hg19Str, "\\."), "[[", 2)), " \\+")
# hg19Str <- str_remove(hg19Str, " \\-")
# hg19Range <- unlist(lapply(str_split(hg19Str, ":"), "[[",  2))
# mm10Str <- sytenyPortalHg19Mm10[grep("mm10\\.", sytenyPortalHg19Mm10)]
# mm10Str <- str_remove(unlist(lapply(str_split(mm10Str, "\\."), "[[", 2)), " \\+")
# mm10Str <- str_remove(mm10Str, " \\-")
# mm10Range <- unlist(lapply(str_split(mm10Str, ":"), "[[",  2))
# sytenyPortalDf <- data.frame("ref_chr" = str_remove(unlist(lapply(str_split(hg19Str, ":"), "[[",  1)), "chr"),
#                              "ref_start_pos" = as.numeric(unlist(lapply(str_split(hg19Range, "-"), "[[",  1))),
#                              "ref_end_pos" = as.numeric(unlist(lapply(str_split(hg19Range, "-"), "[[",  2))),
#                              "comp_chr" = str_remove(unlist(lapply(str_split(mm10Str, ":"), "[[",  1)), "chr"),
#                              "comp_start_pos" = as.numeric(unlist(lapply(str_split(mm10Range, "-"), "[[",  1))),
#                              "comp_end_pos" = as.numeric(unlist(lapply(str_split(mm10Range, "-"), "[[",  2))))
# write.table(sytenyPortalDf, "/mnt/DATA6/kevin_recovery/apps/circos/20230713_mm10_hg19_syntenyDf.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 
# 
# cytohg19 <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/hg19/cytoBand.txt", sep = "\t", header = FALSE)
# cytohg19$V1 <- paste0("h_", cytohg19$V1)
# 
# cytomm10 <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/mm10/cytoBand.txt", sep = "\t", header = FALSE)
# cytomm10$V1 <- paste0("m_", cytomm10$V1)
# 
# hg19mm10cyto <- rbind(cytohg19, cytomm10)
# write.table(hg19mm10cyto, "/mnt/DATA6/kevin_recovery/apps/circos/20230713_mm10_hg19_cytoDf.txt",
#             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 



synteny_hg19mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20230713_mm10_hg19_syntenyDf.txt",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE)

cyto_hg19mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20230713_mm10_hg19_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)


# segInput <- tcga_amp


### ferqs are wrong use ref and comp inputs for correct mapping

# syntenyPlotInputsFreqV3 <- function(segInput, species = "human"){
#   # maps synteny starting from human
#   
#   tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$Chr),
#                                ranges = IRanges(start = as.numeric(segInput$Start),
#                                                 end = as.numeric(segInput$End)))
#   
#   if (species == "human") {
#     synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg19mm10$ref_chr),
#                               ranges = IRanges(start = synteny_hg19mm10$ref_start_pos,
#                                                end = synteny_hg19mm10$ref_end_pos))
#     res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
#     query_hits <- queryHits(res_overlap)
#     subject_hits <- subjectHits(res_overlap)
#     
#     mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg19mm10$comp_chr[query_hits]),
#                             "start" = synteny_hg19mm10$comp_start_pos[query_hits],
#                             "end" = synteny_hg19mm10$comp_end_pos[query_hits],
#                             "freq" = segInput$Freq[query_hits],
#                             stringsAsFactors = FALSE)
#     
#     human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg19mm10$ref_chr[query_hits]),
#                             "start" = synteny_hg19mm10$ref_start_pos[query_hits],
#                             "end" = synteny_hg19mm10$ref_end_pos[query_hits],
#                             stringsAsFactors = FALSE)
#     
#   } else if(species == "mouse"){
#     synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg19mm10$comp_chr),
#                               ranges = IRanges(start = synteny_hg19mm10$comp_start_pos,
#                                                end = synteny_hg19mm10$comp_end_pos))
#     res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
#     query_hits <- queryHits(res_overlap)
#     subject_hits <- subjectHits(res_overlap)
#     
#     mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg19mm10$comp_chr[query_hits]),
#                             "start" = synteny_hg19mm10$comp_start_pos[query_hits],
#                             "end" = synteny_hg19mm10$comp_end_pos[query_hits],
#                             "freq" = segInput$Freq[query_hits],
#                             stringsAsFactors = FALSE)
#     
#     human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg19mm10$ref_chr[query_hits]),
#                             "start" = synteny_hg19mm10$ref_start_pos[query_hits],
#                             "end" = synteny_hg19mm10$ref_end_pos[query_hits],
#                             stringsAsFactors = FALSE)
#     
#   }
#   
# 
#   resList <- list("human_bed" = human_bed,
#                   "mouse_bed" = mouse_bed)
#   
#   return(resList)
# }

# m_amp <- allMouseAneu_bprn_amp_bed
# m_del <- allMouseAneu_bprn_del_bed
# tcga_amp <- armGisticTp53_amp_bed
# tcga_del <- armGisticTp53_del_bed
# geneVar <- "no"
# ref <- "human"
# reduction <- TRUE


circosFreq2 <- function(m_amp, m_del, tcga_amp, tcga_del, filename = "test",
                        geneVar = "no", ref = "human", empty = FALSE, plot = TRUE){
  # 20210802 weird upstream bug of 1bp positions with start > end
  if (length(which(m_amp$Start - m_amp$End > 0)) > 0) {
    m_amp <- m_amp[-which(m_amp$Start - m_amp$End > 0),]
    print(1)
  }
  
  if (length(which(m_del$Start - m_del$End > 0)) > 0) {
    m_del <- m_del[-which(m_del$Start - m_del$End > 0),]
    print(2)
  }
  
  if (length(which(tcga_amp$Start - tcga_amp$End > 0)) > 0) {
    tcga_amp <- tcga_amp[-which(tcga_amp$Start - tcga_amp$End > 0),]
    print(3)
  }
  
  if (length(which(tcga_del$Start - tcga_del$End > 0)) > 0) {
    tcga_del <- tcga_del[-which(tcga_del$Start - tcga_del$End > 0),]
    print(4)
  }
  
  
  
  tcga_del$Freq  <- tcga_del$Freq * -1
  m_del$Freq <- m_del$Freq * -1
  if (geneVar == "yes") {
    tcga_del$Start <- tcga_del$Start - 500000
    tcga_del$End <- tcga_del$End + 500000
    
    tcga_amp$Start <- tcga_amp$Start - 500000
    tcga_amp$End <- tcga_amp$End + 500000
  }
  
  
  mouse_bed <- NULL
  human_bed <- NULL
  
  if (ref == "human") {
    if (nrow(m_amp) > 0) {
      ampSynteny <- syntenyPlotInputsFreqV4(tcga_amp, m_amp, species = "human")
      mouse_bed <- rbind(mouse_bed, ampSynteny$mouse_bed)
      human_bed <- rbind(human_bed, ampSynteny$human_bed)
    }
    
    if(nrow(m_del) >  0){
      delSynteny <- syntenyPlotInputsFreqV4(tcga_del, m_del, species = "human")
      mouse_bed <- rbind(mouse_bed, delSynteny$mouse_bed)
      human_bed <- rbind(human_bed, delSynteny$human_bed)
    }
  } else if(ref == "mouse"){
    if (nrow(m_amp) > 0) {
      ampSynteny <- syntenyPlotInputsFreqV4(m_amp, tcga_amp, species = "mouse")
      mouse_bed <- rbind(mouse_bed, ampSynteny$mouse_bed)
      human_bed <- rbind(human_bed, ampSynteny$human_bed)
    }
    
    if(nrow(m_del) >  0){
      delSynteny <- syntenyPlotInputsFreqV4(m_del, tcga_del, species = "mouse")
      mouse_bed <- rbind(mouse_bed, delSynteny$mouse_bed)
      human_bed <- rbind(human_bed, delSynteny$human_bed)
    }
  }

  

  allSynTable <- cbind(human_bed, mouse_bed)
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "h_freq","m_chr",
                             "m_start", "m_end", "m_freq")
  
  
  # cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
  cyto_interest <- c(paste0("m_chr", 1:19), paste0("h_chr", 1:22))
  cyto_combined_red <- cyto_hg19mm10[which(cyto_hg19mm10$V1 %in% cyto_interest),]
  
  mouse_cn <- trackBedColumns(rbind(m_amp,
                                    m_del))
  mouse_cn$chr <- paste0("m_chr", mouse_cn$chr)
  
  human_cn <- trackBedColumns(rbind(tcga_amp,
                                    tcga_del))
  human_cn$chr <- paste0("h_chr", human_cn$chr)
  
  cn_track_bed <- rbind(mouse_cn, human_cn)
  cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]
  
  track_color <- rep("#000000", nrow(cn_track_bed))
  track_color <- ifelse(cn_track_bed$freq < 0 , "#00008B", "#8B0000") 
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(mouse_bed$chr)
  human_chrs <- unique(human_bed$chr)
  
  colorVector2 <- rep("#000000", nrow(mouse_bed))
  for (j in seq_along(mouse_chrs)) {
    colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
  }
  
  
  ### adds black zero line for all chromosomes
  for (i in unique(cyto_hg19mm10$V1)) {
    if (grepl("X", i)) {
      next()
    } else if(grepl("Y", i)){
      next()
    } else{
      tmp <- cyto_hg19mm10[which(cyto_hg19mm10$V1 == i),]
      cn_track_bed <- rbind(cn_track_bed, data.frame("chr"= tmp$V1[1], "start" = 0, 
                                                     "end" = max(tmp$V3), "freq" = 0, 
                                                     "ybot" = -0.03, "ytop" = 0.03, "col" = "#000000"))
    }
  }
  
  # cn_track_bed$chr <- factor(cn_track_bed$chr, levels = c(paste0("m_chr", 1:19), paste0("h_chr", 1:22)))
  if(plot == TRUE){
    if (empty == TRUE) {
      pdf(paste0("/mnt/DATA5/tmp/kev/misc/", filename, ".pdf"),
          width = 10, height = 10, useDingbats = FALSE)
      circos.par("track.height"= 0.10) +
        circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
                                      plotType = c("ideogram", "labels"),
                                      chromosome.index = cyto_interest) +
        circos.genomicLink(mouse_bed, human_bed,
                           col = colorVector2,
                           border = NA)
      dev.off()
    } else{
      pdf(paste0("/mnt/DATA5/tmp/kev/misc/", filename, ".pdf"),
          width = 10, height = 10, useDingbats = FALSE)
      circos.par("track.height"= 0.10) +
        circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
                                      plotType = c("ideogram", "labels"),
                                      chromosome.index = cyto_interest) +
        circos.genomicTrack(cn_track_bed, ylim = c(-1.00,1.00),
                            panel.fun = function(region, value, ...) {
                              circos.genomicRect(region, value, col = value$col,
                                                 ytop = value$ytop,
                                                 ybottom = value$ybot,
                                                 border = NA,...)
                            }) +
        circos.genomicLink(mouse_bed, human_bed,
                           col = colorVector2,
                           border = NA)
      dev.off()
      return(allSynTable)
    }
  } else{
    return(allSynTable)
  }
}


### below are functions looking to reduce the output of all aneuploidy to genomic areas where changes occur in the same direction
###
###


# refInput <- tcga_amp
# compInput <- m_amp
# speicies <- "human"

syntenyPlotInputsFreqV4 <- function(refInput, compInput, species = "human"){
  # maps synteny starting from human

  tumor_freq_ranges <- GRanges(seqnames = paste0("chr", refInput$Chr),
                               ranges = IRanges(start = as.numeric(refInput$Start),
                                                end = as.numeric(refInput$End)))

  comp_freq_ranges <- GRanges(seqnames = paste0("chr", compInput$Chr),
                              ranges = IRanges(start = as.numeric(compInput$Start),
                                               end = as.numeric(compInput$End)))


  if (species == "human") {
    synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg19mm10$ref_chr),
                              ranges = IRanges(start = synteny_hg19mm10$ref_start_pos,
                                               end = synteny_hg19mm10$ref_end_pos))
    res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
    query_hits <- queryHits(res_overlap)
    subject_hits <- subjectHits(res_overlap)
    
    ### note duplicated regions i.e single query with two subject hits might be interesting
    
    human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg19mm10$ref_chr[query_hits]),
                            "start" = synteny_hg19mm10$ref_start_pos[query_hits],
                            "end" = synteny_hg19mm10$ref_end_pos[query_hits],
                            "freq" = refInput$Freq[subject_hits],
                            stringsAsFactors = FALSE)

    
    
    mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg19mm10$comp_chr[query_hits]),
                            "start" = synteny_hg19mm10$comp_start_pos[query_hits],
                            "end" = synteny_hg19mm10$comp_end_pos[query_hits],
                            stringsAsFactors = FALSE)
    
    mouse_bedGr <- GRanges(seqnames = paste0("chr", synteny_hg19mm10$comp_chr[query_hits]),
                           IRanges(start = synteny_hg19mm10$comp_start_pos[query_hits],
                                   end = synteny_hg19mm10$comp_end_pos[query_hits]))
    
    mouse_bed$freq <- 0
    mouse_bed$freq[queryHits(findOverlaps(mouse_bedGr, comp_freq_ranges))] <- compInput$Freq[subjectHits(findOverlaps(mouse_bedGr, comp_freq_ranges))]
  
    } else if(species == "mouse"){
    synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg19mm10$comp_chr),
                              ranges = IRanges(start = synteny_hg19mm10$comp_start_pos,
                                               end = synteny_hg19mm10$comp_end_pos))
    res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
    query_hits <- queryHits(res_overlap)
    subject_hits <- subjectHits(res_overlap)



    mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg19mm10$comp_chr[query_hits]),
                            "start" = synteny_hg19mm10$comp_start_pos[query_hits],
                            "end" = synteny_hg19mm10$comp_end_pos[query_hits],
                            "freq" = refInput$Freq[subject_hits],
                            stringsAsFactors = FALSE)

    human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg19mm10$ref_chr[query_hits]),
                            "start" = synteny_hg19mm10$ref_start_pos[query_hits],
                            "end" = synteny_hg19mm10$ref_end_pos[query_hits],
                            stringsAsFactors = FALSE)
    
    human_bedGr <- GRanges(seqnames = paste0("chr", synteny_hg19mm10$ref_chr[query_hits]),
                           IRanges(start = synteny_hg19mm10$ref_start_pos[query_hits],
                                   end = synteny_hg19mm10$ref_end_pos[query_hits]))
    
    human_bed$freq <- 0
    human_bed$freq[queryHits(findOverlaps(human_bedGr, comp_freq_ranges))] <- compInput$Freq[subjectHits(findOverlaps(human_bedGr, comp_freq_ranges))]
  }


  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)

  return(resList)
}

### takes filtered allSynTable and plots
# allSynTable <- hgsc_bprn_arm_allSynTable2

circosFreqSimple <- function(allSynTable,  filename = "test")
{
  
  m_amp <- allSynTable[which(allSynTable$m_freq > 0), c("m_chr", "m_start", "m_end", "m_freq")]
  m_del <- allSynTable[which(allSynTable$m_freq < 0), c("m_chr", "m_start", "m_end", "m_freq")]
  tcga_amp <- allSynTable[which(allSynTable$h_freq > 0), c("h_chr", "h_start", "h_end", "h_freq")]
  tcga_del <- allSynTable[which(allSynTable$h_freq < 0), c("h_chr", "h_start", "h_end", "h_freq")]
  
  colnames(m_amp) <- c("Chr", "Start", "End", "Freq")
  colnames(m_del) <- c("Chr", "Start", "End", "Freq")
  colnames(tcga_amp) <- c("Chr", "Start", "End", "Freq")
  colnames(tcga_del) <- c("Chr", "Start", "End", "Freq")
  
  # cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
  cyto_interest <- c(paste0("m_chr", 1:19), paste0("h_chr", 1:22))
  cyto_combined_red <- cyto_hg19mm10[which(cyto_hg19mm10$V1 %in% cyto_interest),]
  
  mouse_cn <- trackBedColumns(rbind(m_amp,
                                    m_del))
  
  human_cn <- trackBedColumns(rbind(tcga_amp,
                                    tcga_del))
  
  cn_track_bed <- rbind(mouse_cn, human_cn)
  cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]
  
  track_color <- rep("#000000", nrow(cn_track_bed))
  track_color <- ifelse(cn_track_bed$freq < 0 , "#00008B", "#8B0000") 
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(allSynTable$m_chr)
  human_chrs <- unique(allSynTable$h_chr)
  
  colorVector2 <- rep("#000000", nrow(mouse_cn[,1:3]))
  for (j in seq_along(mouse_chrs)) {
    colorVector2[which(mouse_cn$chr == mouse_chrs[j])] <- colorVector[j]
  }
  
  
  ### adds black zero line for all chromosomes
  for (i in unique(cyto_hg19mm10$V1)) {
    if (grepl("X", i)) {
      next()
    } else if(grepl("Y", i)){
      next()
    } else{
      tmp <- cyto_hg19mm10[which(cyto_hg19mm10$V1 == i),]
      cn_track_bed <- rbind(cn_track_bed, data_frame("chr"= tmp$V1[1], "start" = 0, 
                                                     "end" = max(tmp$V3), "freq" = 0, 
                                                     "ybot" = -0.03, "ytop" = 0.03, "col" = "#000000"))
    }
  }
  
  # cn_track_bed$chr <- factor(cn_track_bed$chr, levels = c(paste0("m_chr", 1:19), paste0("h_chr", 1:22)))
  
  
  pdf(paste0("/mnt/DATA5/tmp/kev/misc/", filename, ".pdf"),
      width = 10, height = 10, useDingbats = FALSE)
  circos.par("track.height"= 0.10) +
    circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
                                  plotType = c("ideogram", "labels"),
                                  chromosome.index = cyto_interest) +
    circos.genomicTrack(cn_track_bed, ylim = c(-1.00,1.00),
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, col = value$col,
                                             ytop = value$ytop,
                                             ybottom = value$ybot,
                                             border = NA,...)
                        }) +
    circos.genomicLink(mouse_cn[,1:3], human_cn[,1:3],
                       col = colorVector2,
                       border = NA)
  dev.off()
}


### old reducing was over complicated and got some indices wrong 


reducingFreqBed2 <- function(df) {
  
  ### regions with zero frequency may have wrong positions, but it get's filtered out later
  ### script correctly obtains regions with same freq
  reducedDf <- NULL
  
  if (nrow(df) == 1) {
    reducedDf <- data.frame("Chr" = df$chr[1], "Start" = df$pos[1], "End" = df$pos[2], "Freq" = 0)
    return(reducedDf)
  } else {
    for (i in 2:nrow(df)) {
      if (df$pos[i] - df$pos[i - 1] > 0) {
        ### continuting from chromsome. including first indice
        tmpVec <- c("Chr" = df$chr[i], "Start" = df$pos[i - 1], "End" = df$pos[i], "Freq" = df[,3][i])
        reducedDf <- rbind(reducedDf, tmpVec)
      } else if(df$pos[i] - df$pos[i - 1] < 0) {
        ### start of new chromosome
        next()
      }
    }
  }
  reducedDf <- data.frame(reducedDf, stringsAsFactors = FALSE)
  return(reducedDf)
}


### enrichment function 2023
### reducing for only synteny result inputs + hypergeometric test
# syntenyGains <- hgsc_bprn_arm_allSynTable_amp
# syntenyLosses <- hgsc_bprn_arm_allSynTable_del
# pathwayList <- human_hallmarks_list

pathwayHyper <- function(syntenyGains, syntenyLosses, pathwayList){ 
  # only using overlap b/c the frequency + amplitudes was
  # already to be found significant in permutation
  
  pathwayRes <- NULL
  
  gainsGrange <- GRanges(seqnames = syntenyGains$h_chr,
                         IRanges(start = syntenyGains$h_start,
                                 end = syntenyGains$h_end))
  
  lossGrange <- GRanges(seqnames = syntenyLosses$h_chr,
                        IRanges(start = syntenyLosses$h_start,
                                end = syntenyLosses$h_end))
  
  homologDfGrange <- GRanges(seqnames = homologDf$h_chr,
                             IRanges(start = homologDf$h_start,
                                     end = homologDf$h_end))
  
  
  gene_symbol_gains <- homologDf$human_gene[subjectHits(findOverlaps(gainsGrange, homologDfGrange))]
  gene_symbol_losses <- homologDf$human_gene[subjectHits(findOverlaps(lossGrange, homologDfGrange))]
  
  drawnGainL <- length(gene_symbol_gains)
  drawnLossL <- length(gene_symbol_losses)
  allGenesL <- nrow(homologDf)
  
  for(i in 1:length(pathwayList)){
    
    tmpHallList <- unlist(pathwayList[[i]])
    tmpGainL <- length(which(gene_symbol_gains %in% tmpHallList))
    tmpLossL <- length(which(gene_symbol_losses %in% tmpHallList))
    
    gainHyper <- phyper(q = tmpGainL - 1, m  = length(tmpHallList),
                        n = allGenesL - length(tmpHallList), k = drawnGainL,
                        lower.tail = FALSE)
    lossHyper <- phyper(q = tmpLossL - 1, m  = length(tmpHallList),
                        n = allGenesL - length(tmpHallList), k = drawnLossL,
                        lower.tail = FALSE)

    pathw <- names(human_hallmarks_list[[i]])
    pathwayRes <- rbind(pathwayRes, data.frame("pathway" = pathw, 
                                               "gainsP" = gainHyper,
                                               "lossP" = lossHyper))
  }
  pathwayRes$gainsAdj <- p.adjust(pathwayRes$gainsP, method = "BH",
                                  n = length(pathwayRes$gainsP) + length(pathwayRes$lossP))
  pathwayRes$lossesAdj <- p.adjust(pathwayRes$lossP, method = "BH",
                                  n = length(pathwayRes$gainsP) + length(pathwayRes$lossP))

  # pathwayRes$gainsAdj <- p.adjust(pathwayRes$gainsP, method = "BH",
  #                                 n = length(pathwayRes$gainsP))
  # pathwayRes$lossesAdj <- p.adjust(pathwayRes$lossP, method = "BH",
  #                                  n = length(pathwayRes$lossP))
  
  return(pathwayRes)
}


# gainsDf <- armGisticTp53_amp_bed
# lossesDf <- armGisticTp53_del_bed
# speciesType <- "human"

freqPlotv2 <- function(gainsDf, lossesDf, speciesType = "human", main = "no title"){
  require(ggplot2)
  
  if (speciesType == "human") {
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    dummyPoints <- read.table("/mnt/DATA5/tmp/kev/misc/20210803hg19_dummyPoints.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(speciesType == "mouse"){
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    dummyPoints <- read.table("/mnt/DATA5/tmp/kev/misc/20210803mm10_dummyPoints.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  tmpDel <- lossesDf
  tmpDel$Freq <- tmpDel$Freq * -1
  tmpDel$Start <- tmpDel$Start/1e6
  tmpDel$End <- tmpDel$End/1e6
  tmpDel$col <- "#0000FF"
  
  tmpAmp <- gainsDf
  tmpAmp$Start <- tmpAmp$Start/1e6
  tmpAmp$End <- tmpAmp$End/1e6
  tmpAmp$col <- "#FF0000"
  
  
  tmpAmpDel_graph <- rbind(tmpAmp, tmpDel)
  for (i in unique(tmpAmpDel_graph$Chr)) {
    tmpAmpDel_graph$Start[which(tmpAmpDel_graph$Chr == i)] <- tmpAmpDel_graph$Start[which(tmpAmpDel_graph$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    tmpAmpDel_graph$End[which(tmpAmpDel_graph$Chr == i)] <- tmpAmpDel_graph$End[which(tmpAmpDel_graph$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  tmpAmpDel_graph2 <- NULL
  for (i in 1:nrow(tmpAmpDel_graph)) {
    tmpVec1 <- c(tmpAmpDel_graph$Chr[i], tmpAmpDel_graph$Start[i], tmpAmpDel_graph$Freq[i], tmpAmpDel_graph$col[i])
    tmpVec2 <- c(tmpAmpDel_graph$Chr[i], tmpAmpDel_graph$End[i], tmpAmpDel_graph$Freq[i], tmpAmpDel_graph$col[i])
    tmpAmpDel_graph2 <- rbind(tmpAmpDel_graph2, tmpVec1, tmpVec2)
  }
  tmpAmpDel_graph2 <- data.frame(tmpAmpDel_graph2, stringsAsFactors = FALSE)
  tmpAmpDel_graph2[,1:3] <- lapply(tmpAmpDel_graph2[,1:3], as.numeric)
  colnames(tmpAmpDel_graph2) <- c("chrom", "pos", "freq", "col")
  tmpAmpDel_graph2$freq <- tmpAmpDel_graph2$freq/100
  
  tmpAmpDel_graph2 <- rbind(tmpAmpDel_graph2, dummyPoints)
  
  tmpAmpDel_graph2 <- tmpAmpDel_graph
  tmpAmpDel_graph2$Freq <- tmpAmpDel_graph2$Freq/100
  
  if (speciesType == "human") {
    chromBreakNoX <- chromBreak[1:23]
    chromTextdfNoX <- chromTextdf[1:22,]
  } else if(speciesType == "mouse"){
    chromBreakNoX <- chromBreak[1:20]
    chromTextdfNoX <- chromTextdf[1:19,]
  }
  
  ggplot(tmpAmpDel_graph2) + geom_vline(xintercept=chromBreakNoX) + geom_hline(yintercept = 0) + 
    geom_rect(aes(xmin = Start, xmax = End, ymin = 0, ymax=ifelse(Freq>0, Freq, 0)), fill="red") +
    geom_rect(aes(xmin = Start, xmax = End, ymax = 0, ymin=ifelse(Freq<0, Freq,0)), fill="blue") + theme_bw() +
    ylim(c(-1,1)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    geom_text(data = chromTextdfNoX, aes(x = xpos, y = ypos, label = chrom), size = 2) + ggtitle(main)
  
}


# gisticArm <- broadBySampleGisticMelt2
# gisticFocal <- gisticOvFocal
fgaCalculator_tcgaV2 <- function(gisticArm , gisticFocal){
  tmpRes <- NULL
  gisticFocal$sampleID <- str_replace_all(gisticFocal$sampleID, "\\-", "\\.")
  ### for each sample scan the aneuploidy status file
  ### have 2 vectors of base pairs for each chromosome arm, multiple
  tmpIdx <- which(gisticArm$sampleID == unique(gisticArm$sampleID)[1])
  armGrange <- GRanges(seqnames = gisticArm$chrom[tmpIdx],
                       IRanges(start = gisticArm$start.pos[tmpIdx], end = gisticArm$end.pos[tmpIdx]))
  for (i in unique(gisticArm$sampleID)) {
    tmpDfArm <- gisticArm[which(gisticArm$sampleID == i),]
    tmpDfFocal <- gisticFocal[which(gisticFocal$sampleID ==i), ]
    tmpDfFocal <- tmpDfFocal[which(abs(tmpDfFocal$mean) > 0.2), ]
    
    tmpIdx2 <- which(abs(tmpDfArm$mean) < 0.2)
    armSum <- sum(tmpDfArm$length[which(abs(tmpDfArm$mean) > 0.2)])
    
    focalGrange <- GRanges(seqnames = tmpDfFocal$chrom,
                           IRanges(start = tmpDfFocal$start.pos, end = tmpDfFocal$end.pos))
    armGrange <- GRanges(seqnames = tmpDfArm$chrom[tmpIdx2],
                          IRanges(start = tmpDfArm$start.pos[tmpIdx2], end = tmpDfArm$end.pos[tmpIdx2]))
    
    focalSum <- 0
    for (j in 1:length(armGrange)) {
      tmpSum <- sum(tmpDfFocal$length[subjectHits(findOverlaps(armGrange[j], focalGrange))])
      focalSum <- focalSum + tmpSum
    }
    tmpFga <- (armSum + focalSum)/sum(tmpDfArm$length)
    tmpRes <- rbind(tmpRes, data.frame("sample" = i, "fga" = tmpFga))
  }
  return(tmpRes)
}


# ampliconArm <- allMouseAneu_bprn
# ampliconFocal <- allMouseCna_bprn
fgaCalculator_ampV2 <- function(ampliconArm , ampliconFocal){
  tmpRes <- NULL
  ampliconFocal$sampleID <- str_replace_all(ampliconFocal$sampleID, "\\-", "\\.")
  ### for each sample scan the aneuploidy status file
  ### have 2 vectors of base pairs for each chromosome arm, multiple
  tmpIdx <- which(ampliconArm$sampleID == unique(ampliconArm$sampleID)[1])
  armGrange <- GRanges(seqnames = ampliconArm$chrom[tmpIdx],
                       IRanges(start = ampliconArm$start.pos[tmpIdx], end = ampliconArm$end.pos[tmpIdx]))
  for (i in unique(ampliconArm$sampleID)) {
    tmpDfArm <- ampliconArm[which(ampliconArm$sampleID == i),]
    tmpDfFocal <- ampliconFocal[which(ampliconFocal$sampleID ==i), ]
    tmpDfFocal <- tmpDfFocal[which(abs(tmpDfFocal$mean) > 0.2), ]
    
    tmpIdx2 <- which(abs(tmpDfArm$mean) < 0.2)
    armSum <- sum(tmpDfArm$length[which(abs(tmpDfArm$mean) > 0.2)])
    
    focalGrange <- GRanges(seqnames = tmpDfFocal$chrom,
                           IRanges(start = tmpDfFocal$start.pos, end = tmpDfFocal$end.pos))
    armGrange <- GRanges(seqnames = tmpDfArm$chrom[tmpIdx2],
                         IRanges(start = tmpDfArm$start.pos[tmpIdx2], end = tmpDfArm$end.pos[tmpIdx2]))
    
    focalSum <- 0
    if (length(tmpIdx2) > 0) {
      for (j in 1:length(armGrange)) {
        tmpSum <- sum(tmpDfFocal$length[subjectHits(findOverlaps(armGrange[j], focalGrange))])
        focalSum <- focalSum + tmpSum
      }
    }
    
    tmpFga <- (armSum + focalSum)/sum(tmpDfArm$length)
    tmpRes <- rbind(tmpRes, data.frame("sample" = i, "fga" = tmpFga))
  }
  return(tmpRes)
}

### gene plot highlights cancer census gens in long plot format
### while also showing the portions of gains/losses

### only plots

# tab <- coad_adenoCar_arm_allSynTable_simple
# geneList <- c("NRAS", "PIK3CA", "SOX2", "FGFR1", "MYC", "BRCA1", "ERBB2", "NF2")
genePlotSynteny <- function(tab, geneList){
  
  ### cancer genes
  tab2 <- tab[, 1:4]
  tab2$h_chr <- str_remove_all(tab$h_chr, "h\\_chr")
  tmpGr <- GRanges(seqnames = tab2$h_chr, IRanges(start = tab2$h_start, end = tab2$h_end))
  
  cancerGeneCensus <- read.table("/mnt/DATA5/tmp/kev/misc/Census_allWed Jul 12 17_36_03 2023.tsv", sep = "\t",
                                 stringsAsFactors = FALSE, header = TRUE)
  tmp <- str_split(cancerGeneCensus$Genome.Location, ":")
  tmpChr <- unlist(lapply(tmp, '[[', 1))
  tmpRange <- unlist(lapply(tmp, '[[', 2))
  removeIdx <- which(tmpRange == "-")
  tmpChr <- tmpChr[-removeIdx]
  tmpRange <- tmpRange[-removeIdx]
  tmpRange2 <- str_split(tmpRange, "-")
  tmpStart <- unlist(lapply(tmpRange2, '[[', 1))
  tmpEnd <- unlist(lapply(tmpRange2, '[[', 2))
  
  cancerGeneCensus  <- cancerGeneCensus[-removeIdx, ]
  
  cancerGeneCensus$chr <- tmpChr
  cancerGeneCensus$start <- tmpStart
  cancerGeneCensus$end <- tmpEnd  
  cancerGeneCensusGr <- GRanges(seqnames = cancerGeneCensus$chr,
                                IRanges(start = as.numeric(cancerGeneCensus$start), end = as.numeric(cancerGeneCensus$end)))
  
  tmpOverlap <- findOverlaps(tmpGr,cancerGeneCensusGr)
  tmpCensusDf <- cancerGeneCensus[subjectHits(tmpOverlap), c("Gene.Symbol", "chr", "start", "end")]
  tmpCensusDf$ySign <- as.numeric(sign(tab2$h_freq[queryHits(tmpOverlap)])) * 0.7
  tmpCensusDf$start <- (as.numeric(tmpCensusDf$start) - 10000)/1e6
  tmpCensusDf$end <- (as.numeric(tmpCensusDf$end) + 10000)/1e6
  colnames(tmpCensusDf)[2:4] <- c("Chr", "Start", "End")
  for (i in unique(tmpCensusDf$Chr)) {
    tmpCensusDf$Start[which(tmpCensusDf$Chr == i)] <- tmpCensusDf$Start[which(tmpCensusDf$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    tmpCensusDf$End[which(tmpCensusDf$Chr == i)] <- tmpCensusDf$End[which(tmpCensusDf$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  if (any(duplicated(tmpCensusDf$Gene.Symbol))) {
    tmpCensusDf <- tmpCensusDf[-which(duplicated(tmpCensusDf$Gene.Symbol)), ]
  }
  
  ### processing the freq data
  chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  dummyPoints <- read.table("/mnt/DATA5/tmp/kev/misc/20210803hg19_dummyPoints.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  chromBreakNoX <- chromBreak[1:23]
  chromTextdfNoX <- chromTextdf[1:22,]
  
  tab3 <- tab2
  colnames(tab3) <- c("Chr", "Start", "End", "Freq")
  tab3$Freq <- tab3$Freq/100
  tab3$Start <- tab3$Start/1e6
  tab3$End <- tab3$End/1e6
  for (i in unique(tab3$Chr)) {
    tab3$Start[which(tab3$Chr == i)] <- tab3$Start[which(tab3$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    tab3$End[which(tab3$Chr == i)] <- tab3$End[which(tab3$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  ### creating gene text to float about marker lines for gene list
  geneLabels <- geneList
  tmpXStart <- tmpCensusDf$Start[which(tmpCensusDf$Gene.Symbol %in% geneList)]
  tmpXEnd <- tmpCensusDf$End[which(tmpCensusDf$Gene.Symbol %in% geneList)]
  geneLabels_x <- apply(cbind(tmpXStart, tmpXEnd),
                        1, mean)
  geneLabels_y <- tmpCensusDf$ySign[which(tmpCensusDf$Gene.Symbol %in% geneList)]
  geneLabels_y <- ifelse(geneLabels_y > 0, geneLabels_y + 0.05, geneLabels_y - 0.05)
  geneLabelsDf <- data.frame("genes" = geneList, "x" = geneLabels_x, "y" = geneLabels_y)
  
  
  ggplot(tab3) + geom_vline(xintercept=chromBreakNoX) + geom_hline(yintercept = 0) + 
    geom_rect(aes(xmin = Start, xmax = End, ymin = 0, ymax=ifelse(Freq>0, Freq, 0)), fill="red") +
    geom_rect(aes(xmin = Start, xmax = End, ymax = 0, ymin=ifelse(Freq<0, Freq,0)), fill="blue") + theme_bw() +
    ylim(c(-1,1)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    geom_text(data = chromTextdfNoX, aes(x = xpos, y = ypos, label = chrom), size = 2.5) +
  # second set of layering is for the lines for genes
    geom_rect(data = tmpCensusDf[which(tmpCensusDf$ySign > 0), ], inherit.aes = FALSE, 
              aes(xmin = Start, xmax = End, ymin = 0, ymax = ySign)) + 
    geom_rect(data = tmpCensusDf[which(tmpCensusDf$ySign < 0), ], inherit.aes = FALSE, 
              aes(xmin = Start, xmax = End, ymax = 0, ymin= ySign)) + 
    geom_text(data =  geneLabelsDf, aes(x = geneLabels_x, y = geneLabels_y, label = geneLabels),
              inherit.aes = FALSE, size = 1.5) + ylab("cn frequency")
}


### made to graph long plots but intuitively show synteny and what is conserved 

# m_amp <- allMouseAneu_bprn_amp_bed
# m_del <- allMouseAneu_bprn_del_bed
# tcga_amp <- armGisticTp53_amp_bed
# tcga_del <- armGisticTp53_del_bed

syntenyOverlap <- function(m_amp, m_del, tcga_amp, tcga_del){
  ### based on the synteny plot
  ### made to look at which arm aneuploidy changes 
  
  chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  dummyPoints <- read.table("/mnt/DATA5/tmp/kev/misc/20210803hg19_dummyPoints.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  chromBreakNoX <- chromBreak[1:23]
  chromTextdfNoX <- chromTextdf[1:22,]
  
  # 20210802 weird upstream bug of 1bp positions with start > end
  if (length(which(m_amp$Start - m_amp$End > 0)) > 0) {
    m_amp <- m_amp[-which(m_amp$Start - m_amp$End > 0),]
    print(1)
  }
  
  if (length(which(m_del$Start - m_del$End > 0)) > 0) {
    m_del <- m_del[-which(m_del$Start - m_del$End > 0),]
    print(2)
  }
  
  if (length(which(tcga_amp$Start - tcga_amp$End > 0)) > 0) {
    tcga_amp <- tcga_amp[-which(tcga_amp$Start - tcga_amp$End > 0),]
    print(3)
  }
  
  if (length(which(tcga_del$Start - tcga_del$End > 0)) > 0) {
    tcga_del <- tcga_del[-which(tcga_del$Start - tcga_del$End > 0),]
    print(4)
  }
  
  tcga_del$Freq  <- tcga_del$Freq * -1
  m_del$Freq <- m_del$Freq * -1
  
  mouse_bed <- NULL
  human_bed <- NULL
  
  if (nrow(m_amp) > 0) {
    ampSynteny <- syntenyPlotInputsFreqV4(tcga_amp, m_amp, species = "human")
    mouse_bed <- rbind(mouse_bed, ampSynteny$mouse_bed)
    human_bed <- rbind(human_bed, ampSynteny$human_bed)
  }
  
  if(nrow(m_del) >  0){
    delSynteny <- syntenyPlotInputsFreqV4(tcga_del, m_del, species = "human")
    mouse_bed <- rbind(mouse_bed, delSynteny$mouse_bed)
    human_bed <- rbind(human_bed, delSynteny$human_bed)
  }
  
  allSynTable <- cbind(human_bed, mouse_bed)
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "h_freq","m_chr",
                             "m_start", "m_end", "m_freq")
  
  allSynTable$synColor <- ifelse(sign(allSynTable$h_freq) == sign(allSynTable$m_freq),
                                 "#013220", "#5A5A5A")
  allSynTable$h_chr2 <- str_remove_all(allSynTable$h_chr, "h\\_chr")
  allSynTable$m_chr2 <- str_remove_all(allSynTable$m_chr, "m\\_chr")
  
  synGraphTable <- allSynTable[, c("h_chr2", "h_start", "h_end", "h_freq",
                                   "m_freq", "synColor")]
  
  synGraphTable$h_freq <- synGraphTable$h_freq/100
  synGraphTable$m_freq <- synGraphTable$m_freq/100
  synGraphTable[, c("h_start", "h_end")] <- synGraphTable[, c("h_start", "h_end")]/1e6
  
  colnames(synGraphTable)[1:4] <- c("Chr", "Start", "End", "Freq")
  
  for (i in unique(synGraphTable$Chr)) {
    synGraphTable$Start[which(synGraphTable$Chr == i)] <- synGraphTable$Start[which(synGraphTable$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    synGraphTable$End[which(synGraphTable$Chr == i)] <- synGraphTable$End[which(synGraphTable$Chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  
  ggplot(synGraphTable) + geom_vline(xintercept=chromBreakNoX) + geom_hline(yintercept = 0) + 
    geom_rect(aes(xmin = Start, xmax = End, ymin = 0,
                  ymax=ifelse(Freq>0, Freq, 0)), fill = synGraphTable$synColor) +
    geom_rect(aes(xmin = Start, xmax = End,
                  ymax = 0, ymin=ifelse(Freq<0, Freq,0)), fill = synGraphTable$synColor) +
    geom_hline(yintercept = 0) + theme_bw() +
    ylim(c(-1,1)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    geom_text(data = chromTextdfNoX, aes(x = xpos, y = ypos, label = chrom), size = 2.5) + ylab("cn frequency")

}


### different then comparing FGA's 
### technically we'll have precalculated freqplot for human
### then compare by sample "synteny" from mouse to human
### end result would be a

# mouseDf <- allMouseAneu_bprn
# m_amp <- allMouseAneu_bprn_amp_bed
# m_del <- allMouseAneu_bprn_del_bed
# tcga_amp <- armGisticTp53_amp_bed
# tcga_del <- armGisticTp53_del_bed


syntenyFraction <- function(m_amp, m_del, tcga_amp, tcga_del, mouseDf, type = "none"){
  ### based on the synteny plot
  
  # 20210802 weird upstream bug of 1bp positions with start > end
  if (length(which(m_amp$Start - m_amp$End > 0)) > 0) {
    m_amp <- m_amp[-which(m_amp$Start - m_amp$End > 0),]
    print(1)
  }
  
  if (length(which(m_del$Start - m_del$End > 0)) > 0) {
    m_del <- m_del[-which(m_del$Start - m_del$End > 0),]
    print(2)
  }
  
  if (length(which(tcga_amp$Start - tcga_amp$End > 0)) > 0) {
    tcga_amp <- tcga_amp[-which(tcga_amp$Start - tcga_amp$End > 0),]
    print(3)
  }
  
  if (length(which(tcga_del$Start - tcga_del$End > 0)) > 0) {
    tcga_del <- tcga_del[-which(tcga_del$Start - tcga_del$End > 0),]
    print(4)
  }
  
  ### for looking at gains and losses separately
  if (type == "gain") {
    m_del$Freq <- 0
    tcga_del$Freq < - 0
  } else if(type == "loss"){
    m_amp$Freq <- 0
    tcga_amp$Freq <- 0
  }
  
  tcga_del$Freq  <- tcga_del$Freq * -1
  m_del$Freq <- m_del$Freq * -1
  
  mouse_bed <- NULL
  human_bed <- NULL
  
  if (nrow(m_amp) > 0) {
    ampSynteny <- syntenyPlotInputsFreqV4(tcga_amp, m_amp, species = "human")
    mouse_bed <- rbind(mouse_bed, ampSynteny$mouse_bed)
    human_bed <- rbind(human_bed, ampSynteny$human_bed)
  }
  
  if(nrow(m_del) >  0){
    delSynteny <- syntenyPlotInputsFreqV4(tcga_del, m_del, species = "human")
    mouse_bed <- rbind(mouse_bed, delSynteny$mouse_bed)
    human_bed <- rbind(human_bed, delSynteny$human_bed)
  }
  
  allSynTable <- cbind(human_bed, mouse_bed)
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "h_freq","m_chr",
                             "m_start", "m_end", "m_freq")
  allSynTable$m_length <- allSynTable$m_end - allSynTable$m_start
  allSynTable$m_chr2 <- str_remove_all(allSynTable$m_chr, "m\\_chr")
  tmpSynGrange <- GRanges(seqnames = allSynTable$m_chr2, IRanges(start = allSynTable$m_start, end = allSynTable$m_end))
  
  tmpRes <- NULL
  # i <- unique(mouseDf$sampleID)[1]
  for (i in unique(mouseDf$sampleID)) {
    tmpMouse <- mouseDf[which(mouseDf$sampleID == i), ]
    tmpMouseGr <- GRanges(seqnames = tmpMouse$chrom, IRanges(start = tmpMouse$start.pos, end = tmpMouse$end.pos))
    tmpHits <- findOverlaps(tmpSynGrange, tmpMouseGr)
    tmpSynComp <- sign(allSynTable$h_freq[queryHits(tmpHits)])
    tmpMouseComp <- sign(tmpMouse$mean[subjectHits(tmpHits)])
    tmpRes <- c(tmpRes, sum(allSynTable$m_length[which(tmpSynComp == tmpMouseComp & tmpSynComp != 0 & tmpMouseComp != 0)])/sum(allSynTable$m_length[which(tmpSynComp != 0 & tmpMouseComp != 0)]))
  }
  return(tmpRes)
}



# humanDf <- tmpHumanDf
# m_amp <- tmpMouseAmp
# m_del <- tmpMouseDel
# tcga_amp <- tmpAmp
# tcga_del <- tmpDel

### based on fraction
syntenyChromosomeScore <- function(m_amp, m_del, tcga_amp, tcga_del, humanDf, species = "human"){
  ### based on the synteny plot
  
  # 20210802 weird upstream bug of 1bp positions with start > end
  if (length(which(m_amp$Start - m_amp$End > 0)) > 0) {
    m_amp <- m_amp[-which(m_amp$Start - m_amp$End > 0),]
    print(1)
  }
  
  if (length(which(m_del$Start - m_del$End > 0)) > 0) {
    m_del <- m_del[-which(m_del$Start - m_del$End > 0),]
    print(2)
  }
  
  if (length(which(tcga_amp$Start - tcga_amp$End > 0)) > 0) {
    tcga_amp <- tcga_amp[-which(tcga_amp$Start - tcga_amp$End > 0),]
    print(3)
  }
  
  if (length(which(tcga_del$Start - tcga_del$End > 0)) > 0) {
    tcga_del <- tcga_del[-which(tcga_del$Start - tcga_del$End > 0),]
    print(4)
  }
  
  tcga_del$Freq  <- tcga_del$Freq * -1
  m_del$Freq <- m_del$Freq * -1
  
  mouse_bed <- NULL
  human_bed <- NULL
  
  if (nrow(m_amp) > 0) {
    ampSynteny <- syntenyPlotInputsFreqV4(tcga_amp, m_amp, species = "human")
    mouse_bed <- rbind(mouse_bed, ampSynteny$mouse_bed)
    human_bed <- rbind(human_bed, ampSynteny$human_bed)
  }
  
  if(nrow(m_del) >  0){
    delSynteny <- syntenyPlotInputsFreqV4(tcga_del, m_del, species = "human")
    mouse_bed <- rbind(mouse_bed, delSynteny$mouse_bed)
    human_bed <- rbind(human_bed, delSynteny$human_bed)
  }
  
  allSynTable <- cbind(human_bed, mouse_bed)
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "h_freq","m_chr",
                             "m_start", "m_end", "m_freq")
  allSynTable$m_length <- allSynTable$m_end - allSynTable$m_start
  allSynTable$m_chr2 <- str_remove_all(allSynTable$m_chr, "m\\_chr")
  allSynTable$h_chr2 <- str_remove_all(allSynTable$h_chr, "h\\_chr")
  
  allMouseChanges <- rbind(m_amp, m_del)
  allMouseChanges$Chr <- paste0("m_chr", allMouseChanges$Chr)
  allMouseChanges <- allMouseChanges[which(allMouseChanges$Freq != 0),]
  
  tmpChrRes <- NULL
  for (i in 1:nrow(allMouseChanges)) {
    allSynTableRed <- allSynTable[which(allSynTable$m_chr == allMouseChanges$Chr[i]), ]
    allSynTableRed$str <- paste(allSynTableRed$h_chr, allSynTableRed$h_start, allSynTableRed$h_end)
    if (any(duplicated(allSynTableRed$str))) {
      allSynTableRed <- allSynTableRed[-which(duplicated(allSynTableRed$str)),]
    }
    tmpSynGrange <- GRanges(seqnames = allSynTableRed$h_chr2, IRanges(start = allSynTableRed$h_start,
                                                                      end = allSynTableRed$h_end))
    tmpSign <- ifelse(sign(allMouseChanges$Freq)[i] == 1, "gain", "loss")
    tmpRes <- NULL
    
    ### meed to put in humanDf not mouse
    
    # j <- unique(humanDf$sampleID)[1]
    for (j in unique(humanDf$sampleID)) {
      
      allSynTableRed2 <- allSynTableRed
      
      tmpHuman <- humanDf[which(humanDf$sampleID == j), ]
      tmpHumanGr <- GRanges(seqnames = tmpHuman$chrom, IRanges(start = tmpHuman$start.pos, end = tmpHuman$end.pos))
      tmpHits <- findOverlaps(tmpSynGrange, tmpHumanGr)
      if (any(duplicated(queryHits(tmpHits)))) {
        tmpQueryDupeIdx <- queryHits(tmpHits)[which(duplicated(queryHits(tmpHits)))]
        allSynTableRed2 <- allSynTableRed2[-tmpQueryDupeIdx, ]
        tmpSynGrange2 <- GRanges(seqnames = allSynTableRed2$h_chr2, IRanges(start = allSynTableRed2$h_start,
                                                                          end = allSynTableRed2$h_end))
        tmpHits2 <- findOverlaps(tmpSynGrange2, tmpHumanGr)
        tmpHumanComp <- sign(tmpHuman$mean[subjectHits(tmpHits2)])
        tmpSynComp <- rep(sign(allMouseChanges$Freq)[i], length(tmpHumanComp))
        
        tmpFractions <- sum(allSynTableRed2$m_length[which(tmpSynComp == tmpHumanComp & tmpSynComp != 0 & tmpHumanComp != 0)])/sum(allSynTableRed2$m_length)
        tmpRes <- c(tmpRes, tmpFractions)
      } else{
        tmpHumanComp <- sign(tmpHuman$mean[subjectHits(tmpHits)])
        tmpSynComp <- rep(sign(allMouseChanges$Freq)[i], length(tmpHumanComp))
        
        tmpFractions <- sum(allSynTableRed2$m_length[which(tmpSynComp == tmpHumanComp & tmpSynComp != 0 & tmpHumanComp != 0)])/sum(allSynTableRed2$m_length)
        tmpRes <- c(tmpRes, tmpFractions)
      }
      ### can't use duplicated regions or regions that span 2 arms b/c not 1:1 on base pair so don't know how to split synteny
    }
    tmpChrRes <- rbind(tmpChrRes, data.frame("m_chr" = rep(paste0(allSynTableRed$m_chr2[i], "_", tmpSign),
                                                           length(tmpRes)), "fraction" = tmpRes))
  }
  return(tmpChrRes)
}
