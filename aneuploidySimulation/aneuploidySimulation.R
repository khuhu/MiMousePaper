### original scripts on brah: (1 & 2) to run parallel simulations: 20250316aneuploidySimulationAbsolute.R, 20250316aneuploidySimulationAbsolute0to50.R
### (3) to process 20250318onlyMidpointAneuSimAbsolute.R



library(GenomicRanges)
library(stringr)
library(parallel)
library(doParallel)


segResAbs <- read.table("/avatar_data5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_segtabs.fixed.txt",
                        header = TRUE, sep = "\t")

absStates <- read.table("/avatar_data5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                        header = TRUE, sep = "\t")


setOfDiploidSamples <- absStates$array[which(absStates$ploidy > 1.6 & absStates$ploidy < 2.3)]

segResAbs_filt <- segResAbs[which(segResAbs$Sample %in% setOfDiploidSamples), ]
segResAbs_filt$cnChage <- segResAbs_filt$Modal_Total_CN - 2

hg19cyto <- read.table("/br_z1/kevin_storage/CCLE/cytoBand.txt", sep = "\t",
                       col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

hg19ArmLocations <- NULL
for (i in unique(hg19cyto$chrom)) {
  tmp <- hg19cyto[which(hg19cyto$chrom == i),]
  tmpP <- tmp[grep("p", tmp$name), ]
  tmpQ <- tmp[grep("q", tmp$name), ]
  hg19ArmLocations <- rbind(hg19ArmLocations, data.frame("chrom" = rep(i, 2), "arm" = c("p", "q"),
                                                         "start" = c(min(tmpP$chromStart), min(tmpQ$chromStart)),                                                        "end" = c(max(tmpP$chromEnd), max(tmpQ$chromEnd))))
}

hg19ArmLocations$chrStripped <- as.numeric(str_remove(hg19ArmLocations$chrom, "chr"))
hg19ArmLocations <- hg19ArmLocations[which(hg19ArmLocations$chrStripped %in% 1:22),]
hg19ArmLocations <- hg19ArmLocations[order(hg19ArmLocations$chrStripped,decreasing = FALSE),]
hg19ArmLocations$str <- paste0(hg19ArmLocations$chrStripped, hg19ArmLocations$arm)
hg19ArmGr <- GRanges(seqnames = hg19ArmLocations$chrom,
                     IRanges(start = hg19ArmLocations$start,
                             end = hg19ArmLocations$end))



### function for detecting arm level alterations based on a certain percentage of the arm being


detectHgArmChangesRegionV2 <- function(df, perc = 0.799, covGr, tableCov){
  ### input is cns call result after calculating cn status
  df <- df[which(df$Chromosome %in% paste0(1:22)), ]
  tableCov$str <- paste0(tableCov$chrom, tableCov$arm)
  tmpGr <- GRanges(seqnames = paste0("chr", df$Chromosome),
                   IRanges(start = df$Start,
                           end = df$End))
  
  tmpGr_subset <- tmpGr[subjectHits(findOverlaps(covGr, tmpGr))]
  tmpDf <- data.frame(GenomicRanges::disjoin(c(covGr, tmpGr_subset)))
  tmpDf_Gr <- GenomicRanges::disjoin(c(covGr, tmpGr_subset))
  toKeep <- subjectHits(findOverlaps(covGr, tmpDf_Gr))
  armMatch <- queryHits(findOverlaps(covGr, tmpDf_Gr))
  tmpDf <- tmpDf[toKeep, ]
  tmpDf$cnChange <- 0
  tmpDf$cnChange[queryHits(findOverlaps(GRanges(tmpDf), tmpGr))] <- df$cnChage[subjectHits(findOverlaps(GRanges(tmpDf), tmpGr))]
  tmpDf$chrArm <- tableCov$str[armMatch]
  
  ### iterate through each of arm and disjoin, so I can simply find coverage percent
  armChangeStatus <- NULL
  percentageCovered <- NULL
  for (i in unique(tmpDf$chrArm)) {
    tmpDf2 <- tmpDf[which(tmpDf$chrArm == i), ]
    
    tmpDelPercent <- signif(sum(tmpDf2$width[which(tmpDf2$cnChange < 0)])/sum(tmpDf2$width), digits = 2)
    tmpAmpPercent <- signif(sum(tmpDf2$width[which(tmpDf2$cnChange > 0)])/sum(tmpDf2$width), digits = 2)
    
    if (isEmpty(tmpDelPercent)) {
      tmpDelPercent <- 0
    }
    
    if (isEmpty(tmpAmpPercent)) {
      tmpAmpPercent <- 0
    }
    
    if (tmpDelPercent > perc) {
      armChangeStatus <- c(armChangeStatus, -1)
      percentageCovered <-  c(percentageCovered, tmpDelPercent)
    } else if(tmpAmpPercent > perc){
      armChangeStatus <- c(armChangeStatus, 1)
      percentageCovered <- c(percentageCovered, tmpAmpPercent)
    } else{
      armChangeStatus <- c(armChangeStatus, 0)
      percentageCovered <- c(percentageCovered, 0)
    }
  }
  
  return(list(armChangeStatus, percentageCovered))
}


### splitting it into chunks b/c its faster


comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

allDiploidNames <- unique(segResAbs_filt$Sample)


### result of calls from Taylor et. al 2018 (Cell)
ataylorTbl <- readxl::read_xlsx("/br_z1/kevin_storage/misc/Taylor2018_cell_supplement1.xlsx", 
                                sheet = 1, skip = 1)



### need to alter the input information and create limits to the ranges i.e what would theoretically be 
### detected by the markers

for (i in seq(0.5, 1, 0.1)) {
  tmpRes <- NULL
  for (j in 1:nrow(hg19ArmLocations)) {
    tmpMiddle <- hg19ArmLocations$start[j] + (hg19ArmLocations$end[j] - hg19ArmLocations$start[j])/2
    tmpIncrement <- (hg19ArmLocations$end[j] - hg19ArmLocations$start[j]) * i / 2
    tmpStart <- tmpMiddle - tmpIncrement
    tmpEnd <- tmpMiddle + tmpIncrement
    tmpRes <- rbind(tmpRes, data.frame("chrom" = hg19ArmLocations$chrom[j],
                                       "arm" = hg19ArmLocations$arm[j],
                                       "start" = tmpStart, "end" = tmpEnd))
  }
  assign(paste0("tableCoverage", 100 * i), tmpRes)
  
  tmpCovGr <- GRanges(seqnames = tmpRes$chrom, 
                      IRanges(start = tmpRes$start, 
                              end = tmpRes$end))
  
  assign(paste0("cov", 100 * i, "Gr"), tmpCovGr)
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov50Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov50Gr,
                                                                         tableCov = tableCoverage50)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov50Res[[1]])
  assign(paste0("cov50_", z * 100, "ArmDf"), tmpDfRes)
  
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov60Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov60Gr,
                                                                         tableCov = tableCoverage60)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov60Res[[1]])
  assign(paste0("cov60_", z * 100, "ArmDf"), tmpDfRes)
  
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov70Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov70Gr,
                                                                         tableCov = tableCoverage70)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov70Res[[1]])
  assign(paste0("cov70_", z * 100, "ArmDf"), tmpDfRes)
  
}



for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov80Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov80Gr,
                                                                         tableCov = tableCoverage80)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov80Res[[1]])
  assign(paste0("cov80_", z * 100, "ArmDf"), tmpDfRes)
  
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov90Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov90Gr,
                                                                         tableCov = tableCoverage90)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov90Res[[1]])
  assign(paste0("cov90_", z * 100, "ArmDf"), tmpDfRes)
  
}



for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov100Res <- foreach(i=0:29,
                                 .combine = 'comb', .multicombine = TRUE,
                                 .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                   
                                   tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                   if (i == 29) {
                                     tmpIterate <- 6381:6591
                                   }
                                   armStatusDf <- NULL
                                   armPercentDf <- NULL
                                   for (j in tmpIterate) {
                                     tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                     tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                          covGr = cov100Gr,
                                                                          tableCov = tableCoverage100)
                                     armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                     armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                   }
                                   rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                   rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                   return(list(armStatusDf, armPercentDf))
                                 }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov100Res[[1]])
  assign(paste0("cov100_", z * 100, "ArmDf"), tmpDfRes)
  
}


### before I do comparison, need to filter out the NA's

ataylorTbl2 <- readxl::read_xlsx("/br_z1/kevin_storage/misc/NIHMS958047-supplement-1.xlsx", 
                                 sheet = 1, skip = 1)
ataylorTbl_red2 <- data.frame(ataylorTbl[,14:52], check.names = FALSE)
rownames(ataylorTbl_red2) <- ataylorTbl2$Sample
ataylorTbl_red2 <- ataylorTbl_red2[which(rownames(ataylorTbl_red2) %in% rownames(cov100_80ArmDf)), ]


tmpColnames <- str_remove(paste0(tableCoverage100$chrom, tableCoverage100$arm), "chr")
cov100_80ArmDf_red <- cov100_80ArmDf[which(rownames(cov100_80ArmDf) %in% rownames(ataylorTbl_red2)), ]
cov100_80ArmDf_red <- cov100_80ArmDf_red[match(rownames(ataylorTbl_red2), rownames(cov100_80ArmDf_red)), ]
cov100_80ArmDf_red <- cov100_80ArmDf_red[, which(tmpColnames %in% colnames(ataylorTbl_red2))]



for (i in seq(0.5, 1, 0.1)) {
  tmpRes <- NULL
  for (j in 1:nrow(hg19ArmLocations)) {
    tmpMiddle <- hg19ArmLocations$start[j] + (hg19ArmLocations$end[j] - hg19ArmLocations$start[j])/2
    tmpIncrement <- (hg19ArmLocations$end[j] - hg19ArmLocations$start[j]) * i / 2
    tmpStart <- tmpMiddle - tmpIncrement
    tmpEnd <- tmpMiddle + tmpIncrement
    tmpRes <- rbind(tmpRes, data.frame("chrom" = hg19ArmLocations$chrom[j],
                                       "arm" = hg19ArmLocations$arm[j],
                                       "start" = tmpStart, "end" = tmpEnd))
  }
  assign(paste0("tableCoverage", 100 * i), tmpRes)
  
  tmpCovGr <- GRanges(seqnames = tmpRes$chrom, 
                      IRanges(start = tmpRes$start, 
                              end = tmpRes$end))
  
  assign(paste0("cov", 100 * i, "Gr"), tmpCovGr)
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov50Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov50Gr,
                                                                         tableCov = tableCoverage50)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov50Res[[1]])
  assign(paste0("cov50_", z * 100, "ArmDf"), tmpDfRes)
  
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov60Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov60Gr,
                                                                         tableCov = tableCoverage60)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov60Res[[1]])
  assign(paste0("cov60_", z * 100, "ArmDf"), tmpDfRes)
  
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov70Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov70Gr,
                                                                         tableCov = tableCoverage70)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov70Res[[1]])
  assign(paste0("cov70_", z * 100, "ArmDf"), tmpDfRes)
  
}



for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov80Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov80Gr,
                                                                         tableCov = tableCoverage80)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov80Res[[1]])
  assign(paste0("cov80_", z * 100, "ArmDf"), tmpDfRes)
  
}


for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov90Res <- foreach(i=0:29,
                                .combine = 'comb', .multicombine = TRUE,
                                .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                  
                                  tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                  if (i == 29) {
                                    tmpIterate <- 6381:6591
                                  }
                                  armStatusDf <- NULL
                                  armPercentDf <- NULL
                                  for (j in tmpIterate) {
                                    tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                    tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                         covGr = cov90Gr,
                                                                         tableCov = tableCoverage90)
                                    armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                    armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                  }
                                  rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                  rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                  return(list(armStatusDf, armPercentDf))
                                }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov90Res[[1]])
  assign(paste0("cov90_", z * 100, "ArmDf"), tmpDfRes)
  
}



for (z in  seq(0.1, 1, 0.1)) {
  cl <- makeCluster(30)
  registerDoParallel(cl)
  start <- Sys.time()
  simulationCov100Res <- foreach(i=0:29,
                                 .combine = 'comb', .multicombine = TRUE,
                                 .init = list(list(), list()), .packages = c("GenomicRanges")) %dopar% {
                                   
                                   tmpIterate <- (i * 220 + 1):((i + 1) * 220)
                                   if (i == 29) {
                                     tmpIterate <- 6381:6591
                                   }
                                   armStatusDf <- NULL
                                   armPercentDf <- NULL
                                   for (j in tmpIterate) {
                                     tmpDf <- segResAbs_filt[which(segResAbs_filt$Sample %in% allDiploidNames[j]), ]
                                     tmpRes <- detectHgArmChangesRegionV2(tmpDf, perc = z,
                                                                          covGr = cov100Gr,
                                                                          tableCov = tableCoverage100)
                                     armStatusDf <- rbind(armStatusDf, tmpRes[[1]])
                                     armPercentDf <- rbind(armPercentDf, tmpRes[[2]])
                                   }
                                   rownames(armStatusDf) <- allDiploidNames[tmpIterate]
                                   rownames(armPercentDf) <- allDiploidNames[tmpIterate]
                                   return(list(armStatusDf, armPercentDf))
                                 }
  
  
  stopCluster(cl)
  print( Sys.time() - start )
  
  tmpDfRes <- do.call(rbind, simulationCov100Res[[1]])
  assign(paste0("cov100_", z * 100, "ArmDf"), tmpDfRes)
  
}


### before I do comparison, need to filter out the NA's

ataylorTbl2 <- readxl::read_xlsx("/br_z1/kevin_storage/misc/Taylor2018_cell_supplement1.xlsx", 
                                 sheet = 1, skip = 1)
ataylorTbl_red2 <- data.frame(ataylorTbl[,14:52], check.names = FALSE)
rownames(ataylorTbl_red2) <- ataylorTbl2$Sample
ataylorTbl_red2 <- ataylorTbl_red2[which(rownames(ataylorTbl_red2) %in% rownames(cov100_80ArmDf)), ]


tmpColnames <- str_remove(paste0(tableCoverage100$chrom, tableCoverage100$arm), "chr")
cov100_80ArmDf_red <- cov100_80ArmDf[which(rownames(cov100_80ArmDf) %in% rownames(ataylorTbl_red2)), ]
cov100_80ArmDf_red <- cov100_80ArmDf_red[match(rownames(ataylorTbl_red2), rownames(cov100_80ArmDf_red)), ]
cov100_80ArmDf_red <- cov100_80ArmDf_red[, which(tmpColnames %in% colnames(ataylorTbl_red2))]
### make another table where I just have aneuploidy instead of gain and loss + no change

confMatRes <- NULL
confMatResV2 <- NULL
i <- 60
j <- 10
for (i in c(60, 70, 80, 90)) {
  print(i)
  for (j in seq(10, 90, 10)) {
    print(j)
    tmpTable <- tryCatch(eval(parse(text = paste0("cov", i, "_", j, "ArmDf"))), error = function(x) return(NULL))
    if (is.null(tmpTable)) {
      next()
    }
    colnames(tmpTable) <- tmpColnames
    tmpTable_red <- tmpTable[which(rownames(tmpTable) %in% rownames(ataylorTbl_red2)), ]
    tmpTable_red <- tmpTable_red[match(rownames(ataylorTbl_red2), rownames(tmpTable_red)), ]
    colnames(tmpTable_red) <- hg19ArmLocations$str
    tmpTable_red <- tmpTable_red[, which(colnames(tmpTable_red) %in% colnames(ataylorTbl_red2))]
    tmpRef2 <- NULL
    tmpCov2 <- NULL
    for (k in 1:ncol(ataylorTbl_red2)) {
      removeIdx <- which(is.na(ataylorTbl_red2[,k]))
      tmpRef <- cov100_80ArmDf_red[-removeIdx, k]
      tmpCov <- tmpTable_red[-removeIdx,k]
      tmpRef2 <- c(tmpRef2, tmpRef)
      tmpCov2 <- c(tmpCov2, tmpCov)
    }
    
    tmpRef3 <- tmpRef2
    tmpCov3 <- tmpCov2
    tmpRef3[tmpRef3 != 0] <- 1
    tmpCov3[tmpCov3 != 0] <- 1
    tmpRef3 <- factor(tmpRef3)
    tmpCov3 <- factor(tmpCov3)
    
    tmpRef2 <- factor(tmpRef2)
    tmpCov2 <- factor(tmpCov2)
    tmpConf <- caret::confusionMatrix(tmpCov2, reference = tmpRef2,
                                      mode = "prec_recall")
    
    
    tmpConf2 <- tmpConf$byClass
    tmpRes <- data.frame("type" = c("loss", "none", "gain"), "genomeCov" = i, "contigCov" = j,
                         "PPV" = tmpConf2[ ,("Pos Pred Value")], "NPV" = tmpConf2[, c("Neg Pred Value")],
                         "Sens" = tmpConf2[ ,("Sensitivity")], "Spec" = tmpConf2[ ,("Specificity")])
    rownames(tmpRes) <- NULL
    confMatRes <- rbind(confMatRes, tmpRes)
    
    ### only one class of anueploidy or no anueploidy
    tmpRef3 <- factor(tmpRef3, levels = c(1, 0))
    tmpCov3 <- factor(tmpCov3, levels = c(1, 0))
    tmpConfV2 <- caret::confusionMatrix(tmpCov3, reference = tmpRef3,
                                        mode = "prec_recall")
    
    tmpConfV2_2 <- tmpConfV2$byClass
    tmpResV2 <- rbind(data.frame("type" = "aneuploidy", "genomeCov" = i, "contigCov" = j,
                                 "PPV" = tmpConfV2_2["Pos Pred Value"], "NPV" = tmpConfV2_2["Neg Pred Value"],
                                 "Sens" = tmpConfV2_2["Sensitivity"], "Spec" = tmpConfV2_2["Specificity"])
                      , tmpRes[2, ])
    
    rownames(tmpResV2) <- NULL
    confMatResV2 <- rbind(confMatResV2, tmpResV2)
    
  }
}

### calculate performance metric for simulations

tableIds <- read.table("/br_z1/kevin_storage/misc/20250318cbioportal_combined_study_clinical_data.tsv",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE,
                       fill = TRUE, quote="")

tableIdsRed <- tableIds[which(tableIds$`Sample ID` %in% rownames(ataylorTbl_red2)), ]
table(tableIdsRed$`TCGA PanCanAtlas Cancer Type Acronym`)

confMatRes <- NULL
confMatResV2byC <- NULL
i <- 10
j <- 10
for (i in seq(10, 100, 10)) {
  print(i)
  for (j in seq(10, 100, 10)) {
    print(j)
    tmpTable <- tryCatch(eval(parse(text = paste0("cov", i, "_", j, "ArmDf"))), error = function(x) return(NULL))
    if (is.null(tmpTable)) {
      next()
    }
    
    colnames(tmpTable) <- tmpColnames
    tmpTable_red <- tmpTable[which(rownames(tmpTable) %in% rownames(ataylorTbl_red2)), ]
    tmpTable_red <- tmpTable_red[match(rownames(ataylorTbl_red2), rownames(tmpTable_red)), ]
    colnames(tmpTable_red) <- hg19ArmLocations$str
    tmpTable_red <- tmpTable_red[, which(colnames(tmpTable_red) %in% colnames(ataylorTbl_red2))]
    z <- unique(tableIdsRed$`TCGA PanCanAtlas Cancer Type Acronym`)[1] 
    for(z in unique(tableIdsRed$`TCGA PanCanAtlas Cancer Type Acronym`)){
      # print(z)
      tmpIds <- tableIdsRed$`Sample ID`[which(tableIdsRed$`TCGA PanCanAtlas Cancer Type Acronym` %in% z)]
      tmpTableCancer <- tmpTable_red[which(rownames(tmpTable_red) %in% tmpIds), ]
      
      tmpTaylorTbl_red2 <- ataylorTbl_red2[which(rownames(ataylorTbl_red2) %in% rownames(tmpTableCancer)), ]
      tmpTaylorTbl_red2 <- tmpTaylorTbl_red2[match(rownames(tmpTableCancer), rownames(tmpTaylorTbl_red2)), ]
      
      cov100_80ArmDf_red2 <- cov100_80ArmDf_red[which(rownames(cov100_80ArmDf_red) %in% rownames(tmpTaylorTbl_red2)), ]
      cov100_80ArmDf_red2 <- cov100_80ArmDf_red2[match(rownames(tmpTaylorTbl_red2), rownames(cov100_80ArmDf_red2)), ]
      
      tmpRef2 <- NULL
      tmpCov2 <- NULL
      for (k in 1:ncol(tmpTaylorTbl_red2)) {
        removeIdx <- which(is.na(tmpTaylorTbl_red2[,k]))
        if(length(removeIdx) < 1){
          tmpRef <- cov100_80ArmDf_red2[ , k]
          tmpCov <- tmpTableCancer[ ,k]
        } else{
          tmpRef <- cov100_80ArmDf_red2[-removeIdx, k]
          tmpCov <- tmpTableCancer[-removeIdx,k]
        }
        tmpRef2 <- c(tmpRef2, tmpRef)
        tmpCov2 <- c(tmpCov2, tmpCov)
      }
      tmpRef3 <- tmpRef2
      tmpCov3 <- tmpCov2
      tmpRef3 <- factor(tmpRef3, levels = c(-1, 0, 1))
      tmpCov3 <- factor(tmpCov3, levels = c(-1, 0, 1))
      

      
      ### only one class of anueploidy or no anueploidy
      tmpRef3 <- factor(tmpRef3)
      tmpCov3 <- factor(tmpCov3)
      tmpConfV2 <- caret::confusionMatrix(tmpCov3, reference = tmpRef3,
                                          mode = "prec_recall")
      
      ppvVar <- (tmpConfV2$table[1,1] + tmpConfV2$table[3,3])/(sum(tmpConfV2$table[1,] + tmpConfV2$table[3,]))
      npvVar <- sum(tmpConfV2$table[2,2])/sum(tmpConfV2$table[2,])
      tprVar <- (tmpConfV2$table[1,1] + tmpConfV2$table[3,3])/(sum(tmpConfV2$table[,1]) + sum(tmpConfV2$table[,3]))
      eventsVar <- sum(tmpConfV2$table[,1]) + sum(tmpConfV2$table[,3])
      tmpConfV2_2 <- tmpConfV2$byClass
      tmpResV2 <- data.frame("type" = "aneuploidy", "genomeCov" = i, "contigCov" = j,
                                   "PPV" = ppvVar, "NPV" = npvVar, "TPR" = tprVar, 
                                   "Cancer" = z, "Events" = eventsVar 
                               )
      
      rownames(tmpResV2) <- NULL
      confMatResV2byC <- rbind(confMatResV2byC, tmpResV2) 
    }
    
  }
}

rownames(confMatResV2byC) <- NULL
confMatRes_filtV2byC <- confMatResV2byC
confMatRes_filtV2byC$confSize <- 1


confMatRes_filtV2byC <- confMatRes_filtV2byC[which(confMatRes_filtV2byC$type == "aneuploidy"), ]
confMatRes_filtV2byC$Cancer2 <- paste(confMatRes_filtV2byC$Cancer,
                                      paste0("n=", confMatRes_filtV2byC$Events))




confMatRes_filtV2byC <- confMatRes_filtV2byC[-which(confMatRes_filtV2byC$genomeCov >= 90 &
                             confMatRes_filtV2byC$contigCov >= 90),]

# pdf("/br_z1/kevin_storage/misc/20250318ppvArmCoverageByCancer.pdf", height = 25, width = 25)
ggplot(confMatRes_filtV2byC, aes(x = genomeCov/100, y = PPV, color = contigCov)) +
  ggtitle("ppv:arm coverage and coverage percentages") + 
  geom_beeswarm() + facet_wrap(~Cancer2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  scale_y_continuous(breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1)) +
  xlab("Percent Arm Coverage") 
dev.off()


# pdf("/br_z1/kevin_storage/misc/20250418ArmCoverageByCancerTpr.pdf", height = 25, width = 25)
ggplot(confMatRes_filtV2byC, aes(x = genomeCov/100, y = TPR, color = contigCov)) +
  ggtitle("ppv:arm coverage and coverage percentages") + 
  geom_beeswarm() + facet_wrap(~Cancer2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  scale_y_continuous(breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1)) +
  xlab("Percent Arm Coverage") 
dev.off()



### doing this for regular ppv together

confMatRes <- NULL
confMatResV2 <- NULL
i <- 50
j <- 50
for (i in seq(10, 100, 10)) {
  print(i)
  for (j in seq(10, 100, 10)) {
    if (i >= 100 & j >= 90) {
      next()
    }
    print(j)
    tmpTable <- tryCatch(eval(parse(text = paste0("cov", i, "_", j, "ArmDf"))),
                         error = function(x) return(NULL))
    if (is.null(tmpTable)) {
      print(paste("i", i))
      print(paste("j", j))
      next()
    }
    colnames(tmpTable) <- tmpColnames
    tmpTable_red <- tmpTable[which(rownames(tmpTable) %in% rownames(ataylorTbl_red2)), ]
    tmpTable_red <- tmpTable_red[match(rownames(ataylorTbl_red2), rownames(tmpTable_red)), ]
    colnames(tmpTable_red) <- hg19ArmLocations$str
    tmpTable_red <- tmpTable_red[, which(colnames(tmpTable_red) %in% colnames(ataylorTbl_red2))]
    tmpRef2 <- NULL
    tmpCov2 <- NULL
    for (k in 1:ncol(ataylorTbl_red2)) {
      removeIdx <- which(is.na(ataylorTbl_red2[,k]))
      if(length(removeIdx) < 1){
        tmpRef <- cov100_80ArmDf_red[ , k]
        tmpCov <- tmpTable_red[ ,k]
      } else{
        tmpRef <- cov100_80ArmDf_red[-removeIdx, k]
        tmpCov <- tmpTable_red[-removeIdx,k]
      }
      ### using a taylor as reference vs our own cov analyis
      
      tmpRef2 <- c(tmpRef2, tmpRef)
      tmpCov2 <- c(tmpCov2, tmpCov)
    }
    
    tmpRef3 <- tmpRef2
    tmpCov3 <- tmpCov2
    tmpRef3 <- factor(tmpRef3)
    tmpCov3 <- factor(tmpCov3)
    
    tmpRef3 <- factor(tmpRef3, levels = c(-1, 0, 1))
    tmpCov3 <- factor(tmpCov3, levels = c(-1, 0, 1))
    tmpConfV2 <- caret::confusionMatrix(tmpCov3, reference = tmpRef3,
                                        mode = "prec_recall")
    
    ppvVar <- (tmpConfV2$table[1,1] + tmpConfV2$table[3,3])/(sum(tmpConfV2$table[1,]) + sum(tmpConfV2$table[3,]))
    npvVar <- sum(tmpConfV2$table[2,2])/sum(tmpConfV2$table[2,])
    tprVar <- (tmpConfV2$table[1,1] + tmpConfV2$table[3,3])/(sum(tmpConfV2$table[,1]) + sum(tmpConfV2$table[,3]))
    
    tmpConfV2_2 <- tmpConfV2$byClass
    tmpResV2 <- data.frame("type" = "aneuploidy", "genomeCov" = i, "contigCov" = j,
                                 "PPV" = ppvVar, "NPV" = npvVar, "TPR" = tprVar)
    
    rownames(tmpResV2) <- NULL
    confMatResV2 <- rbind(confMatResV2, tmpResV2)
    
  }
}


rownames(confMatResV2) <- NULL
confMatRes_filtV2 <- confMatResV2
confMatRes_filtV2$confSize <- 1

confMatRes_filtV2 <- confMatRes_filtV2[which(confMatResV2$type == "aneuploidy"), ]
confMatRes_filtV2[which(confMatRes_filtV2$genomeCov == 90 &
                          confMatRes_filtV2$contigCov == 100), ] <- NA



### plots for visualizing the difference performance metrics


library(ggbeeswarm)

dev.off()
# pdf("/br_z1/kevin_storage/misc/20250318ppvArmCoverage.pdf", height = 10, width = 10)
ggplot(confMatRes_filtV2, aes(x = genomeCov/100, y = PPV, color = contigCov)) +
  ggtitle("ppv:arm coverage and coverage percentages") + 
  geom_beeswarm() +
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  scale_y_continuous(breaks = seq(0.7, 1, 0.1), limits = c(0.7, 1)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1)) +
  xlab("Percent Arm Coverage") 
dev.off()

ggplot(confMatRes_filtV2, aes(x = genomeCov/100, y = TPR, color = contigCov)) +
  ggtitle("ppv:arm coverage and coverage percentages") + 
  geom_beeswarm() +
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  scale_y_continuous(breaks = seq(0.7, 1, 0.1), limits = c(0.7, 1)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1)) +
  xlab("Percent Arm Coverage") 
dev.off()




### this is Rc as an expected threshold
confMatRes_filtV2$Rc <- confMatRes_filtV2$contigCov/100
confMatRes_filtV2$genomeCov2 <- factor(confMatRes_filtV2$genomeCov)

cor(confMatRes_filtV2$Rc[-which(is.nan(confMatRes_filtV2$PPV))],
    confMatRes_filtV2$PPV[-which(is.nan(confMatRes_filtV2$PPV))])

cor(confMatRes_filtV2$Rc[-which(is.na(confMatRes_filtV2$PPV))],
    confMatRes_filtV2$PPV[-which(is.na(confMatRes_filtV2$PPV))])


pdf("/br_z1/kevin_storage/misc/20250318ppvRcCor.pdf", height = 10, width = 10)
ggplot(confMatRes_filtV2, aes(x = Rc, y = PPV)) +
  ggtitle("ppv:arm coverage and coverage percentages") + 
  geom_point() + geom_smooth(method = "lm", color = "red") +
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank())  +
  xlab("Percent Arm Coverage") +
  scale_y_continuous(limits = c(0.75, 1), breaks = seq(0.75, 1, 0.05)) +
  scale_x_continuous(limits = c(0.1, 1), breaks = seq(0.1, 0.9, 0.1))
dev.off()



