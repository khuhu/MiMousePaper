




load("/mnt/DATA4/test_nextflow/20250207crcPaperPortion.RData")

### synteny permutation tests for humans and mice

broadBySampleGisticCoad <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_COADREAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/broad_values_by_arm.txt",
                                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

library(stringr)
library(circlize)
library(doFuture)
library(future)
library(reshape2)
library(GenomicRanges)
library(doParallel)

broadBySampleGisticCoad[25, 4:ncol(broadBySampleGisticCoad)]


comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

### need to fix the empty flag output
coad_adenoCar_arm_allSynTable <- circosFreq2(allMouseAneu_adenoCar_amp_bed, allMouseAneu_adenoCar_del_bed, armGisticCoadApc_amp_bed,
                                             armGisticCoadApc_del_bed, filename = "20231025fearon_aneuAll_coad_adenoCar", ref = "human")


coad_adenoCar_arm_allSynTable$str <- paste(coad_adenoCar_arm_allSynTable$h_chr,
                                           coad_adenoCar_arm_allSynTable$h_start,
                                           coad_adenoCar_arm_allSynTable$h_end)
apcCoadArms <- broadBySampleGisticCoad2
tmpSynTable_adenoCar <- coad_adenoCar_arm_allSynTable[-which(duplicated(coad_adenoCar_arm_allSynTable$str)), ]
tmpSynTable_adenoCar$h_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$h_chr, "h_chr"))
tmpSynTable_adenoCar$m_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$m_chr, "m_chr"))
tmpSynTable_adenoCar <- tmpSynTable_adenoCar[order(tmpSynTable_adenoCar$h_chr, tmpSynTable_adenoCar$h_start, decreasing = FALSE),]
apcCoadArmsGr <- GRanges(seqnames = apcCoadArms$chrStripped, IRanges(start = apcCoadArms$start, end = apcCoadArms$end))
coadSyntenyGrange <- GRanges(seqnames = tmpSynTable_adenoCar$h_chr, IRanges(start = tmpSynTable_adenoCar$h_start, end = tmpSynTable_adenoCar$h_end))
coadMm10SyntenyGrange <- GRanges(seqnames = tmpSynTable_adenoCar$m_chr, IRanges(start = tmpSynTable_adenoCar$m_start, end = tmpSynTable_adenoCar$m_end))




calcSyntenyMatMm <- function(df, tmpSynTable, syntenyGr){
  
  resTable <- data.frame("chr" = tmpSynTable$m_chr, "start" = tmpSynTable$m_start, "end" = tmpSynTable$m_end)
  
  ### given df with sample and chromosome aneuploidy status
  ### return matrix of sample x syntenic block
  for (i in seq_along(unique(df$sampleID))) {
    
    tmp <- df[which(df$sampleID == unique(df$sampleID)[i]),]
    # tmpGrange <- GRanges(seqnames = paste0("m_chr", tmp$chrom), IRanges(start = tmp$start.pos, end = tmp$end.pos))
    tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
    
    tmpIdx <- subjectHits(findOverlaps(query = syntenyGr, subject = tmpGrange))
    testQuery <- queryHits(findOverlaps(query = syntenyGr, subject = tmpGrange))
    
    if (length(which(duplicated(testQuery))) > 0) {
      tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
      testQuery <- testQuery[-which(duplicated( testQuery))]
    }
    ### error probably happens b/c there is no match for one of the queries
    
    resTable[ , i + 3] <- NA
    resTable[testQuery, i + 3] <- tmp$mean[tmpIdx]
  }
  
  resTable2 <- resTable[ , 4:ncol(resTable)]
  return(resTable2)
}

df <- apcCoadArms
tmpSynTable <- tmpSynTable_adenoCar
syntenyGr <- coadSyntenyGrange
syntenyGr2 <- apcCoadArmsGr

calcSyntenyMatHg <- function(df, tmpSynTable, syntenyGr, syntenyGr2){
  
  resTable <- data.frame("chr" = tmpSynTable$h_chr,
                         "start" = tmpSynTable$h_start,
                         "end" = tmpSynTable$h_end)
  
  ### given df with sample and chromosome aneuploidy status
  ### return matrix of sample x syntenic block
  for (i in 4:ncol(df)) {
    
    syntenyGr2Df <- data.frame(syntenyGr2)
    syntenyGr2_1 <- GRanges(seqnames = syntenyGr2Df$seqnames,
                            IRanges(start = syntenyGr2Df$start, 
                                    end = syntenyGr2Df$end))
    tmpIdx <- subjectHits(findOverlaps(query = syntenyGr,
                                       subject = syntenyGr2_1))
    testQuery <- queryHits(findOverlaps(query = syntenyGr,
                                        subject = syntenyGr2_1))
    
    
    if (length(which(duplicated(testQuery))) > 0) {
      tmpIdx <-  tmpIdx[-which(duplicated(testQuery))]
      testQuery <- testQuery[-which(duplicated(testQuery))]
    }
    
    ### error probably happens b/c there is no match for one of the queries
    resTable[ ,i] <- NA
    resTable[ ,i] <- df[tmpIdx, i]
  }
  
  resTable2 <- resTable[ , 4:ncol(resTable)]
  return(resTable2)
}



mouseCoadRegionMat_adenoCar2 <- calcSyntenyMatMm(allMouseAneu_adenoCar, tmpSynTable_adenoCar,
                                                 coadMm10SyntenyGrange)

humanCrcRegionMat2 <- calcSyntenyMatHg(apcCoadArms, tmpSynTable_adenoCar,
                                       coadSyntenyGrange, apcCoadArmsGr)


### creating the test statistic for permutation test

twoSampPropZTest <- function(x, y, pooled = TRUE){
  if(pooled == TRUE){
    n1 <- y[1]
    n2 <- y[2]
    p1 <- x[1]/n1
    p2 <- x[2]/n2
    p0 <- (x[1] + x[2])/(n1 + n2)
    tmpZ <- (p1 - p2)/sqrt(p0 * (1 - p0) * (1/n1 + 1/n2))
  } else{
    n1 <- y[1]
    n2 <- y[2]
    p1 <- x[1]/n1
    p2 <- x[2]/n2
    tmpZ <- (p1 - p2)/sqrt((p1 * (1 - p1))/n1 + (p2 * (1 - p2))/n2)
  }
  return(tmpZ)
}



libs <- .libPaths()

doFuture::registerDoFuture()
future::plan("multisession")

cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synMmCrcPerm <- foreach(i=1:25,
                        .combine = 'comb', .multicombine = TRUE,
                        .init = list(list(), list(), list(), list(), list(), list(), list())) %dopar% {
                          
                          .libPaths(libs)
                          library(reshape2)
                          library(GenomicRanges)
                          
                          tmpAmpAll <- NULL
                          tmpDelAll <- NULL
                          tmpAmpFreq <- NULL
                          tmpDelFreq <- NULL
                          j <- 0
                          
                          tmpAmpStatAll <- NULL
                          tmpDelStatAll <- NULL
                          tmpAmpPAll <- NULL
                          tmpDelPAll <- NULL
                          while (j < 200) {
                            
                            ### resampling across samples for chromosome
                            
                            tmpDcast <- reshape2::dcast(allMouseAneu_adenoCar,
                                                        formula = sampleID ~ chrom,
                                                        value.var = "mean") 
                            tmpDcast2 <- data.frame("sampleID" = tmpDcast$sampleID, 
                                                    t(apply(tmpDcast[, 2:ncol(tmpDcast)], 1,
                                                            function(x) sample(x, (ncol(tmpDcast) - 1)))),
                                                    check.names = FALSE);
                            tmpDcast3 <- melt(tmpDcast2)
                            colnames(tmpDcast3) <- c("sampleID", "chrom", "mean")
                            tmpDcast3$start.pos <- 1
                            tmpMmLocs <- allMouseAneu_adenoCar[-which(duplicated(allMouseAneu_adenoCar$end.pos)), ]
                            tmpDcast3$end.pos <- tmpMmLocs$end.pos[match(tmpDcast3$chrom,tmpMmLocs$chrom)]
                            tmp <- calcSyntenyMatMm(tmpDcast3, tmpSynTable_adenoCar, coadMm10SyntenyGrange)
                            
                            
                            # tmp <- apply(mouseCoadRegionMat_adenoCar2, 2, function(x) sample(x, nrow(mouseCoadRegionMat_adenoCar2)))
                            
                            ampFreq <- apply(tmp, 1, function(x) length(which(x > 0.2))/length(x))
                            delFreq <- apply(tmp, 1, function(x) length(which(x < -0.2))/length(x))
                            tmpAmpFreq <- cbind(tmpAmpFreq, ampFreq)
                            tmpDelFreq <- cbind(tmpDelFreq, delFreq)
                            
                            ampStatRes <- NULL
                            delStatRes <- NULL
                            pValAmpRes <- NULL
                            pValDelRes <- NULL
                            
                            
                            ### calculating frequencies of co-occurrence for mouse arms sytenic to 
                            ### target human arm:
                            
                            checkMatrix <- tmp[c(3, 5, 15, 17, 20, 56), ]
                            count1 <- length(which(apply(checkMatrix, 2, sum) == nrow(checkMatrix)))
                            
                            checkMatrix <- tmp[c(3, 5, 20), ]
                            count2 <- length(which(apply(checkMatrix, 2, sum) == nrow(checkMatrix)))
                            
                            checkMatrix <- tmp[c(4, 5, 9, 16, 53), ]
                            count3 <- length(which(apply(checkMatrix, 2, sum) == nrow(checkMatrix)))
                            
                            checkMatrix <- tmp[32, ]
                            count4 <- length(which(apply(checkMatrix, 2, sum) == nrow(checkMatrix)))
                            
                            tmpComboArmRes <- c(count1, count2, count3, count4)
                            
                            for (q in 1:396) {
                              
                              ampProp <- c(length(which(tmp[q, ] > 0.2)), length(which(tmp[-q, ] > 0.2)))
                              delProp <- c(length(which(tmp[q, ] < -0.2)), length(which(tmp[-q, ] < -0.2)))
                              totalProp <- c(length(unlist(tmp[q, ])), length(unlist(tmp[-q, ])))
                              
                              ampTestStat <- twoSampPropZTest(ampProp, totalProp, pooled = FALSE)
                              delTestStat <- twoSampPropZTest(delProp, totalProp, pooled = FALSE)
                              
                              # prop.test(x = ampProp, n = totalProp, correct = FALSE,
                              #           alternative = "less")
                              
                              ampStatRes <- c(ampStatRes, ampTestStat)
                              delStatRes <- c(delStatRes, delTestStat)
                              pValAmpRes <- c(pValAmpRes, pnorm(ampTestStat, lower.tail = FALSE))
                              pValDelRes <- c(pValDelRes, pnorm(delTestStat, lower.tail = FALSE))
                            }
                            
                            tmpAmpStatAll <- cbind(tmpAmpStatAll, ampStatRes)
                            tmpDelStatAll <- cbind(tmpDelStatAll, delStatRes)
                            tmpAmpPAll <- cbind(tmpAmpPAll, pValAmpRes) 
                            tmpDelPAll <- cbind(tmpDelPAll, pValDelRes)
                            
                            j <- j + 1
                          }
                          
                          return(list(tmpAmpFreq, tmpDelFreq, tmpAmpStatAll, tmpDelStatAll, tmpAmpPAll, tmpDelPAll, tmpComboArmRes))
                        }


stopCluster(cl)
print( Sys.time() - start )





### using the results of from permutaton test for different statistical comparisons


mm10ChrIdx <- which(!duplicated(tmpSynTable_adenoCar$m_chr))
tmpDf <- tmpSynTable_adenoCar
tmpDf$arm <- "q"
mm10ChrPermDf_yes <- tmpDf[mm10ChrIdx, c("m_chr", "arm")]
mm10ChrPermDf_no <- tmpDf[mm10ChrIdx, c("m_chr", "arm")]

### need to make this for chromosome arms in humans
hg19cyto <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/hg19/cytoBand.txt", sep = "\t",
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
hg19ArmLocations <- hg19ArmLocations[which(hg19ArmLocations$str %in% broadBySampleGisticCoad$`Chromosome Arm`), ]


### human matrix is correct, but somehow this mapping is off

hg19ArmLocationsGr <- GRanges(seqnames = hg19ArmLocations$chrom,
                              IRanges(start = hg19ArmLocations$start,
                                      hg19ArmLocations$end))

syntenyGr <- GRanges(seqnames = paste0("chr",tmpSynTable_adenoCar$h_chr),
                     IRanges(start = tmpSynTable_adenoCar$h_start,
                             end = tmpSynTable_adenoCar$h_end))

### for some odd reason, this isn't doing what I want it to do
# tmpIdx <- which(!duplicated(queryHits(findOverlaps(hg19ArmLocationsGr, syntenyGr))))
# hg19ChrIdx <- subjectHits(findOverlaps(hg19ArmLocationsGr, syntenyGr))[tmpIdx]
tmpIdx <- NULL
for (i in 1:length(hg19ArmLocationsGr)) {
  res <- subjectHits(findOverlaps(hg19ArmLocationsGr[i], syntenyGr))[1]
  tmpIdx <- c(tmpIdx, res)
}
hg19ChrIdx <- tmpIdx
hg19ChrPermDf_yes <- hg19ArmLocations[c("chrom", "arm")]
hg19ChrPermDf_no <- hg19ArmLocations[c("chrom", "arm")]


hg19_yesCount <- NULL
hg19_noCount <- NULL
i <- 148
for (i in hg19ChrIdx) {
  
  
  ### directionality only works b/c we're only look at gains
  yes_hgAlt <- length(which(humanCrcRegionMat2[i, ] > 0.2))
  no_hgAlt <- length(which(humanCrcRegionMat2[i, ] < 0.2))
  
  hg19_yesCount  <- c(hg19_yesCount, yes_hgAlt)
  hg19_noCount <- c(hg19_noCount, no_hgAlt)
  
}

hg19ChrPermDf_yes$refCount <- hg19_yesCount
hg19ChrPermDf_no$refCount <- hg19_noCount

### now I need to know the sytenic regions between human to mouse chrs for
### alterations of interest
hg19ArmLocations_red <- hg19ArmLocations[which(hg19ArmLocations$str %in%
                                                 c("7p", "7q", "13q", "20q")), ]
hg19ArmLocations_redGr <- GRanges(seqnames = hg19ArmLocations_red$chrom,
                                  IRanges(start = hg19ArmLocations_red$start,
                                          end = hg19ArmLocations_red$end))
tmpGr <- GRanges(seqnames = tmpSynTable_adenoCar$h_chr,
                 IRanges(start = tmpSynTable_adenoCar$h_start,
                         end = tmpSynTable_adenoCar$h_end))

i <- 2
hg19ArmLocations_red$syntenic <- NA
for (i in 1:length(hg19ArmLocations_redGr)) {
  tmpIdx <- subjectHits(findOverlaps(hg19ArmLocations_redGr[i], syntenyGr))
  hg19ArmLocations_red$syntenic[i] <- paste(unique(tmpSynTable_adenoCar$m_chr[tmpIdx]),
                                            collapse = " ")
}

### calculating regulat chi-squared test for one vs all
### then just one vs at least one
i <- 1

### make sure the proportions look right - before I used the wrong
### synteny chart to start it

### got my ways of calculating it; now to load in the the 10,000 permutations and calculate the p-values
### starting making table here, append it with later code, one table per arm

allSyntenicArms <- NULL
allSyntenicArmsOr <- NULL
allSyntenicArmsX <- NULL
for (i in 1:nrow(hg19ArmLocations_red)) {
  
  hgIdx <- which(paste(hg19ChrPermDf_yes$chrom, hg19ChrPermDf_yes$arm) == 
                   paste(hg19ArmLocations_red$chrom[i], hg19ArmLocations_red$arm[i]))
  
  syntenicChrs <- unlist(strsplit(hg19ArmLocations_red$syntenic[i], " "))
  
  tmpMmIdx <- which(tmpSynTable_adenoCar$m_chr %in% syntenicChrs)
  tmpMmIdx2 <- tmpMmIdx[which(!duplicated(tmpSynTable_adenoCar$m_chr[tmpMmIdx]))]
  
  
  
  checkMatrix <- mouseCoadRegionMat_adenoCar2[tmpMmIdx2, ]
  colnames(checkMatrix) <- NULL
  yes_mmAlt <- length(which(apply(checkMatrix, 2, sum) == nrow(checkMatrix)))
  no_mmAlt <- length(which(apply(checkMatrix, 2, sum) != nrow(checkMatrix)))
  
  yes_hg_alt <- hg19ChrPermDf_yes$refCount[hgIdx]
  no_hg_alt <- hg19ChrPermDf_no$refCount[hgIdx]
  
  
  tmpConting <- t(matrix(c(yes_mmAlt, no_mmAlt, yes_hg_alt, no_hg_alt),
                         nrow=2,ncol=2, 
                         dimnames=list(c("aneuploidy","no_aneuploidy"),
                                       c("mm10","hg19"))))

  chiRes <- chisq.test(tmpConting)
  
  allSyntenicArms <- rbind(allSyntenicArms, 
                           data.frame("arm" = hg19ArmLocations_red$str[i],
                                      "pval" = signif(chiRes$p.value, digits = 3)))
  
  
  # tmpOddsRatio <- log(tmpConting[1,1] * tmpConting[2,2]/(tmpConting[1,2] * tmpConting[2,1]))
  tmpNa <- tmpConting[1,1] + tmpConting[1,2]
  tmpNb <- tmpConting[2,1] + tmpConting[2,2]
  # tmpOddsRatio <- log((((tmpConting[1,1]+0.5)/tmpNa) * ((tmpConting[2,2] + 0.5)/tmpNa))/(((tmpConting[1,2] + 0.5)/tmpNb) * ((tmpConting[2,1] + 0.5)/tmpNb)))
  tmpOddsRatio <- log((((tmpConting[1,1]+0.5)) * ((tmpConting[2,2] + 0.5)))/(((tmpConting[1,2] + 0.5)) * ((tmpConting[2,1] + 0.5))))
  
  allSyntenicArmsOr <- rbind(allSyntenicArmsOr,
                             data.frame("arm" = hg19ArmLocations_red$str[i],
                                        "or" = signif( tmpOddsRatio , digits = 3)))
  
  
  allSyntenicArmsX <- rbind(allSyntenicArmsX,
                             data.frame("arm" = hg19ArmLocations_red$str[i],
                                        "chi" = signif(chiRes$statistic, digits = 3)))
  
}


### code for singular genes



singularSyntenicArms <- NULL
singularSyntenicArmsOr <- NULL
singularSyntenicArmsX <- NULL
for (i in 1:nrow(hg19ArmLocations_red)) {
  
  hgIdx <- which(paste(hg19ChrPermDf_yes$chrom, hg19ChrPermDf_yes$arm) == 
                   paste(hg19ArmLocations_red$chrom[i], hg19ArmLocations_red$arm[i]))
  
  syntenicChrs <- unlist(strsplit(hg19ArmLocations_red$syntenic[i], " "))
  
  tmpMmIdx <- which(tmpSynTable_adenoCar$m_chr %in% syntenicChrs)
  tmpMmIdx2 <- tmpMmIdx[which(!duplicated(tmpSynTable_adenoCar$m_chr[tmpMmIdx]))]
  
  
  ### for permutation to return 
  print(tmpMmIdx2)
  
  ### just need to do check differently
  checkMatrix <- mouseCoadRegionMat_adenoCar2[tmpMmIdx2, ]
  colnames(checkMatrix) <- NULL
  checkMatrix <- matrix(as.numeric(unlist(checkMatrix )),nrow=nrow(checkMatrix ))
  
  
  yes_hg_alt <- hg19ChrPermDf_yes$refCount[hgIdx]
  no_hg_alt <- hg19ChrPermDf_no$refCount[hgIdx]
  for(j in 1:nrow(checkMatrix)){
    yes_mmAlt <- length(which(checkMatrix[j,] == 1))
    no_mmAlt <- length(which(checkMatrix[j,] != 1))
    
    tmpConting <- t(matrix(c(yes_mmAlt, no_mmAlt, yes_hg_alt, no_hg_alt),
                           nrow=2,ncol=2, 
                           dimnames=list(c("aneuploidy","no_aneuploidy"),
                                         c("mm10","hg19"))))
    # 
    chiRes <- chisq.test(tmpConting)

    
    singularSyntenicArms <- rbind(singularSyntenicArms,
                                  data.frame("hg19_arm" = hg19ArmLocations_red$str[i],
                                             "mm10_arm" = tmpSynTable_adenoCar$m_chr[tmpMmIdx2[j]],
                                             "pval" = chiRes$p.value))
    
    
    tmpNa <- tmpConting[1,1] + tmpConting[1,2]
    tmpNb <- tmpConting[2,1] + tmpConting[2,2]
    tmpOddsRatio <- log((((tmpConting[1,1]+0.5)) * ((tmpConting[2,2] + 0.5)))/(((tmpConting[1,2] + 0.5)) * ((tmpConting[2,1] + 0.5))))
    
    singularSyntenicArmsOr <- rbind(singularSyntenicArmsOr,
                                    data.frame("hg19_arm" = hg19ArmLocations_red$str[i],
                                               "mm10_arm" = tmpSynTable_adenoCar$m_chr[tmpMmIdx2[j]],
                                               "or" = tmpOddsRatio ))

    singularSyntenicArmsX <- rbind(singularSyntenicArmsX,
                              data.frame("hg19_arm" = hg19ArmLocations_red$str[i],
                                         "mm10_arm" = tmpSynTable_adenoCar$m_chr[tmpMmIdx2[j]],
                                         "chi" = signif(chiRes$statistic, digits = 3)))
  }
}


load("/mnt/DATA4/test_nextflow/20250214synMmCrcPerm.RData")
load("/mnt/DATA4/test_nextflow/20250214synHgCrcPerm.RData")

mousePermAmpFreq_Crc <- do.call(cbind, synMmCrcPerm[[1]])
mousePermDelFreq_Crc <- do.call(cbind, synMmCrcPerm[[2]])

humanPermAmpCrcFreq <- do.call(cbind, synCrcPerm[[1]])
humanPermDelCrcFreq <- do.call(cbind, synCrcPerm[[2]])

mousePermAmpFreq_CrcRed <- mousePermAmpFreq_Crc[, ]


### running these permutations through sample test and calculating p-value

comboPvalMatrix <- NULL
comboOrMatrix <- NULL
comboResidMatrix <- NULL

for (i in 1:nrow(hg19ArmLocations_red)) {
  
  hgIdx <- which(paste(hg19ChrPermDf_yes$chrom, hg19ChrPermDf_yes$arm) == 
                   paste(hg19ArmLocations_red$chrom[i], hg19ArmLocations_red$arm[i]))
  
  tmpHgIdx <- which(tmpSynTable_adenoCar$h_chr %in% str_remove(hg19ArmLocations_red$chrom[i], "chr"))
  tmpHgIdx2 <- tmpHgIdx[which(!duplicated(tmpSynTable_adenoCar$h_chr[tmpHgIdx]))]
  
  syntenicChrs <- unlist(strsplit(hg19ArmLocations_red$syntenic[i], " "))
  
  tmpMmIdx <- which(tmpSynTable_adenoCar$m_chr %in% syntenicChrs)
  tmpMmIdx2 <- tmpMmIdx[which(!duplicated(tmpSynTable_adenoCar$m_chr[tmpMmIdx]))]
  
  checkMatrix <- mousePermAmpFreq_Crc[tmpMmIdx2, ]
  colnames(checkMatrix) <- NULL
  
  tmpArmPvals <- NULL
  tmpArmOr <- NULL
  tmpArmResid <- NULL
  ### like above said, for now set it to 0
  for (j in 1:5000) {
    yes_mmAlt <- 0
    no_mmAlt <- 28
    
    yes_hg_alt <- round(humanPermAmpCrcFreq[tmpHgIdx2, j] * 619) # frequency times total
    no_hg_alt <- pmin((619 - no_hg_alt), 619)
    
    tmpConting <- t(matrix(c(yes_mmAlt, no_mmAlt, yes_hg_alt, no_hg_alt),
                           nrow=2,ncol=2, 
                           dimnames=list(c("aneuploidy","no_aneuploidy"),
                                         c("mm10","hg19"))))
    f.test <- fisher.test(tmpConting)
    chiRes <- chisq.test(tmpConting)
    tmpOddsRatio <- (tmpConting[1,1] * tmpConting[2,2])/(tmpConting[1,2] * tmpConting[2,1])
    tmpArmPvals <- c(tmpArmPvals, chiRes$p.value)
    tmpArmOr <- c(tmpArmOr, tmpOddsRatio)
    tmpArmOr <- c(tmpArmOr, chiRes$statistic)

  }
  
  comboPvalMatrix  <- rbind(comboPvalMatrix , tmpArmPvals)
  comboOrMatrix <- rbind(comboOrMatrix, tmpArmOr)
}

allSyntenicArms2 <- cbind(allSyntenicArms, comboPvalMatrix)
allSyntenicArmsOr2 <- cbind(allSyntenicArmsOr, comboOrMatrix)

### calculatin permutaiton pval
for (i in 1:nrow(allSyntenicArms2)) {
  ### find pvalues less than or equal too
  permPval <- (length(which(allSyntenicArms2[i, 3:ncol(allSyntenicArms2)] <= allSyntenicArms2[i, 2])))/5000
  print(permPval)
}

for (i in 1:nrow(allSyntenicArmsOr2)) {
  ### find pvalues less than or equal too
  permPval <- (length(which(allSyntenicArmsOr2[i, 3:ncol(allSyntenicArmsOr2)] >= allSyntenicArmsOr2[i, 2])))/5000
  print(permPval)
}

### for singular comps
###
###

singularPvalMatrix <- NULL
singularOrMatrix <- NULL
for (i in 1:nrow(hg19ArmLocations_red)) {
  
  hgIdx <- which(paste(hg19ChrPermDf_yes$chrom, hg19ChrPermDf_yes$arm) == 
                   paste(hg19ArmLocations_red$chrom[i], hg19ArmLocations_red$arm[i]))
  
  tmpHgIdx <- which(tmpSynTable_adenoCar$h_chr %in% str_remove(hg19ArmLocations_red$chrom[i], "chr"))
  tmpHgIdx2 <- tmpHgIdx[which(!duplicated(tmpSynTable_adenoCar$h_chr[tmpHgIdx]))]
  
  syntenicChrs <- unlist(strsplit(hg19ArmLocations_red$syntenic[i], " "))
  
  tmpMmIdx <- which(tmpSynTable_adenoCar$m_chr %in% syntenicChrs)
  tmpMmIdx2 <- tmpMmIdx[which(!duplicated(tmpSynTable_adenoCar$m_chr[tmpMmIdx]))]
  
  if (length(tmpMmIdx2) == 1) {
    checkMatrix <- matrix(mousePermAmpFreq_Crc[tmpMmIdx2, ],nrow = 1)
    
  } else{
    checkMatrix <- mousePermAmpFreq_Crc[tmpMmIdx2, ]
  }
  
  colnames(checkMatrix) <- NULL
  
  ### like above said, for now set it to 0
  for(k in 1:nrow(checkMatrix)){
    
    tmpArmPvals <- NULL
    tmpArmOr <- NULL
    for (j in 1:5000) {
      yes_mmAlt <- round(checkMatrix[k, j] * 28)
      no_mmAlt <- 28 - yes_mmAlt
      
      yes_hg_alt <- round(humanPermAmpCrcFreq[tmpHgIdx2, j] * 619) # frequency times total
      no_hg_alt <- pmin((619 - no_hg_alt), 619)
      
      tmpConting <- t(matrix(c(yes_mmAlt, no_mmAlt, yes_hg_alt, no_hg_alt),
                             nrow=2,ncol=2, 
                             dimnames=list(c("aneuploidy","no_aneuploidy"),
                                           c("mm10","hg19"))))
      
      
      f.test <- fisher.test(tmpConting)
      chiRes <- chisq.test(tmpConting)
      # tmpArmPvals <- c(tmpArmPvals, f.test$p.value)
      tmpArmPvals <- c(tmpArmPvals, chiRes$p.value)
      tmpArmOr <- c(tmpArmOr, chiRes$statistic)
    }
    singularPvalMatrix  <- rbind(singularPvalMatrix , tmpArmPvals)
    singularOrMatrix <- rbind(singularOrMatrix, tmpArmOr)
  }
}

singularSyntenicArms2 <- cbind(singularSyntenicArms, singularPvalMatrix)

### calculatin permutation pval
for (i in 1:nrow(singularSyntenicArms2)) {
  ### find pvalues less than or equal too
  permPval <- (5000 - length(which(singularSyntenicArms2[i, 4:ncol(singularSyntenicArms2)] > singularSyntenicArms2[i, 3])))/5000
  print(paste(singularSyntenicArms2$hg19_arm[i],
              singularSyntenicArms2$mm10_arm[i], permPval))
}


tmp <- broadBySampleGisticCoad2
rownames(tmp) <- NULL
length(which(tmp[13, 4:ncol(tmp)] > 0.2 & tmp[14, 4:ncol(tmp)] > 0.2))/(ncol(broadBySampleGisticCoad2) - 3)
length(which(tmp[13, 4:ncol(tmp)] > 0.2 & tmp[25, 4:ncol(tmp)] > 0.2))/(ncol(broadBySampleGisticCoad2) - 3)
length(which(tmp[13, 4:ncol(tmp)] > 0.2 & tmp[14, 4:ncol(tmp)] > 0.2 & tmp[25, 4:ncol(tmp)] > 0.2))/(ncol(broadBySampleGisticCoad2) - 3)

### having singular human arm to syntenic mouse arms works best here
### creating function to make specific circos plots

circosFreq2 <- function(m_amp, m_del, tcga_amp, tcga_del, filename = "test",
                        geneVar = "no", ref = "human", empty = FALSE, reduced = FALSE,
                        show_lower = TRUE){
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
  
  
  tcga_del$Freq  <- tcga_del$Freq
  m_del$Freq <- m_del$Freq
  # tcga_del$Freq  <- tcga_del$Freq * -1
  # m_del$Freq <- m_del$Freq * -1
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
  
  
  if (reduced == TRUE) {
    tmpMouse <- unique(mouse_bed$chr)
    tmpMouse <- tmpMouse[order(as.numeric(str_remove(tmpMouse, "m_chr")))]
    tmpHuman <- unique(human_bed$chr)
    tmpHuman <- tmpHuman[order(as.numeric(str_remove(tmpHuman, "h_chr")))]
    cyto_interest <- c(tmpMouse, tmpHuman)
    cyto_combined_red <- cyto_hg19mm10[which(cyto_hg19mm10$V1 %in% cyto_interest),]
  } else{
    cyto_interest <- c(paste0("m_chr", 1:19), paste0("h_chr", 1:22))
    cyto_combined_red <- cyto_hg19mm10[which(cyto_hg19mm10$V1 %in% cyto_interest),]
  }
  
  
  
  mouse_cn <- trackBedColumns(rbind(m_amp,
                                    m_del))
  mouse_cn$chr <- paste0("m_chr", mouse_cn$chr)
  
  human_cn <- trackBedColumns(rbind(tcga_amp,
                                    tcga_del))
  human_cn$chr <- paste0("h_chr", human_cn$chr)
  
  if (show_lower == TRUE) {
    human_cn$ytop[which(human_cn$freq != 0 )]  <- human_cn$ytop[which(human_cn$freq != 0 )] + 0.05
    mouse_cn$ytop[which(mouse_cn$freq != 0 )]  <- mouse_cn$ytop[which(mouse_cn$freq != 0 )] + 0.05
  }
  
  cn_track_bed <- rbind(mouse_cn, human_cn)
  cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]
  
  track_color <- rep("#000000", nrow(cn_track_bed))
  track_color <- ifelse(cn_track_bed$freq < 0 , "#00008B", "#8B0000") 
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(mouse_bed$chr)
  human_chrs <- unique(human_bed$chr)
  
  if (ref == "human") {
    colorVector2 <- rep("#000000", nrow(human_bed))
    for (j in seq_along(human_chrs)) {
      colorVector2[which(human_bed$chr == human_chrs[j])] <- colorVector[j]
    }
  } else if(ref == "mouse"){
    colorVector2 <- rep("#000000", nrow(mouse_bed))
    for (j in seq_along(mouse_chrs)) {
      colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
    }
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
      ### lowered the ceilings as frequencies never reach those and the
      ### values look smaller than they actually are
                   }) +
      circos.genomicTrack(cn_track_bed, ylim = c(-0.85, 0.85),
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
}

crcSynTable <- circosFreq2(allMouseAneu_adenoCar_amp_bed,
                                              allMouseAneu_adenoCar_del_bed,
                                              armGisticCoadApc_amp_bed, 
                                              armGisticCoadApc_del_bed,
                                              filename = "20250214crc", ref = "human", reduced = TRUE)



crcRedMm <- crcSynTable[, c("m_chr", "m_start", "m_end", "m_freq")]
crcRedHg <- crcSynTable[, c("h_chr", "h_start", "h_end", "h_freq")]
crcRedMm$m_chr <- str_remove(crcRedMm$m_chr, "m_chr")
crcRedHg$h_chr <- str_remove(crcRedHg$h_chr, "h_chr")
colnames(crcRedHg) <- c("Chr", "Start", "End", "Freq")
colnames(crcRedMm) <- c("Chr", "Start", "End", "Freq")


crcRedHg_h20 <- crcRedHg[which(crcRedHg$Chr == "20"), ]
crcRedMm_h20 <- crcRedMm[which(crcRedHg$Chr == "20"), ]

crcRedHg_h20_amp <- crcRedHg_h20[which(crcRedHg_h20$Freq > 0),]
crcRedHg_h20_del <- crcRedHg_h20[which(crcRedHg_h20$Freq < 0),]
crcRedHg_h20_del$Freq <- 0

crcRedMm_h20_amp <- crcRedMm_h20[1:(nrow(crcRedMm_h20)/2), ]
crcRedMm_h20_del <- crcRedMm_h20[(nrow(crcRedMm_h20)/2):nrow(crcRedMm_h20), ]
crcRedMm_h20_del$Freq <- 0


crcExampleSingleH20 <- circosFreq2(crcRedMm_h20_amp,
                                            crcRedMm_h20_del,
                                            crcRedHg_h20_amp,
                                            crcRedHg_h20_del,
                                            filename = "20250214syntenyH20",
                                            ref = "human", reduced = TRUE)


crcRedHg_h13 <- crcRedHg[which(crcRedHg$Chr == "13"), ]
crcRedMm_h13 <- crcRedMm[which(crcRedHg$Chr == "13"), ]

crcRedHg_h13_amp <- crcRedHg_h13[which(crcRedHg_h13$Freq > 0),]
crcRedHg_h13_del <- crcRedHg_h13[which(crcRedHg_h13$Freq > 0),]
crcRedHg_h13_del$Freq <- 0

crcRedMm_h13_amp <- crcRedMm_h13[1:(nrow(crcRedMm_h13)/2), ]
crcRedMm_h13_del <- crcRedMm_h13[(nrow(crcRedMm_h13)/2):nrow(crcRedMm_h13), ]
crcRedMm_h13_del$Freq <- 0


crcExampleSingleH13 <- circosFreq2(crcRedMm_h13_amp,
                                   crcRedMm_h13_del,
                                   crcRedHg_h13_amp,
                                   crcRedHg_h13_del,
                                   filename = "20250214syntenyH13",
                                   ref = "human", reduced = TRUE)


crcRedMm_h13_amp2 <- crcRedMm_h13_amp
crcRedMm_h13_amp2$Freq <- 100
crcRedHg_h13_amp2 <- crcRedHg_h13_amp
crcRedHg_h13_amp2$Freq <- 100

crcExampleSingleH13 <- circosFreq2(crcRedMm_h13_amp2,
                                   crcRedMm_h13_del,
                                   crcRedHg_h13_amp2,
                                   crcRedHg_h13_del,
                                   filename = "20250313syntenyH13Example",
                                   ref = "human", reduced = TRUE)

### 7p and 7q do seprately

crcSynTableGr <- GRanges(seqnames = str_replace_all(crcSynTable$h_chr, "h_chr", "chr"),
                         IRanges(crcSynTable$h_start, crcSynTable$h_end))
index7p <- subjectHits(findOverlaps(hg19ArmLocations_redGr[1], crcSynTableGr))
index7q <- subjectHits(findOverlaps(hg19ArmLocations_redGr[2], crcSynTableGr))


crcRedHg_h7p <- crcRedHg[index7p, ]
crcRedMm_h7p <- crcRedMm[index7p, ]

crcRedHg_h7p_amp <- crcRedHg_h7p[which(crcRedHg_h7p$Freq > 0),]
crcRedHg_h7p_del <- crcRedHg_h7p[which(crcRedHg_h7p$Freq < 0),]
crcRedHg_h7p_del$Freq <- 0

crcRedMm_h7p_amp <- crcRedMm_h7p[1:(nrow(crcRedMm_h7p)/2), ]
crcRedMm_h7p_del <- crcRedMm_h7p[(nrow(crcRedMm_h7p)/2):nrow(crcRedMm_h7p), ]
crcRedMm_h7p_del$Freq <- 0


crcExampleSingleh7p <- circosFreq2(crcRedMm_h7p_amp,
                                   crcRedMm_h7p_del,
                                   crcRedHg_h7p_amp,
                                   crcRedHg_h7p_del,
                                   filename = "20250214syntenyh7p",
                                   ref = "human", reduced = TRUE)



crcRedHg_h7q <- crcRedHg[index7q, ]
crcRedMm_h7q <- crcRedMm[index7q, ]

crcRedHg_h7q_amp <- crcRedHg_h7q[which(crcRedHg_h7q$Freq > 0),]
crcRedHg_h7q_del <- crcRedHg_h7q[which(crcRedHg_h7q$Freq < 0),]
crcRedHg_h7q_del$Freq <- 0

crcRedMm_h7q_amp <- crcRedMm_h7q[1:(nrow(crcRedMm_h7q)/2), ]
crcRedMm_h7q_del <- crcRedMm_h7q[(nrow(crcRedMm_h7q)/2):nrow(crcRedMm_h7q), ]
crcRedMm_h7q_del$Freq <- 0


crcExampleSingleh7q <- circosFreq2(crcRedMm_h7q_amp,
                                   crcRedMm_h7q_del,
                                   crcRedHg_h7q_amp,
                                   crcRedHg_h7q_del,
                                   filename = "20250214syntenyh7q",
                                   ref = "human", reduced = TRUE)


### same but for 3p: bap1 example

index3p <- subjectHits(findOverlaps(hg19ArmLocationsGr[5], crcSynTableGr))

crcRedHg_h3p <- crcRedHg[index3p, ]
crcRedMm_h3p <- crcRedMm[index3p, ]

crcRedHg_h3p_amp <- crcRedHg_h3p[which(crcRedHg_h3p$Freq > 0),]
crcRedHg_h3p_del <- crcRedHg_h3p[which(crcRedHg_h3p$Freq < 0),]
crcRedHg_h3p_del$Freq <- 0

crcRedMm_h3p_amp <- crcRedMm_h3p[1:(nrow(crcRedMm_h3p)/2), ]
crcRedMm_h3p_del <- crcRedMm_h3p[(nrow(crcRedMm_h3p)/2):nrow(crcRedMm_h3p), ]
crcRedMm_h3p_del$Freq <- 0


crcExampleSingleh3p <- circosFreq2(crcRedMm_h3p_amp,
                                   crcRedMm_h3p_del,
                                   crcRedHg_h3p_amp,
                                   crcRedHg_h3p_del,
                                   filename = "20250528syntenyh3p",
                                   ref = "human", reduced = TRUE)

### make one for mouse to show smad4 and apc

crcRedHg_m18 <- crcRedHg[which(crcRedMm$Chr == "18"), ]
crcRedMm_m18 <- crcRedMm[which(crcRedMm$Chr == "18"), ]

crcRedHg_m18_amp <- crcRedHg_m18[which(crcRedHg_m18$Freq > 0),]
crcRedHg_m18_del <- crcRedHg_m18[which(crcRedHg_m18$Freq < 0),]
# crcRedHg_m18_del$Freq <- 0

crcRedMm_m18_amp <- crcRedMm_m18[1:(nrow(crcRedMm_m18)/2), ]
crcRedMm_m18_del <- crcRedMm_m18[(nrow(crcRedMm_m18)/2):nrow(crcRedMm_m18), ]
# crcRedMm_m18_del$Freq <- 0


crcExampleSingleM18 <- circosFreq2(crcRedMm_m18_amp,
                                   crcRedMm_m18_del,
                                   crcRedHg_m18_amp,
                                   crcRedHg_m18_del,
                                   filename = "20250214syntenyM18",
                                   ref = "mouse", reduced = TRUE)


### getting stats for paper i.e % of mouse to human


crcSynTable2 <- crcSynTable
crcSynTable2$str2 <- paste(crcSynTable2$h_chr, crcSynTable2$h_start,crcSynTable2$h_end,
                           crcSynTable2$m_chr, crcSynTable2$m_start, crcSynTable2$m_end)

crcSynTable2 <- crcSynTable2[-which(duplicated(crcSynTable2$str2)), ]
crcSynTableGr <- GRanges(seqnames = str_remove(crcSynTable2$h_chr, "h_"),
                         IRanges(start = crcSynTable2$h_start, end = crcSynTable2$h_end))



tmpOv <- findOverlaps(hg19ArmLocationsGr, crcSynTableGr)

crcSynTable2$hg19Arm[subjectHits(tmpOv)] <- hg19ArmLocations$str[queryHits(tmpOv)]

crcSynTable2$str <- paste(crcSynTable2$hg19Arm, crcSynTable2$m_chr)
crcSynTable2$mm10Length <- crcSynTable2$m_end - crcSynTable2$m_start

sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7p m_chr5")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr5")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7p m_chr6")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr6")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7p m_chr9")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr9")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7p m_chr11")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr11")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7p m_chr12")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr12")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7p m_chr13")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr13")])

sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7q m_chr5")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr5")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7q m_chr6")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr6")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "7q m_chr12")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr12")])


sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "13q m_chr1")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr1")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "13q m_chr3")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr3")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "13q m_chr5")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr5")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "13q m_chr8")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr8")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "13q m_chr14")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr14")])


sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "20p m_chr2")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr2")])
sum(crcSynTable2$mm10Length[which(crcSynTable2$str == "20q m_chr2")])/sum(mm10ChromArms$length[which(mm10ChromArms$chrom == "chr2")])

