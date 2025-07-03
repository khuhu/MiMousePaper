### Dev version in /home/kevhu/scripts/20201014mouseAaConversion.R
### note some functions may be depreciated, or improved for final version

### required libraries and tables

library(liftOver)
library(GenomicRanges)
library(liftOver)
library(Biostrings)
library(msa)

data("BLOSUM100")

h38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

hg38biomartTable <- read.table("/home/kevhu/data/20201020hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20201020Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")


geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")

hg38Peptide <- read.table("/home/kevhu/data/20201022proteinHg.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")

mm10Peptide <- read.table("/home/kevhu/data/20201022proteinMm.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")



### function to convert genomic position to relative cds

cdsConversion <- function(df){
  if(df$strand[1] == 1){
    df$exon_chrom_start_strand <- df$exon_chrom_start
    df$exon_chrom_end_strand <- df$exon_chrom_end
  } else{
    df$exon_chrom_start_strand <- df$exon_chrom_end
    df$exon_chrom_end_strand <- df$exon_chrom_start
  }
  
  strand_cds_start_conv <- NULL
  strand_cds_end_conv <- NULL
  lastExon <- max(df$rank)
  for (i in 1:nrow(df)) {
    if(df$rank[i] == 1){
      strand_cds_start_conv <- c(strand_cds_start_conv,
                                 df$exon_chrom_start_strand[i] - 1 + df$cdna_coding_start[i])
      strand_cds_end_conv <- c(strand_cds_end_conv, df$exon_chrom_end_strand[i])
    } else if(df$rank[i] == lastExon){
      strand_cds_start_conv <- c(strand_cds_start_conv, df$exon_chrom_start_strand[i])
      strand_cds_end_conv <- c(strand_cds_end_conv,
                               df$exon_chrom_start_strand[i] - 1 + df$cdna_coding_end[i])
    } else{
      strand_cds_start_conv <- c(strand_cds_start_conv, df$exon_chrom_start_strand[i])
      strand_cds_end_conv <- c(strand_cds_end_conv, df$exon_chrom_end_strand[i])
    }
  }
  df$strand_cds_start_conv<- strand_cds_start_conv
  df$strand_cds_end_conv <- strand_cds_end_conv
  return(df)
}


### converts position of amino acid for canonical transcript into genomic coordinates
aaToGenome <- function(gene, position){
  tmpTable <- hg38biomartTable[which(hg38biomartTable$external_gene_name == gene),]
  tmpTable <- cdsConversion(tmpTable)
  aaCdsEnd <- position * 3
  aaCdsStart <- aaCdsEnd - 2
  aaRange <- IRanges(start = aaCdsStart, end = aaCdsEnd)
  cdsRange <- IRanges(start = tmpTable$cds_start, end = tmpTable$cds_end)
  aaExon <- tmpTable[subjectHits(findOverlaps(aaRange, cdsRange)),]
  if(aaExon$strand == 1){
    tmpPos <- aaExon$strand_cds_start_conv + (aaCdsStart - aaExon$cds_start) - 1
  } else{
    tmpPos <- aaExon$strand_cds_end_conv - (aaCdsEnd - aaExon$cds_end) - 1
  }
  
  genomeCoords <- list("chrom" = aaExon$chromosome_name,"start" = tmpPos,
                       "end" = tmpPos + 2, "strand" = aaExon$strand, 
                       "ensembl" = aaExon$ensembl_transcript_id,
                       "cdna_coding_start" = aaExon$cdna_coding_start,
                       "cdna_coding_end" = aaExon$cdna_coding_end,
                       "gene" = aaExon$external_gene_name)
  
  #genomeCoords <- GRanges(seqnames = paste0("chr", aaExon$chromosome_name),
  #                        ranges = IRanges(start = tmpPos, end = tmpPos + 2),
  #                        strand = aaExon$strand)
  
  return(genomeCoords)
}

### reverse process of genomic coordinates to amino acid position

genomeToAA <- function(liftOverRes){
  ### make this work for any species to any species at a later date
  tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == liftOverRes$gene)]
  genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name == tmpGeneName),]
  genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
  
  cdsRangeBp <- IRanges(start = liftOverRes$start, end = liftOverRes$end)
  cdsExonRangeBp <- IRanges(start = genomeCdsTmpTable$exon_chrom_start, end = genomeCdsTmpTable$exon_chrom_end)
  cdsExon <- genomeCdsTmpTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]
  
  ### the subtraction of end/start depends on strand
  tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end))) - 1
  tmpAaPos <- ceiling(tmpPos/3)
  
  ### might need to return more
  return(tmpAaPos)
}


### lifting over cds regions between hg19 and mm10

cdsLiftOver <- function(cdsToGenomeOut, chainFile){
  out <- liftOver(GRanges(seqnames = paste0("chr", unlist(cdsToGenomeOut$chrom)),
                          ranges = IRanges(start = unlist(cdsToGenomeOut$start), end = unlist(cdsToGenomeOut$end)),
                          strand = unlist(cdsToGenomeOut$strand)),
                  chain = chainFile)
  
  
  out2 <- list("chrom" = as.character(out@unlistData@seqnames[1]),
               "start" = out@unlistData@ranges@start,
               "end" = out@unlistData@ranges@start + 2,
               "strand" = as.character(out@unlistData@strand@values),
               "gene" = unlist(cdsToGenomeOut$gene))
  
  return(out2)
}




substrMiddlePep <- function(prt1, prt2, adjacentLength, aaPos1, aaPos2, aaArea){
  if (aaArea == "end") {
    sub_prt1 <- substr(prt1, aaPos1 - adjacentLength, aaPos1)
    sub_prt2 <- substr(prt2, aaPos2 - adjacentLength, aaPos2)
  } else if(aaArea == "start"){
    sub_prt1 <- substr(prt1, aaPos1, aaPos1 + adjacentLength)
    sub_prt2 <- substr(prt2, aaPos2, aaPos2 + adjacentLength)
  } else{
    sub_prt1 <- substr(prt1, aaPos1 - ceiling(adjacentLength/2), aaPos1 + ceiling(adjacentLength/2))
    sub_prt2 <- substr(prt2, aaPos2 - ceiling(adjacentLength/2), aaPos2 + ceiling(adjacentLength/2))
  }
  
  res <- list("sub_prt1" = sub_prt1, "sub_prt2" = sub_prt2)
  return(res)
}

### check for same amino acid at the converted position

aaConversionCheck <- function(coordinatesList){
  ### get proteins then subset based on the length of the shortest peptide - with a minimum of 10 probably
  convertedAaPos <- genomeToAA(checkInput[[3]])
  originalAaPos <- coordinatesList[[2]]
  originalGene <- coordinatesList[[1]]
  
  newGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == checkInput[[3]]$gene)]
  
  originalPeptide <- hg38Peptide$peptide[which(hg38Peptide$external_gene_name == originalGene)]
  otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == newGeneName)]
  
  aaLength <- min(nchar(originalPeptide), nchar(otherPeptide))
  aaSurround <- aaLength * 0.05
  ### also need to check if they amino acids are on the ends/start of the peptide
  if (aaSurround < 10) {
    if((originalAaPos - 10) < 1){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                     originalAaPos, convertedAaPos, "start")
    } else if ((originalAaPos + 10) < aaLength){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                     originalAaPos, convertedAaPos, "end")
    } else{
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10,
                                     originalAaPos, convertedAaPos, "middle")
    }
  } else{
    if((originalAaPos - aaSurround) < 1){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "start")
    } else if ((originalAaPos + aaSurround) < aaLength){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "end")
    } else{
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "middle")
    }
  }
  
  aaCheckOriginal <- substr(originalPeptide, originalAaPos,originalAaPos)
  aaCheckConverted <- substr(otherPeptide, convertedAaPos, convertedAaPos)
  matchingCheck <- ifelse(aaCheckConverted == aaCheckOriginal, "yes", "no")
  
  msaRes <- msa(AAStringSet(c(adjacentAas$sub_prt1,adjacentAas$sub_prt2)), type = "protein")
  conScore <- msaConservationScore(msaRes, substitutionMatrix = BLOSUM100)
  ConSeq <- msaConsensusSequence(msaRes)
  
  return(list("seq1" = adjacentAas$sub_prt1,"seq2"=adjacentAas$sub_prt1, "conSeq"=ConSeq,
              "check" = matchingCheck, "conservationScore"=conScore))
}



