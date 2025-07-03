### original script at: /home/kevhu/scripts/20230221greedyHoyspot.R 


### previous version only looks at canonical transcript
### greedy version compares all transcript from one to another


source("hotspotConversionFunctions.R")

library(stringr)

hg19toHg38Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg19ToHg38.over.chain")

h38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

hg38biomartTable <- read.table("/home/kevhu/data/20230320hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20210203Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")


geneNameDf <- read.table("/home/kevhu/data/20230315geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")

# hg38Peptide <- read.table("/home/kevhu/data/20201030proteinHg.txt",  header = TRUE,
#                           stringsAsFactors = FALSE, sep = "\t")

hg38Peptide <- read.table("/mnt/DATA5/tmp/kev/misc/20230301allPeptideQuery.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")

mm10Peptide <- read.table("/mnt/DATA5/tmp/kev/misc/20230320allPeptideQueryMm10.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")


uniprotWsTableMm10 <- read.table("/mnt/DATA5/tmp/kev/misc/20230320mm10UniprotWsTable.txt",
                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

uniprotWsTableHg38 <- read.table("/mnt/DATA5/tmp/kev/misc/20230315hg38UniprotWsTable.txt",
                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)


small_cosmic <- read.table("/home/kevhu/data/20190225mouseConvertedOncogeneHotspots.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hotspots <- unique(small_cosmic$Hotspot.Residue)

hotspotDf <- NULL
for (i in hotspots) {
  tmpVar <- strsplit(i, "\\|")
  tmpGene <- tmpVar[[1]][1]
  tmpVar2 <- tmpVar[[1]][2]
  tmpAa <- substr(tmpVar2, 1, 1)
  tmpPos <- substr(tmpVar2, 2, nchar(tmpVar2))
  tmpLine <- c(tmpGene, tmpAa, tmpPos)
  hotspotDf <- rbind(hotspotDf, tmpLine)
}

rownames(hotspotDf) <- NULL
colnames(hotspotDf) <- c("Gene", "Aa", "Position")
hotspotDf <- data.frame(hotspotDf, stringsAsFactors = FALSE)


hotspotDf <- hotspotDf[-which(hotspotDf == "GNAS"),] # different transcript ID
# hotspotDf <- hotspotDf[-which(hotspotDf == "GNA11"),]
hotspotDf <- hotspotDf[-which(hotspotDf == "GNAQ"),] # the name searching method is off -----------------------fixed
hotspotDf <- hotspotDf[-which(hotspotDf == "H3F3A"),] # Lower case is mouse name H3F3A, hg38 name is H3-3A
hotspotDf <- hotspotDf[-which(hotspotDf == "DNMT3A"),] # not canonical isoform? - liftOver empty
hotspotDf <- hotspotDf[-which(hotspotDf == "FOXL2"),] # msa error - no seq2 - missing mouse peptide
hotspotDf <- hotspotDf[-which(hotspotDf == "HIST1H3B"),] # not the gene name

### so to redesign portions. as mentioned before, first (1) the gene name change should be easy enough to implement
### (2) fix other issues like msa matching and what-not (3) the wrong transcript-isoform being useds
### the last portion is the hardest to deal with since, I need to get all the different isoforms ... and figure out
### an efficient way to query them - for the time being just skip ones like this



### testing lollipop plot scirpt for betaBatenin in humans and mouse overview

library(drawProteins)
betaCateninMm <- getBM(mart = ensemblMm, attributes = c("external_gene_name", "ensembl_transcript_id", 'uniprot_gn_symbol', 'uniprot_gn_id'),
      filters = "ensembl_transcript_id", values = c("ENSMUST00000108658"))

betaCateninHg <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                        "ensembl_transcript_id_version", 'uniprot_gn_symbol',
                                                        'uniprot_gn_id', "uniprot_isoform","transcript_mane_select"),
                       filters = "ensembl_transcript_id", values = geneNameDf$ensembl_transcript_id[which(geneNameDf$external_gene_name == "BRCA1")])

betaCateninHg <- betaCateninHg[-which(betaCateninHg$transcript_mane_select == ""),]

betaCateninHg <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                        "ensembl_transcript_id_version",'uniprot_gn_id',
                                                        "uniprotswissprot", "uniprotsptrembl"),
                       filters = "ensembl_transcript_id", values = geneNameDf$ensembl_transcript_id[which(geneNameDf$external_gene_name == "BRCA1")])

### custom function for lollipop
### need another function to add features to the object

### I need to redo the functions b/c for some odd reason uniprot canonical protein
### doesn't include 

drawProteins::get_features(paste(betaCateninHg$uniprot_gn_id[1:21], collapse = " ")) -> rel_json_hg
# drawProteins::get_features(paste(betaCateninHg$uniprot_gn_id, collapse = " ")) -> rel_json_hg
drawProteins::feature_to_dataframe(rel_json_hg) -> rel_data_hg
draw_canvas(rel_data_hg) -> p
p <- draw_chains(p, rel_data_hg)
p <- draw_domains(p, rel_data_hg)
# p <- draw_phospho(p, rel_data_hg)
p <- draw_motif(p, rel_data_hg)
p <- p + theme_bw(base_size = 20) + # white background
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p


drawProteins::get_features(paste(betaCateninMm$uniprot_gn_id, collapse = " ")) -> rel_json_mm
drawProteins::feature_to_dataframe(rel_json_mm) -> rel_data_mm
draw_canvas(rel_data_mm) -> p
p <- draw_chains(p, rel_data_mm)
p <- draw_domains(p, rel_data_mm)
p <- draw_phospho(p, rel_data_mm)
p <- p + theme_bw(base_size = 20) + # white background
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p

### using both human and mosue and then fake phospho sites for marker

lollipop_site_info <- function(features)
{
  features <- features[features$type == "LOLLIPOP", ]
  lolli_list <- grep("Phospho", features$description)
  lolli_features <- features[phospho_list, ]
  return(phospho_features)
}

drawLollipop <- function (p, data = data, size = 2, fill = "yellow", show.legend = FALSE) {
  begin = end = description = NULL
  p <- p + ggplot2::geom_point(data = drawProteins::phospho_site_info(data), 
                               ggplot2::aes(x = begin, y = order + 0.25), shape = 21, 
                               colour = "black", fill = fill, size = size, show.legend = show.legend)
  return(p)
}


drawProteins::get_features(paste(c("K7PPA8", "Q549C9"), collapse = " ")) -> rel_json_hg
drawProteins::feature_to_dataframe(rel_json_hg) -> rel_data_hg

### can't just rbind  it, need to add it to the s3 layer
rel_data_hg <- rbind(rel_data_hg, data.frame("MOD_RES", "Hotspot mutation", 273, 273, 1,"K7PPA8", "K7PPA8_HUMAN", 9606, 1))

draw_canvas(rel_data_hg) -> p 
p <- draw_chains(p, rel_data_hg)
p <- draw_domains(p, rel_data_hg)
p <- p + theme_bw(base_size = 20) + # white background
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p



### heuristic  for create table of specific uniprot id, one used is longest with phosphorylation data

tmpTable <- read.table("/mnt/DATA5/tmp/kev/misc/iprscan5-R20230221-202508-0275-70018558-p2m.tsv", sep = "\t",
                       stringsAsFactors = FALSE, header = FALSE, fill = TRUE)



tmpTable <- read.table("/mnt/DATA5/tmp/kev/interpro/20230223hotspotGenes.fasta.tsv", sep = "\t",
                       stringsAsFactors = FALSE, header = FALSE, fill = TRUE)


### need to create file of fasta files for all proteins for hotspot genes for now
### can get the "canonical" annotation from uniprot based on entry name
library(seqinr)
unique(hotspotDf$Gene)

length(unique(geneNameDf$external_gene_name))

# hotspotPepQuery5k <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id", "uniprot_gn_id",
#                                                           "transcript_mane_select", "peptide"),
#                          filters = "external_gene_name", values = unique(geneNameDf$external_gene_name)[1:500])


tmpVector <- c(1:196)
finalHotspotPepQuery <- NULL
for (i in tmpVector) {
  tmpStart <- (i - 1) * 100 + 1
  tmpEnd <- i * 100
  print(paste0(tmpStart, ":", tmpEnd))
  if (i == 196) {
    tmpQuery <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id", "uniprot_gn_id",
                                                       "transcript_mane_select", "peptide"),
                      filters = "external_gene_name", values = unique(geneNameDf$external_gene_name)[19500:19564])
    finalHotspotPepQuery <- rbind(finalHotspotPepQuery, tmpQuery)
  } else{
    tmpQuery <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id", "uniprot_gn_id",
                                                       "transcript_mane_select", "peptide", "mmusculus_homolog_ensembl_peptide"),
                      filters = "external_gene_name", values = unique(geneNameDf$external_gene_name)[tmpStart:tmpEnd])
    finalHotspotPepQuery <- rbind(finalHotspotPepQuery, tmpQuery)
  }
}


finalHotspotPepQuery$tmp <- paste(finalHotspotPepQuery$external_gene_name,
                                  finalHotspotPepQuery$ensembl_transcript_id, finalHotspotPepQuery$uniprot_gn_id, sep = "|")

finalHotspotPepQuery$peptide <- stringr::str_remove_all(finalHotspotPepQuery$peptide, "\\*")

write.fasta(sequences = as.list(finalHotspotPepQuery$peptide), names = finalHotspotPepQuery$tmp,
            file.out = "/mnt/DATA5/tmp/kev/misc/20230301allHotspotGenes.fasta")

write.table(finalHotspotPepQuery, "/mnt/DATA5/tmp/kev/misc/20230301allPeptideQuery.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

### needed to do the same thing for mouse 

# tmpVectorMm <- c(1:187)
# finalHotspotPepQueryMm <- NULL
# for (i in tmpVectorMm) {
#   tmpStart <- (i - 1) * 100 + 1
#   tmpEnd <- i * 100
#   print(paste0(tmpStart, ":", tmpEnd))
#   if (i == 187) {
#     tmpQuery <- getBM(mart = ensemblMm, attributes = c("external_gene_name", "ensembl_transcript_id",
#                                                        "uniprot_gn_id", "peptide", "ensembl_peptide_id"),
#                       filters = "external_gene_name", values = unique(geneNameDf$mmusculus_homolog_associated_gene_name)[18600:18681])
#     finalHotspotPepQueryMm <- rbind(finalHotspotPepQueryMm, tmpQuery)
#   } else{
#     tmpQuery <- getBM(mart = ensemblMm, attributes = c("external_gene_name", "ensembl_transcript_id",
#                                                        "uniprot_gn_id", "peptide", "ensembl_peptide_id"),
#                       filters = "external_gene_name", values = unique(geneNameDf$mmusculus_homolog_associated_gene_name)[tmpStart:tmpEnd])
#     finalHotspotPepQueryMm <- rbind(finalHotspotPepQueryMm, tmpQuery)
#   }
# }
allPossibleNames <- unique(c(unique(geneNameDf$mmusculus_homolog_associated_gene_name),  unique(geneNameDf$external_gene_name)))


tmpVectorMm <- c(1:383)
finalHotspotPepQueryMm <- NULL
for (i in tmpVectorMm) {
  tmpStart <- (i - 1) * 100 + 1
  tmpEnd <- i * 100
  print(paste0(tmpStart, ":", tmpEnd))
  if (i == 383) {
    tmpQuery <- getBM(mart = ensemblMm, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                       "uniprot_gn_id", "peptide", "ensembl_peptide_id"),
                      filters = "external_gene_name", values = allPossibleNames[38200:38209])
    finalHotspotPepQueryMm <- rbind(finalHotspotPepQueryMm, tmpQuery)
  } else{
    tmpQuery <- getBM(mart = ensemblMm, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                       "uniprot_gn_id", "peptide", "ensembl_peptide_id"),
                      filters = "external_gene_name", values = allPossibleNames[tmpStart:tmpEnd])
    finalHotspotPepQueryMm <- rbind(finalHotspotPepQueryMm, tmpQuery)
  }
}


finalHotspotPepQueryMm$tmp <- paste(finalHotspotPepQueryMm$external_gene_name,
                                  finalHotspotPepQueryMm$ensembl_transcript_id, finalHotspotPepQueryMm$uniprot_gn_id, sep = "|")

finalHotspotPepQueryMm$peptide <- stringr::str_remove_all(finalHotspotPepQueryMm$peptide, "\\*")

finalHotspotPepQueryMm <- finalHotspotPepQueryMm[-which(duplicated(finalHotspotPepQueryMm$tmp)),]

write.fasta(sequences = as.list(finalHotspotPepQueryMm$peptide), names = finalHotspotPepQueryMm$tmp,
            file.out = "/mnt/DATA5/tmp/kev/misc/20230320allHotspotGenesMm10.fasta")

write.table(finalHotspotPepQueryMm, "/mnt/DATA5/tmp/kev/misc/20230320allPeptideQueryMm10.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)


### recreating the geneDf table but with two more columns

tmpVector <- c(1:196)
geneQuery <- NULL
for (i in tmpVector) {
  tmpStart <- (i - 1) * 100 + 1
  tmpEnd <- i * 100
  print(paste0(tmpStart, ":", tmpEnd))
  if (i == 196) {
    tmpQuery <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                       "mmusculus_homolog_ensembl_peptide", "mmusculus_homolog_associated_gene_name"),
                      filters = "external_gene_name", values = unique(hg38Peptide$external_gene_name)[19500:19564])
    geneQuery<- rbind(geneQuery, tmpQuery)
  } else{
    tmpQuery <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                       "mmusculus_homolog_ensembl_peptide", "mmusculus_homolog_associated_gene_name"),
                      filters = "external_gene_name", values = unique(geneNameDf$external_gene_name)[tmpStart:tmpEnd])
    geneQuery <- rbind(geneQuery, tmpQuery)
  }
}


geneQuery$transcript_mane_select <- ""
geneQuery$transcript_mane_select <- hg38Peptide$transcript_mane_select[match(geneQuery$ensembl_transcript_id,
                                                                             hg38Peptide$ensembl_transcript_id)]

write.table(geneQuery, "/home/kevhu/data/20230315geneNameBiomart.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")





### probably should redo above by with all genes so I can make a diagram  for p53
### doing this for all genes in oncokb

interproDomanins <- read.table("/mnt/DATA5/tmp/kev/interpro/20230301allHotspotGenes.fasta.tsv", header = FALSE,
                               stringsAsFactors = FALSE, sep = "\t", quote="")


annoTable2 <- read.table("/mnt/DATA5/tmp/kev/oncokbAnnoOut/mskTestGRCh38.maf", sep = "\t",
                         stringsAsFactors = FALSE, header = TRUE, fill = TRUE)


### idea would 

# egfrPositions <- annoTable2[which(annoTable2$Hugo_Symbol == "EGFR"), ]
# egfrPositions <- egfrPositions[-grep("dup|del|ins", egfrPositions$cDNA_change), ]
# egfrPositions <- egfrPositions[-grep("splic", egfrPositions$HGVSp_Short), ]
# egfrPositions <- egfrPositions[which(egfrPositions$MUTATION_EFFECT %in% c("Gain-of-function", "Likely Gain-of-function")),]


# oncokbPositions <- annoTable2
# oncokbPositions <- oncokbPositions[-grep("dup|del|ins", oncokbPositions$cDNA_change), ]
# oncokbPositions <- oncokbPositions[-grep("splic", oncokbPositions$HGVSp_Short), ]

oncokbPositions <- annoTable2[which(annoTable2$Variant_Type == "SNP"),]



### using all and then I will label end matrix with the mutation effect
# oncokbPositions <- oncokbPositions[which(oncokbPositions$MUTATION_EFFECT %in% c("Gain-of-function", "Likely Gain-of-function")),]


# origAa <- NULL
# for (i in 1:nrow(egfrPositions)) {
#   tmpRes <- substr(egfrPositions$HGVSp_Short[i], 3, 3)
#   origAa  <- c(origAa , tmpRes)
# }
# origAa
# 
# egfrPositions$AminoAcid <- origAa
# uniquePositions <- unique(egfrPositions$Protein_position)
# egfrPositionsRed <- egfrPositions[-which(duplicated(egfrPositions$Protein_position)),]
# 
# hotspotDf <- data.frame("Gene" = egfrPositionsRed$Hugo_Symbol,
#                              "Aa" = egfrPositionsRed$AminoAcid, "Position" = egfrPositionsRed$Protein_position)


origAa <- NULL
for (i in 1:nrow(oncokbPositions)) {
  tmpRes <- substr(oncokbPositions$HGVSp_Short[i], 3, 3)
  origAa  <- c(origAa , tmpRes)
}
origAa

oncokbPositions$AminoAcid <- origAa
oncokbPositionsHg38 <- oncokbPositions[which(oncokbPositions$NCBI_Build == "GRCh38"),]
oncokbPositionsHg19 <- oncokbPositions[which(oncokbPositions$NCBI_Build == "GRCh37"),]
colnames(oncokbPositionsHg19)[c(1, 5:8)] <- c("gene", "chrom", "start", "end", "strand")

### redo later
# oncokbPositionsLift <- NULL
# for (j in 1:nrow(oncokbPositionsHg19)) {
#   print(j)
#   tmpLiftOverRes <- cdsLiftOver(oncokbPositionsHg19[j, ],
#                                 chainFile = hg19toHg38Chain)
# 
#   if (is.null(tmpLiftOverRes)) {
#     tmpDf <- data.frame("chrom" = NA,
#                         "start" = NA,
#                         "end" = NA,
#                         "strand" = NA,
#                         "gene" = oncokbPositionsHg19$gene[j])
#     oncokbPositionsLift <- rbind(oncokbPositionsLift,tmpDf)
#     next()
#   }
# 
#   oncokbPositionsLift <- rbind(oncokbPositionsLift,tmpLiftOverRes)
# }

library(foreach)
library(doParallel)
cl <- makeCluster(30)
registerDoParallel(cl)
oncokbPositionsLift <- NULL
oncokbPositionsLift <- foreach(i = 1:nrow(oncokbPositionsHg19), .combine = "rbind", .packages = c('GenomicRanges', 'liftOver')) %dopar% {
  
  
  cdsLiftOver <- function(cdsToGenomeOut, chainFile){
    out <- liftOver(GRanges(seqnames = paste0("chr", unlist(cdsToGenomeOut$chrom)),
                            ranges = IRanges(start = unlist(cdsToGenomeOut$start), end = unlist(cdsToGenomeOut$end)),
                            strand = unlist(cdsToGenomeOut$strand)),
                    chain = chainFile)
    
    if (isEmpty(out)) {
      return()
    }
    
    out2 <- data.frame("chrom" = as.character(out@unlistData@seqnames[1]),
                       "start" = out@unlistData@ranges@start,
                       "end" = out@unlistData@ranges@start + 2,
                       "strand" = as.character(out@unlistData@strand@values),
                       "gene" = unlist(cdsToGenomeOut$gene), stringsAsFactors = FALSE)
    
    return(out2)
  }
  
  tmpLiftOverRes <- cdsLiftOver(oncokbPositionsHg19[i, ],
                                chainFile = hg19toHg38Chain)
  
  if (is.null(tmpLiftOverRes)) {
    tmpLiftOverRes <- data.frame("chrom" = NA,
                        "start" = NA,
                        "end" = NA,
                        "strand" = NA,
                        "gene" = oncokbPositionsHg19$gene[i])
    }
  tmpLiftOverRes
}

stopCluster(cl)

oncokbPositionsHg19$chrom <- oncokbPositionsLift$chrom
oncokbPositionsHg19$start <- oncokbPositionsLift$start
oncokbPositionsHg19$end <- oncokbPositionsLift$end
oncokbPositionsHg19 <- oncokbPositionsHg19[-which(is.na(oncokbPositionsHg19$chrom)),]

colnames(oncokbPositionsHg19) <- colnames(oncokbPositionsHg38)
oncokbPositionsAllHg38 <- rbind(oncokbPositionsHg38, oncokbPositionsHg19)

oncokbPositionsAllHg38$tmp <- paste0(oncokbPositionsAllHg38$Hugo_Symbol, oncokbPositionsAllHg38$Protein_position)
oncokbPositionsRed <- oncokbPositionsAllHg38[-which(duplicated(oncokbPositionsAllHg38$tmp)),]
# oncokbPositionsRed <- oncokbPositionsAllHg38[-which(duplicated(oncokbPositionsAllHg38$Protein_position)),]

hotspotDf <- data.frame("Gene" = oncokbPositionsRed$Hugo_Symbol,"Aa" = oncokbPositionsRed$AminoAcid,
                        "Position" = oncokbPositionsRed$Protein_position, "MutEff" = oncokbPositionsRed$MUTATION_EFFECT,
                        "Consequence" = oncokbPositionsRed$Consequence, "OncokbListed" = oncokbPositionsRed$VARIANT_IN_ONCOKB)

hotspotDf <- hotspotDf[hotspotDf$OncokbListed == "True", ]
rownames(hotspotDf) <- NULL

table(hotspotDf$Consequence)


i <- 20
### single, non parallel
# resTbl <- NULL
# for (i in 1:nrow(hotspotDf)) {
#   print(i)
#   testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
#                                           position = as.numeric(hotspotDf$Position[i])),
#                                error = function(x) return(NULL))
#   #print(testAaConversion)
#   #if (testAaConversion == "empty") {
#   #  next()
#   #}
#   if (is.null(testAaConversion)) {
#     next()
#   }
#   
#   liftOverOut <- NULL
#   for (j in 1:nrow(testAaConversion)) {
#     print(j)
#     tmpLiftOverRes <- cdsLiftOver(testAaConversion[j, ],
#                                   chainFile = h38toMm10Chain)
#     
#     if (is.null(tmpLiftOverRes)) {
#       next()
#     }
#     tmpLiftOverRes$ensembl <- testAaConversion$ensembl[j]
#     tmpLiftOverRes$aa <- hotspotDf$Aa[i]
#     liftOverOut <- rbind(liftOverOut,tmpLiftOverRes)
#   }
#   
#   # liftOverOut$ensembl <- testAaConversion$ensembl
#   # liftOverOut$aa <- hotspotDf$Aa[i]
#   #print(liftOverOut)
#   if (is.null(liftOverOut)) {
#     next()
#   }
#   
#   if(any(liftOverOut$strand == "+" | liftOverOut$strand == "-")){
#     liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\+", "1")
#     liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\-", "-1")
#     liftOverOut$strand <- as.numeric(liftOverOut$strand)
#   }
#   
#   
#   
#   ### when this is made into a function or a wrapper
#   ### make sure I have a variable for MANE. if true reduce set
#   mane_filt <- geneNameDf$transcript_mane_select[match(unique(liftOverOut$ensembl), geneNameDf$ensembl_transcript_id)]
#   if(grepl("NM", paste0(mane_filt, collapse = "|"))){
#     mane_transcript <- unique(liftOverOut$ensembl)[grep("NM", mane_filt)]
#     liftOverOut <- liftOverOut[which(liftOverOut$ensembl == mane_transcript),]
#     liftOverOut$homo_pep <- geneNameDf$mmusculus_homolog_ensembl_peptide[match(liftOverOut$ensembl,
#                                                                                geneNameDf$ensembl_transcript_id)]
#   } else{
#     liftOverOut$homo_pep <- ""
#   }
#   
#   ### iterating over all if MANE is not selected
#   for (y in 1:nrow(liftOverOut)) {
#     print(y)
#     checkInput <- list(hotspotDf$Gene[i],
#                        as.numeric(hotspotDf$Position[i]) ,
#                        "liftOverOut" = liftOverOut[y,],
#                        hotspotDf$MutEff[i])
#     #print(checkInput)
#     conversionCheckOut <- aaConversionCheck(checkInput)
#     #print(conversionCheckOut)
#     resTbl <- rbind(resTbl, conversionCheckOut)
#   }
# }


### below  is to look at all msk mutations, not just what was listed as recurrent on oncokb

resTbl <- NULL
cl <- makeCluster(30)
registerDoParallel(cl)
resTbl <- foreach(i = 1:10000, .combine = "rbind", .packages = c("liftOver", "GenomicRanges", "Biostrings", "msa", "stringr")) %dopar% {
  
  chroms <- c(1:22, "X", "Y")
  # source("/home/kevhu/scripts/20201020hotspotConversionFunctions.R", local = TRUE)
  testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
                                          position = as.numeric(hotspotDf$Position[i])),
                               error = function(x) return(NULL))
  testAaConversion <- testAaConversion[which(testAaConversion$chrom %in% chroms),]
  if (is.null(testAaConversion)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  } else if (nrow(testAaConversion) == 0) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  liftOverOut <- NULL
  for (j in 1:nrow(testAaConversion)) {
    print(j)
    tmpLiftOverRes <- cdsLiftOver(testAaConversion[j, ],
                                  chainFile = h38toMm10Chain)
    
    if (is.null(tmpLiftOverRes)) {
      next()
    }
    tmpLiftOverRes$ensembl <- testAaConversion$ensembl[j]
    tmpLiftOverRes$aa <- hotspotDf$Aa[i]
    liftOverOut <- rbind(liftOverOut,tmpLiftOverRes)
  }
  
  if (is.null(liftOverOut)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  if(any(liftOverOut$strand == "+" | liftOverOut$strand == "-")){
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\+", "1")
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\-", "-1")
    liftOverOut$strand <- as.numeric(liftOverOut$strand)
  }
  
  mane_filt <- geneNameDf$transcript_mane_select[match(unique(liftOverOut$ensembl), geneNameDf$ensembl_transcript_id)]
  if(grepl("NM", paste0(mane_filt, collapse = "|"))){
    mane_transcript <- unique(liftOverOut$ensembl)[grep("NM", mane_filt)]
    liftOverOut <- liftOverOut[which(liftOverOut$ensembl == mane_transcript),]
    liftOverOut$homo_pep <- geneNameDf$mmusculus_homolog_ensembl_peptide[match(liftOverOut$ensembl,
                                                                               geneNameDf$ensembl_transcript_id)]
  } else{
    liftOverOut$homo_pep <- ""
  }
  
  for (y in 1:nrow(liftOverOut)) {
    print(y)
    checkInput <- list(hotspotDf$Gene[i],
                       as.numeric(hotspotDf$Position[i]) ,
                       "liftOverOut" = liftOverOut[y,],
                       hotspotDf$MutEff[i])
    conversionCheckOut <- aaConversionCheck(checkInput)
  }
  conversionCheckOut 
}

closeAllConnections()


### 2 

resTbl2 <- NULL
cl <- makeCluster(30)
registerDoParallel(cl)
resTbl2 <- foreach(i = 10001:20000, .combine = "rbind", .packages = c("liftOver", "GenomicRanges", "Biostrings", "msa", "stringr")) %dopar% {
  
  chroms <- c(1:22, "X", "Y")
  # source("/home/kevhu/scripts/20201020hotspotConversionFunctions.R", local = TRUE)
  testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
                                          position = as.numeric(hotspotDf$Position[i])),
                               error = function(x) return(NULL))
  testAaConversion <- testAaConversion[which(testAaConversion$chrom %in% chroms),]
  if (is.null(testAaConversion)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  } else if (nrow(testAaConversion) == 0) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  liftOverOut <- NULL
  for (j in 1:nrow(testAaConversion)) {
    print(j)
    tmpLiftOverRes <- cdsLiftOver(testAaConversion[j, ],
                                  chainFile = h38toMm10Chain)
    
    if (is.null(tmpLiftOverRes)) {
      next()
    }
    tmpLiftOverRes$ensembl <- testAaConversion$ensembl[j]
    tmpLiftOverRes$aa <- hotspotDf$Aa[i]
    liftOverOut <- rbind(liftOverOut,tmpLiftOverRes)
  }
  
  if (is.null(liftOverOut)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  if(any(liftOverOut$strand == "+" | liftOverOut$strand == "-")){
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\+", "1")
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\-", "-1")
    liftOverOut$strand <- as.numeric(liftOverOut$strand)
  }
  
  mane_filt <- geneNameDf$transcript_mane_select[match(unique(liftOverOut$ensembl), geneNameDf$ensembl_transcript_id)]
  if(grepl("NM", paste0(mane_filt, collapse = "|"))){
    mane_transcript <- unique(liftOverOut$ensembl)[grep("NM", mane_filt)]
    liftOverOut <- liftOverOut[which(liftOverOut$ensembl == mane_transcript),]
    liftOverOut$homo_pep <- geneNameDf$mmusculus_homolog_ensembl_peptide[match(liftOverOut$ensembl,
                                                                               geneNameDf$ensembl_transcript_id)]
  } else{
    liftOverOut$homo_pep <- ""
  }
  
  for (y in 1:nrow(liftOverOut)) {
    print(y)
    checkInput <- list(hotspotDf$Gene[i],
                       as.numeric(hotspotDf$Position[i]) ,
                       "liftOverOut" = liftOverOut[y,],
                       hotspotDf$MutEff[i])
    conversionCheckOut <- aaConversionCheck(checkInput)
  }
  conversionCheckOut 
}

closeAllConnections()

### 3

resTbl3 <- NULL
cl <- makeCluster(30)
registerDoParallel(cl)
resTbl3 <- foreach(i = 20001:30000, .combine = "rbind", .packages = c("liftOver", "GenomicRanges", "Biostrings", "msa", "stringr")) %dopar% {
  
  chroms <- c(1:22, "X", "Y")
  # source("/home/kevhu/scripts/20201020hotspotConversionFunctions.R", local = TRUE)
  testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
                                          position = as.numeric(hotspotDf$Position[i])),
                               error = function(x) return(NULL))
  testAaConversion <- testAaConversion[which(testAaConversion$chrom %in% chroms),]
  if (is.null(testAaConversion)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  } else if (nrow(testAaConversion) == 0) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  liftOverOut <- NULL
  for (j in 1:nrow(testAaConversion)) {
    print(j)
    tmpLiftOverRes <- cdsLiftOver(testAaConversion[j, ],
                                  chainFile = h38toMm10Chain)
    
    if (is.null(tmpLiftOverRes)) {
      next()
    }
    tmpLiftOverRes$ensembl <- testAaConversion$ensembl[j]
    tmpLiftOverRes$aa <- hotspotDf$Aa[i]
    liftOverOut <- rbind(liftOverOut,tmpLiftOverRes)
  }
  
  if (is.null(liftOverOut)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  if(any(liftOverOut$strand == "+" | liftOverOut$strand == "-")){
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\+", "1")
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\-", "-1")
    liftOverOut$strand <- as.numeric(liftOverOut$strand)
  }
  
  mane_filt <- geneNameDf$transcript_mane_select[match(unique(liftOverOut$ensembl), geneNameDf$ensembl_transcript_id)]
  if(grepl("NM", paste0(mane_filt, collapse = "|"))){
    mane_transcript <- unique(liftOverOut$ensembl)[grep("NM", mane_filt)]
    liftOverOut <- liftOverOut[which(liftOverOut$ensembl == mane_transcript),]
    liftOverOut$homo_pep <- geneNameDf$mmusculus_homolog_ensembl_peptide[match(liftOverOut$ensembl,
                                                                               geneNameDf$ensembl_transcript_id)]
  } else{
    liftOverOut$homo_pep <- ""
  }
  
  for (y in 1:nrow(liftOverOut)) {
    print(y)
    checkInput <- list(hotspotDf$Gene[i],
                       as.numeric(hotspotDf$Position[i]) ,
                       "liftOverOut" = liftOverOut[y,],
                       hotspotDf$MutEff[i])
    conversionCheckOut <- aaConversionCheck(checkInput)
  }
  conversionCheckOut 
}

closeAllConnections()

### 4 

resTbl4 <- NULL
cl <- makeCluster(30)
registerDoParallel(cl)
resTbl4 <- foreach(i = 30001:42420, .combine = "rbind", .packages = c("liftOver", "GenomicRanges", "Biostrings", "msa", "stringr")) %dopar% {
  
  chroms <- c(1:22, "X", "Y")
  # source("/home/kevhu/scripts/20201020hotspotConversionFunctions.R", local = TRUE)
  testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
                                          position = as.numeric(hotspotDf$Position[i])),
                               error = function(x) return(NULL))
  testAaConversion <- testAaConversion[which(testAaConversion$chrom %in% chroms),]
  if (is.null(testAaConversion)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  } else if (nrow(testAaConversion) == 0) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  liftOverOut <- NULL
  for (j in 1:nrow(testAaConversion)) {
    print(j)
    tmpLiftOverRes <- cdsLiftOver(testAaConversion[j, ],
                                  chainFile = h38toMm10Chain)
    
    if (is.null(tmpLiftOverRes)) {
      next()
    }
    tmpLiftOverRes$ensembl <- testAaConversion$ensembl[j]
    tmpLiftOverRes$aa <- hotspotDf$Aa[i]
    liftOverOut <- rbind(liftOverOut,tmpLiftOverRes)
  }
  
  if (is.null(liftOverOut)) {
    conversionCheckOut <- data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                     "check" = NA, "conservationScore"= NA,
                                     "originalEns" = originalEns, "ENSID" = NA,
                                     "Hotspot" = paste0(hotspotDf$Gene[i], hotspotDf$Position[i], hotspotDf$Aa[i]),
                                     "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                     "liftOverChrom" = NA,
                                     "liftOverStart" = NA,
                                     "liftOverEnd" = NA,
                                     "annoMane" = "", stringsAsFactors = FALSE)
    return(conversionCheckOut)
  }
  
  if(any(liftOverOut$strand == "+" | liftOverOut$strand == "-")){
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\+", "1")
    liftOverOut$strand <- str_replace_all(liftOverOut$strand, "\\-", "-1")
    liftOverOut$strand <- as.numeric(liftOverOut$strand)
  }
  
  mane_filt <- geneNameDf$transcript_mane_select[match(unique(liftOverOut$ensembl), geneNameDf$ensembl_transcript_id)]
  if(grepl("NM", paste0(mane_filt, collapse = "|"))){
    mane_transcript <- unique(liftOverOut$ensembl)[grep("NM", mane_filt)]
    liftOverOut <- liftOverOut[which(liftOverOut$ensembl == mane_transcript),]
    liftOverOut$homo_pep <- geneNameDf$mmusculus_homolog_ensembl_peptide[match(liftOverOut$ensembl,
                                                                               geneNameDf$ensembl_transcript_id)]
  } else{
    liftOverOut$homo_pep <- ""
  }
  
  for (y in 1:nrow(liftOverOut)) {
    print(y)
    checkInput <- list(hotspotDf$Gene[i],
                       as.numeric(hotspotDf$Position[i]) ,
                       "liftOverOut" = liftOverOut[y,],
                       hotspotDf$MutEff[i])
    conversionCheckOut <- aaConversionCheck(checkInput)
  }
  conversionCheckOut 
}

closeAllConnections()

allResTbl <- rbind(resTbl, resTbl2, resTbl3, resTbl4)

### combine all 4 tables, each was about 131 gbs so needed to do them in parts

### make the graph before of length by msa score - reduce it by unique coordinates

### then can make the draw proteins graph - I think I might need my own functions 
### because I want to draw lollipops based on pop count

### filter by MANE annotation, then for each unique mouse transcript create graph
### filtered out b/c old version is all vs current version using only mutations with likely gain/loss of function

# rownames(resTbl) <- NULL
# resTblMane <- resTbl[-which(resTbl$annoMane == ""),]
rownames(allResTbl) <- NULL
resTblMane <- allResTbl[-which(allResTbl$annoMane == ""),]
rownames(resTblMane) <- NULL
# resTblMane <- rbind(resTblMane, resTblMane[61, ])
# resTblMane$check[66] <- "no"
resTblMane$check <- factor(resTblMane$check, levels = c("yes", "no"))
ggplot(resTblMane, aes(x = ENSID, y = conservationScore, fill = check)) + geom_boxplot() +
  geom_point(aes(fill = check), position = position_jitterdodge(), alpha = 1) +
  scale_fill_manual(values = c("darkgreen", "darkblue")) + ylab("Alignment score (VTML10)") +
  xlab("Ensembl ID") +   theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# resTblMane <- resTbl[-which(resTbl$annoMane == ""),]
# resTblMane$check <- factor(resTblMane$check, levels = c("yes", "no"))
# ggplot(resTblMane, aes(x = ENSID, y = conservationScore, fill = check)) + geom_boxplot() + 
#   geom_point(aes(fill = check), position = position_jitterdodge(), alpha = 1) + 
#   scale_fill_manual(values = c("darkgreen", "darkblue")) + ylab("MSA score (BLOSUM 100)") + 
#   xlab("Ensembl ID")


### second lollipop graphs


### lollipop plots will have to graph multiple graphs differently because y axis counts number of occurence
### the height should be 5 if count not available, but get rid of side labels

### need to create table of color palette for a domains

exampleDf <- rel_data_hg[which(rel_data_hg$accession == "K7EPC7"),]
exampleDf$order <- 1
exampleDf <- rbind(exampleDf, data.frame("type" = "LOLLIPOP", "description" = "5", "begin" = 200, "end" =  200, "length" = 1,
                                         "accession" = "K7EPC7", "entryName" = "K7EPC7_HUMAN", "taxid" = 9606, order = 1))
exampleDf <- rbind(exampleDf, data.frame("type" = "LOLLIPOP", "description" = "10", "begin" = 300, "end" =  300, "length" = 1,
                                         "accession" = "K7EPC7", "entryName" = "K7EPC7_HUMAN", "taxid" = 9606, order = 1))
exampleDf <- rbind(exampleDf, data.frame("type" = "LOLLIPOP", "description" = "200", "begin" = 310, "end" =  310, "length" = 1,
                                         "accession" = "K7EPC7", "entryName" = "K7EPC7_HUMAN", "taxid" = 9606, order = 1))
df <- exampleDf
draw_lolli <- function(df, outline_chain = "black", fill_chain = "grey", label_chains = TRUE,
                       label_size_chain = 4, point_color = NULL){
  require(ggplot2)
  require(ggfittext)
  require(gridExtra)
  require(datawizard)
  
  tmpDf <- df[df$type == "CHAIN", ]
  p <- ggplot() + ylim(c(0,1))
  p <- p + xlim(-max(df$end, na.rm = TRUE) * 0.2, 
                         max(df$end, na.rm = TRUE) + max(df$end, na.rm = TRUE) * 
                           0.1)
  p <- p + geom_rect(aes(xmin = tmpDf$begin, xmax = tmpDf$end, ymin = 0.3, ymax = 0.7),
                                     color = outline_chain, fill = fill_chain) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_blank(),
                       axis.ticks = element_blank(), axis.text = element_blank()) + xlab("") + ylab("")
  
  if (label_chains == TRUE) {
    p <- p + annotate("text", x = -1, y = 0.5, label = tmpDf$entryName,
                      hjust = 1, size = label_size_chain)
  }
  
  tmpDf3 <- df[df$type == "LOLLIPOP", ]
  ### adding one, since graph of yaxis of the lollipop technically is 1-based because of the domain figure and not 0-based
  ### for lollipop I'm going to add the count to description as character, then make it as numeric. easiest way without adding another column
  
  if (as.numeric(tmpDf3$description[1]) == 0) {
    tmpCount <- rep(as.numeric(tmpDf3$description[1]), nrow(tmpDf3))
  } else{
    tmpCount <- as.numeric(tmpDf3$description)
  }
  
  ### need so factor so that the protein chain + domains all stay the same size between different graphs
  ### the y-axis changes each depending on the data so it shrinks the protein graph portion
  tmpCount2 <- rescale(tmpCount, to= c(1,11), from = c(0,max(tmpCount)))
  if(!is.null(point_color)){
    tmpCount2_color <- point_color
  } else{
    tmpCount2_color <- rep("black", length(tmpCount2))
  }
  
  p <- p + geom_segment(aes(x = tmpDf3$begin, xend = tmpDf3$end, y = rep(0.7, length(tmpCount2)), yend = tmpCount2)) +
    geom_point(aes(x = tmpDf3$begin, y = tmpCount2), size = 3, color = tmpCount2_color) + 
    scale_y_continuous(limits=c(0, 11))
  
  
  tmpDf2 <- df[df$type == "DOMAIN", ]
  p <- p + geom_rect(mapping = aes(xmin = tmpDf2$begin, xmax = tmpDf2$end, ymin = 0.1, ymax = 0.9, fill = tmpDf2$description),
                     show.legend = FALSE)
  
  if (label_domains == TRUE) {
    p <- p + geom_fit_text(aes(xmin = tmpDf2$begin, xmax = tmpDf2$end, ymin = rep(0.21, length(tmpDf2$end)), ymax = rep(0.79, length(tmpDf2$end)), 
                               label = tmpDf2$description))
  }
  
  
  ### adding custom x and y axes, not using different grobs b/c this is inset in graph not bordering it
  p <- p + geom_segment(aes(x = 0, xend = 0, y = 0.9, yend = max(tmpCount2))) +
    geom_segment(aes(x = -max(df$end, na.rm = TRUE) * 0.01, xend = 0, y = 0.9, yend = 0.9)) + 
    geom_segment(aes(x = -max(df$end, na.rm = TRUE) * 0.01, xend = 0, y = max(tmpCount2), yend = max(tmpCount2))) +
    annotate("text", x = -max(df$end, na.rm = TRUE) * 0.01, y = max(tmpCount2), label = max(tmpCount), vjust = -1, hjust = 1) + 
    annotate("text", x = -max(df$end, na.rm = TRUE) * 0.01, y = 0.9, label = 0, vjust = -1, hjust = 2)
  
  
  p <- p + geom_segment(aes(x = 0, xend = max(tmpDf$end), y = 0, yend = 0)) +
    annotate("text", x = 0, y = 0, label = paste0(0, "aa"), hjust = 1.5) + 
    annotate("text", x = max(tmpDf$end), y = 0, label = paste0(max(tmpDf$end), "aa"), hjust = -0.5)
  
  
  return(p)
}


### looks like everything works so now time to build 1 graph for human and 3 for mouse - using EGFR as an example
###
###

interproDomaninsHg <- read.table("/mnt/DATA5/tmp/kev/interpro/20230301allHotspotGenes.fasta.tsv", header = FALSE,
                                 stringsAsFactors = FALSE, sep = "\t", quote="")


interproDomaninsMm <- read.table("/mnt/DATA5/tmp/kev/interpro/20230320allHotspotGenesMm10.fasta.tsv", header = FALSE,
                                 stringsAsFactors = FALSE, sep = "\t", quote="")


lollipopSearchHg <- getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id",
                                                        "ensembl_transcript_id_version",'uniprot_gn_id',
                                                        "uniprotswissprot", "uniprotsptrembl"),
                       filters = "ensembl_transcript_id", values = "ENST00000275493")

hgProtId <- hg38Peptide$uniprot_gn_id[grep("ENST00000275493", hg38Peptide$ensembl_transcript_id)]


### creating custom df to graph lollipops
drawProteins::get_features(paste(unlist(strsplit(hgProtId, ";")), collapse = " ")) -> rel_json_hg
drawProteins::feature_to_dataframe(rel_json_hg) -> rel_data_hg
rel_data_hg_filt <- rel_data_hg[grep("EGFR", rel_data_hg$entryName),]
rel_data_hg_filt <- rel_data_hg_filt[which(rel_data_hg_filt$type %in% "CHAIN"),]
rel_data_hg_filt$begin <- 1


annoTable2 <- read.table("/mnt/DATA5/tmp/kev/oncokbAnnoOut/mskTestGRCh38.maf", sep = "\t",
                         stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

# resTblMane_filt <- resTblMane[which(resTblMane$check == "yes"),]
resTblMane_filt <- resTblMane
varCountTbl <- egfrPositions[,c("Hugo_Symbol", "Protein_position")]


hgDomainMatch <- interproDomaninsHg[grep("ENST00000275493", interproDomaninsHg$V1),]

rel_data_hg_filt <- rbind(rel_data_hg_filt, data.frame("type" = rep("DOMAIN", nrow(hgDomainMatch)),
                                                       "description" = hgDomainMatch$V6,
                                                       "begin" = hgDomainMatch$V7, "end" = hgDomainMatch$V8,
                                                       "length" = hgDomainMatch$V8 - hgDomainMatch$V7,
                                                       "accession" = rep(rel_data_hg_filt$accession, nrow(hgDomainMatch)),
                                                       "entryName" = rep(rel_data_hg_filt$entryName, nrow(hgDomainMatch)),
                                                       "taxid" = rep(rel_data_hg_filt$taxid, nrow(hgDomainMatch)),
                                                       "order" = rep(rel_data_hg_filt$order, nrow(hgDomainMatch))))

hotspotCountTable <- table(egfrPositions$Protein_position)

for (i in 1:length(hotspotCountTable)) {
  rel_data_hg_filt <- rbind(rel_data_hg_filt, data.frame("type" = "LOLLIPOP", "description" = as.character(hotspotCountTable[i]),
                                                         "begin" = as.numeric(names(hotspotCountTable)[i]), "end" =  as.numeric(names(hotspotCountTable)[i]),
                                                         "length" = 1, "accession" = rel_data_hg_filt$accession[1], "entryName" = rel_data_hg_filt$entryName[1],
                                                         "taxid" = rel_data_hg_filt$taxid[1], order = rel_data_hg_filt$order[1]))
}


hgLolli <- draw_lolli(rel_data_hg_filt)


### for the mouse variants of the lollipop plot, I need to have mapped over amino acid coordinates
unique(resTblMane_filt$ENSID)

mmProtId <- mm10Peptide$uniprot_gn_id[grep(unique(resTblMane_filt$ENSID)[1], mm10Peptide$ensembl_transcript_id)]

mmDomainMatch <- interproDomaninsMm[grep(unique(resTblMane_filt$ENSID)[1], interproDomaninsMm$V1),]

### out of the three the first transcript matches best with Q01279
drawProteins::get_features("Q01279") -> rel_json_mm
drawProteins::feature_to_dataframe(rel_json_mm) -> rel_data_mm
rel_data_mm_filt <- rel_data_mm[which(rel_data_mm$type %in% "CHAIN"),]
rel_data_mm_filt$begin <- 1

rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = rep("DOMAIN", nrow(mmDomainMatch)),
                                                       "description" = mmDomainMatch$V6,
                                                       "begin" = mmDomainMatch$V7, "end" = mmDomainMatch$V8,
                                                       "length" = mmDomainMatch$V8 - mmDomainMatch$V7,
                                                       "accession" = rep(rel_data_mm_filt$accession, nrow(mmDomainMatch)),
                                                       "entryName" = rep(rel_data_mm_filt$entryName, nrow(mmDomainMatch)),
                                                       "taxid" = rep(rel_data_mm_filt$taxid, nrow(mmDomainMatch)),
                                                       "order" = rep(rel_data_mm_filt$order, nrow(mmDomainMatch))))

### need to map the counts from the human spots to mouse - interesting b/c the mouse mapped positions are different for each transcript
resTblMane_filt_mm <- resTblMane_filt[which(resTblMane_filt$ENSID == unique(resTblMane_filt$ENSID)[1]),
                                      c("conservationScore", "OriginalPosition", "ConvertedPosition", "check")]

hotspotCountTable_mm <- hotspotCountTable[which(names(hotspotCountTable) %in% resTblMane_filt_mm$OriginalPosition)]
names(hotspotCountTable_mm) <- resTblMane_filt_mm$ConvertedPosition[match(names(hotspotCountTable_mm), resTblMane_filt_mm$OriginalPosition)]

for (i in 1:length(hotspotCountTable_mm)) {
  rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = "LOLLIPOP", "description" = as.character(hotspotCountTable_mm[i]),
                                                         "begin" = as.numeric(names(hotspotCountTable_mm)[i]), "end" =  as.numeric(names(hotspotCountTable_mm)[i]),
                                                         "length" = 1, "accession" = rel_data_mm_filt$accession[1], "entryName" = rel_data_mm_filt$entryName[1],
                                                         "taxid" = rel_data_mm_filt$taxid[1], order = rel_data_mm_filt$order[1]))
}


colVector <- ifelse(resTblMane_filt_mm$check[match(names(hotspotCountTable_mm), resTblMane_filt_mm$ConvertedPosition)] == "yes", "darkgreen", "darkblue")
mmLolli_1 <- draw_lolli(rel_data_mm_filt, point_color = colVector)

### 2
###
###

mmDomainMatch <- interproDomaninsMm[grep(unique(resTblMane_filt$ENSID)[2], interproDomaninsMm$V1),]

### out of the three the first transcript matches best with Q5SVE7
drawProteins::get_features("Q9WVF5") -> rel_json_mm
drawProteins::feature_to_dataframe(rel_json_mm) -> rel_data_mm
rel_data_mm_filt <- rel_data_mm[which(rel_data_mm$type %in% "CHAIN"),]
rel_data_mm_filt$begin <- 1

rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = rep("DOMAIN", nrow(mmDomainMatch)),
                                                       "description" = mmDomainMatch$V6,
                                                       "begin" = mmDomainMatch$V7, "end" = mmDomainMatch$V8,
                                                       "length" = mmDomainMatch$V8 - mmDomainMatch$V7,
                                                       "accession" = rep(rel_data_mm_filt$accession, nrow(mmDomainMatch)),
                                                       "entryName" = rep(rel_data_mm_filt$entryName, nrow(mmDomainMatch)),
                                                       "taxid" = rep(rel_data_mm_filt$taxid, nrow(mmDomainMatch)),
                                                       "order" = rep(rel_data_mm_filt$order, nrow(mmDomainMatch))))

### need to map the counts from the human spots to mouse - interesting b/c the mouse mapped positions are different for each transcript
resTblMane_filt_mm <- resTblMane_filt[which(resTblMane_filt$ENSID == unique(resTblMane_filt$ENSID)[2]),
                                      c("conservationScore", "OriginalPosition", "ConvertedPosition", "check")]

hotspotCountTable_mm <- hotspotCountTable[which(names(hotspotCountTable) %in% resTblMane_filt_mm$OriginalPosition)]
names(hotspotCountTable_mm) <- resTblMane_filt_mm$ConvertedPosition[match(names(hotspotCountTable_mm), resTblMane_filt_mm$OriginalPosition)]

for (i in 1:length(hotspotCountTable_mm)) {
  rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = "LOLLIPOP", "description" = as.character(hotspotCountTable_mm[i]),
                                                         "begin" = as.numeric(names(hotspotCountTable_mm)[i]), "end" =  as.numeric(names(hotspotCountTable_mm)[i]),
                                                         "length" = 1, "accession" = rel_data_mm_filt$accession[1], "entryName" = rel_data_mm_filt$entryName[1],
                                                         "taxid" = rel_data_mm_filt$taxid[1], order = rel_data_mm_filt$order[1]))
}


colVector <- ifelse(resTblMane_filt_mm$check[match(names(hotspotCountTable_mm), resTblMane_filt_mm$ConvertedPosition)] == "yes", "darkgreen", "darkblue")
mmLolli2 <- draw_lolli(rel_data_mm_filt, point_color = colVector)

### 3
###
###

mmDomainMatch <- interproDomaninsMm[grep(unique(resTblMane_filt$ENSID)[3], interproDomaninsMm$V1),]

### out of the three the first transcript matches best with Q5SVE7
drawProteins::get_features("Q5SVE7") -> rel_json_mm
drawProteins::feature_to_dataframe(rel_json_mm) -> rel_data_mm
rel_data_mm_filt <- rel_data_mm[which(rel_data_mm$type %in% "CHAIN"),]
rel_data_mm_filt$begin <- 1

rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = rep("DOMAIN", nrow(mmDomainMatch)),
                                                       "description" = mmDomainMatch$V6,
                                                       "begin" = mmDomainMatch$V7, "end" = mmDomainMatch$V8,
                                                       "length" = mmDomainMatch$V8 - mmDomainMatch$V7,
                                                       "accession" = rep(rel_data_mm_filt$accession, nrow(mmDomainMatch)),
                                                       "entryName" = rep(rel_data_mm_filt$entryName, nrow(mmDomainMatch)),
                                                       "taxid" = rep(rel_data_mm_filt$taxid, nrow(mmDomainMatch)),
                                                       "order" = rep(rel_data_mm_filt$order, nrow(mmDomainMatch))))

### need to map the counts from the human spots to mouse - interesting b/c the mouse mapped positions are different for each transcript
resTblMane_filt_mm <- resTblMane_filt[which(resTblMane_filt$ENSID == unique(resTblMane_filt$ENSID)[3]),
                                      c("conservationScore", "OriginalPosition", "ConvertedPosition", "check")]

hotspotCountTable_mm <- hotspotCountTable[which(names(hotspotCountTable) %in% resTblMane_filt_mm$OriginalPosition)]
names(hotspotCountTable_mm) <- resTblMane_filt_mm$ConvertedPosition[match(names(hotspotCountTable_mm), resTblMane_filt_mm$OriginalPosition)]

for (i in 1:length(hotspotCountTable_mm)) {
  rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = "LOLLIPOP", "description" = as.character(hotspotCountTable_mm[i]),
                                                         "begin" = as.numeric(names(hotspotCountTable_mm)[i]), "end" =  as.numeric(names(hotspotCountTable_mm)[i]),
                                                         "length" = 1, "accession" = rel_data_mm_filt$accession[1], "entryName" = rel_data_mm_filt$entryName[1],
                                                         "taxid" = rel_data_mm_filt$taxid[1], order = rel_data_mm_filt$order[1]))
}


colVector <- ifelse(resTblMane_filt_mm$check[match(names(hotspotCountTable_mm), resTblMane_filt_mm$ConvertedPosition)] == "yes", "darkgreen", "darkblue")
mmLolli3 <- draw_lolli(rel_data_mm_filt, point_color = colVector)

### create the 4 individual plots + boxplot for the figure, on the figure have proportions that matched vs unmatched
### I need to standardize colors .... for now I'll do it in illustrator

### this type of mapping solves a majority of the process ..... could always see if highest number of maps is always the homologue?
### for egfr: use mane transcript, then homologue peptide into homologue transcript id
# getBM(mart = ensemblHg, attributes = c("external_gene_name", "ensembl_transcript_id",
#                                        "mmusculus_homolog_ensembl_peptide"),
#       filters = "ensembl_transcript_id", values = "ENST00000275493")
# 
# 
# getBM(mart = ensemblMm, attributes = c("external_gene_name", "ensembl_transcript_id",
#                                        "ensembl_peptide_id"),
#       filters = "ensembl_peptide_id", values = "ENSMUSP00000020329")

tsgOncoList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/volgelstein2012oncoTsgList.xlsx", sheet = 6, skip = 1)
tsgList <- tsgOncoList$`Gene Symbol`[which(tsgOncoList$`Classification*` == "TSG")]
oncoList <- tsgOncoList$`Gene Symbol`[which(tsgOncoList$`Classification*` == "Oncogene")]

tmpGraphTbl <- resTblMane
tmpGraphTbl$tmp <- paste(resTblMane$check, resTblMane$MutEff)
ggboxplot(tmpGraphTbl, y = "conservationScore", x = "check")  + stat_compare_means(method = "t.test") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tmpGraphTbl$gene <- geneNameDf$external_gene_name[match(tmpGraphTbl$originalEns, geneNameDf$ensembl_transcript_id)]
tmpGraphTbl <- tmpGraphTbl[order(tmpGraphTbl$gene), ]
### quick script to check for lying in domain for both human and mouse

domainMatchMm <- NULL
domainMatchHg <- NULL
i <- unique(tmpGraphTbl$gene)[5]
for (i in unique(tmpGraphTbl$gene)) {
  print(i)
  # for each protein get domain information for human and mouse
  tmpTable <- tmpGraphTbl[which(tmpGraphTbl$gene == i),]
  tmpMmEnsId <- unique(tmpTable$ENSID)
  tmpMmProtId <- uniprotWsTableMm10$Entry[which(uniprotWsTableMm10$From == tmpMmEnsId)]
  tmpHgEnsId <- unique(tmpTable$originalEns)
  tmpHgProtId <- uniprotWsTableHg38$Entry[which(uniprotWsTableHg38$From == tmpHgEnsId)]
  
  tmpMmDomainMatch <- interproDomaninsMm[grep(tmpMmEnsId, interproDomaninsMm$V1),]
  tmpHgDomainMatch <- interproDomaninsHg[grep(tmpHgEnsId, interproDomaninsHg$V1),]
  
  ### dummy chromosome since these are amino acid positions and I'm comparing protein by protein
  tmpMmProtGrange <- GRanges(seqnames = rep(1, nrow(tmpMmDomainMatch)),
                             IRanges(start = tmpMmDomainMatch$V7, end = tmpMmDomainMatch$V8))
  tmpHgProtGrange <- GRanges(seqnames = rep(1, nrow(tmpHgDomainMatch)),
                             IRanges(start = tmpHgDomainMatch$V7, end = tmpHgDomainMatch$V8))
  
  tmpMmHotspotGrange <- GRanges(seqnames = rep(1, nrow(tmpTable)),
                                IRanges(start = tmpTable$ConvertedPosition, end = tmpTable$ConvertedPosition))
  tmpHgHotspotGrange <- GRanges(seqnames = rep(1, nrow(tmpTable)),
                                IRanges(start = tmpTable$OriginalPosition, end = tmpTable$OriginalPosition))
  
  tmpResMm <- rep("no", nrow(tmpTable))
  tmpResHg <- rep("no", nrow(tmpTable))
  
  tmpResMm[subjectHits(findOverlaps(tmpMmProtGrange, tmpMmHotspotGrange))] <- "yes"
  tmpResHg[subjectHits(findOverlaps(tmpHgProtGrange, tmpHgHotspotGrange))] <- "yes"
  
  domainMatchMm <- c(domainMatchMm , tmpResMm)
  domainMatchHg <- c(domainMatchHg, tmpResHg)
}

tmpGraphTbl$HgDomainMatch <- domainMatchHg
tmpGraphTbl$MmDomainMatch <- domainMatchMm

tmpGraphTbl2 <- tmpGraphTbl[which(tmpGraphTbl$MutEff %in% c("Unknown", "Loss-of-function",
                                                            "Gain-of-function", "Likely Gain-of-function",
                                                            "Likely Loss-of-function")),]

tmpGraphTbl2$geneDesig <- NA
tmpGraphTbl2$geneDesig[which(tmpGraphTbl2$gene %in% tsgInTable)] <- "tsg"
tmpGraphTbl2$geneDesig[-which(tmpGraphTbl2$gene %in% tsgInTable)] <- "onco"

tmpGraphTbl2$domainMatch <- "mismatch"
tmpGraphTbl2$domainMatch[which(tmpGraphTbl2$HgDomainMatch == "no" & tmpGraphTbl2$MmDomainMatch == "no")]  <- "no"
tmpGraphTbl2$domainMatch[which(tmpGraphTbl2$HgDomainMatch == "yes" & tmpGraphTbl2$MmDomainMatch == "yes")]  <- "yes"

### dummy coding muteff unknown is baseline vs known functions, domain match baseline will be no, and tsgs will be baseline

tmpGraphTbl2$dummyMutEff <- tmpGraphTbl2$MutEff
tmpGraphTbl2$dummyMutEff <- str_replace_all(tmpGraphTbl2$dummyMutEff, "Unknown", "0")
tmpGraphTbl2$dummyMutEff <- str_replace_all(tmpGraphTbl2$dummyMutEff, "Likely Gain-of-function", "1")
tmpGraphTbl2$dummyMutEff <- str_replace_all(tmpGraphTbl2$dummyMutEff, "Likely Loss-of-function", "2")
tmpGraphTbl2$dummyMutEff <- str_replace_all(tmpGraphTbl2$dummyMutEff, "Gain-of-function", "3")
tmpGraphTbl2$dummyMutEff <- str_replace_all(tmpGraphTbl2$dummyMutEff, "Loss-of-function", "4")
tmpGraphTbl2$dummyMutEff <- factor(tmpGraphTbl2$dummyMutEff)

tmpGraphTbl2$dummyDomain <- tmpGraphTbl2$domainMatch
tmpGraphTbl2$dummyDomain <- str_replace_all(tmpGraphTbl2$dummyDomain, "no", "0")
tmpGraphTbl2$dummyDomain <- str_replace_all(tmpGraphTbl2$dummyDomain, "yes", "1")
tmpGraphTbl2$dummyDomain <- str_replace_all(tmpGraphTbl2$dummyDomain, "mismatch", "2")
tmpGraphTbl2$dummyDomain <- factor(tmpGraphTbl2$dummyDomain, levels = c(2, 0, 1))

tmpGraphTbl2$dummyDesig <- tmpGraphTbl2$geneDesig
tmpGraphTbl2$dummyDesig <- str_replace_all(tmpGraphTbl2$dummyDesig, "onco", "0")
tmpGraphTbl2$dummyDesig <- str_replace_all(tmpGraphTbl2$dummyDesig, "tsg", "1")
tmpGraphTbl2$dummyDesig <- factor(tmpGraphTbl2$dummyDesig)

### second cutoff was needs matching and also greater than certain conservation score - modified z score cause of skewness
# tmpGraphTbl2$check2 <- ifelse(tmpGraphTbl2$conservationScore > 15.8 & tmpGraphTbl2$check == "yes", "yes", "no")

tmpGraphTbl2$check2 <- tmpGraphTbl2$check

### swapped baseline from no, so OR increases make more sense 

tmpGraphTbl2$check <- factor(tmpGraphTbl2$check, levels = c("no", "yes"))
tmpGraphTbl2$check2 <- factor(tmpGraphTbl2$check2, levels = c("no", "yes"))

oncokbPositionsAllHg38$string <- paste0(oncokbPositionsAllHg38$Hugo_Symbol, oncokbPositionsAllHg38$Protein_position, oncokbPositionsAllHg38$AminoAcid)
oncokbPositionsRed$string <- paste0(oncokbPositionsRed$Hugo_Symbol, oncokbPositionsRed$Protein_position, oncokbPositionsRed$AminoAcid)
oncokbPositionsRed$count <- 0

for (i in seq_along(oncokbPositionsRed$string)) {
  oncokbPositionsRed$count[i] <- length(which(oncokbPositionsAllHg38$string == oncokbPositionsRed$string[i]))
}

tmpGraphTbl2$count <- oncokbPositionsRed$count[match(tmpGraphTbl2$Hotspot, oncokbPositionsRed$string)]

# glmRes <- glm(data = tmpGraphTbl2, formula = check ~ conservationScore + dummyMutEff + dummyDomain + dummyDesig, family = binomial(link = "logit"))
glmRes <- glm(data = tmpGraphTbl2, formula = check2 ~ conservationScore + dummyMutEff + dummyDesig, family = binomial(link = "logit"))
glmRes <- glm(data = tmpGraphTbl2, formula = check ~ conservationScore + dummyDesig, family = binomial(link = "logit"))
glmRes <- glm(data = tmpGraphTbl2, formula = check ~ dummyDesig, family = binomial(link = "logit"))

summary(glmRes)

effectsize::standardize_parameters(glmRes, exponentiate = TRUE, two_sd = TRUE, robust = TRUE)

### will need this to plugin into forest_model
glmResRefit <- effectsize::standardize(glmRes, two_sd = TRUE, robust = TRUE)
parameters::model_parameters(glmResRefit, exponentiate = TRUE)

library(forestmodel)
# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231205forestPlot.pdf", width = 12, height = 8, useDingbats = FALSE)
forest_model(glmResRefit, exponentiate = TRUE)
dev.off()
### interesting results from looking at both p-values and effect size for whether hotspots converted properly
### out of all the variables, how well the local area is conserved increases the likely hood of conversion the most i.e 4.36
### the next highest surprisingly is if the the mutation is gain of function with 2.78
### interpreting the confidence intervals showing true mean of both conservation score and being a gain of function mutation
### does always increase the chance of whether the hotspot is conserved between the two species
### at least more likely than loss of function which is interesting. both likely and loss of function 95% CI range
### includes both >1 and < 1


### now for the figures, probably just boxplots dividing data by type of mutations using conservation score as boxplot
### too many dots for jitter, use fill by check show diff
### only remaking three to reorder check, not wanting to change reference for multiple regression previously done
tmpGraphTbl3 <- tmpGraphTbl2
tmpGraphTbl3$check <- factor(tmpGraphTbl3$check, levels = c("yes", "no"))

ggboxplot(tmpGraphTbl3, y = "conservationScore", x = "MutEff", outlier.shape=NA, fill = "check") +
  scale_fill_manual(values=c("#8B0000", "#00008B")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))


ggboxplot(tmpGraphTbl3, y = "conservationScore", x = "domainMatch", outlier.shape=NA, fill = "check") +
  scale_fill_manual(values=c("#8B0000", "#00008B")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))


ggboxplot(tmpGraphTbl3, y = "conservationScore", x = "geneDesig", outlier.shape=NA, fill = "check") +
  scale_fill_manual(values=c("#8B0000", "#00008B")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))


# load("/mnt/DATA5/tmp/kev/misc/20230402greedyHotspot.RData")
# save.image(file = "/mnt/DATA5/tmp/kev/misc/20230402greedyHotspot.RData")
### I think boxplots splitting it by function for regular graph
### split it by oncogenes/tsg and whether it lies in a domain in supplementary
### have p-values for these too




moreThan4 <- names(table(tmpGraphTbl3$gene))[which(table(tmpGraphTbl3$gene) > 4)]
tmpGraphTblFilt <- tmpGraphTbl3[which(tmpGraphTbl3$gene %in%moreThan4),]
tsgInTable <- unique(tmpGraphTblFilt$gene)[which(unique(tmpGraphTblFilt$gene) %in% tsgList)]
tmpGraphTblFilt_tsg <- tmpGraphTblFilt[which(tmpGraphTblFilt$gene %in% tsgInTable),]
byMedTsg <- with(tmpGraphTblFilt_tsg, reorder(gene, conservationScore, median))
tmpGraphTblFilt_onco <- tmpGraphTblFilt[-which(tmpGraphTblFilt$gene %in% tsgInTable),]
byMedOnco <- with(tmpGraphTblFilt_onco, reorder(gene, conservationScore, median))

customLevels <- c(rev(levels(byMedTsg)), rev(levels(byMedOnco)))
tmpGraphTblFilt$gene <- factor(tmpGraphTblFilt$gene, levels = customLevels)

tmpGraphTblFilt$color <- ifelse(tmpGraphTblFilt$check == "yes", "darkblue", "darkred")

# tmpGraphTblFilt <- tmpGraphTblFilt[order(tmpGraphTblFilt$gene),]
ggboxplot(tmpGraphTblFilt, y = "conservationScore", x = "gene", outlier.shape=NA)  + geom_jitter(color = tmpGraphTblFilt$color, size = 0.5, alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))



### drawing the two plots for notch1 and tp53
###
###

notch1Positions <- oncokbPositions[which(oncokbPositions$Hugo_Symbol == "NOTCH3"),]

hotspotCountTable <- table(notch1Positions$Protein_position)

mmNotchEnsId <- unique(tmpGraphTblFilt$ENSID[which(tmpGraphTblFilt$gene == "NOTCH3")])

mmProtId <- uniprotWsTableMm10$Entry[which(uniprotWsTableMm10$From == mmNotchEnsId)]
# mmProtId <- mm10Peptide$uniprot_gn_id[grep(mmNotchEnsId, mm10Peptide$ensembl_transcript_id)]

mmDomainMatch <- interproDomaninsMm[grep(mmNotchEnsId, interproDomaninsMm$V1),]

### out of the three the first transcript matches best with Q01279
drawProteins::get_features(mmProtId) -> rel_json_mm
drawProteins::feature_to_dataframe(rel_json_mm) -> rel_data_mm
rel_data_mm_filt <- rel_data_mm[which(rel_data_mm$type %in% "CHAIN"),]
rel_data_mm_filt$begin <- 1

rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = rep("DOMAIN", nrow(mmDomainMatch)),
                                                       "description" = mmDomainMatch$V6,
                                                       "begin" = mmDomainMatch$V7, "end" = mmDomainMatch$V8,
                                                       "length" = mmDomainMatch$V8 - mmDomainMatch$V7,
                                                       "accession" = rep(rel_data_mm_filt$accession, nrow(mmDomainMatch)),
                                                       "entryName" = rep(rel_data_mm_filt$entryName, nrow(mmDomainMatch)),
                                                       "taxid" = rep(rel_data_mm_filt$taxid, nrow(mmDomainMatch)),
                                                       "order" = rep(rel_data_mm_filt$order, nrow(mmDomainMatch))))

### need to map the counts from the human spots to mouse - interesting b/c the mouse mapped positions are different for each transcript
resTblMane_filt_mm <- tmpGraphTblFilt[which(tmpGraphTblFilt$ENSID == mmNotchEnsId),
                                      c("conservationScore", "OriginalPosition", "ConvertedPosition", "check")]

hotspotCountTable_mm <- hotspotCountTable[which(names(hotspotCountTable) %in% resTblMane_filt_mm$OriginalPosition)]
names(hotspotCountTable_mm) <- resTblMane_filt_mm$ConvertedPosition[match(names(hotspotCountTable_mm), resTblMane_filt_mm$OriginalPosition)]

for (i in 1:length(hotspotCountTable_mm)) {
  rel_data_mm_filt <- rbind(rel_data_mm_filt, data.frame("type" = "LOLLIPOP", "description" = as.character(hotspotCountTable_mm[i]),
                                                         "begin" = as.numeric(names(hotspotCountTable_mm)[i]), "end" =  as.numeric(names(hotspotCountTable_mm)[i]),
                                                         "length" = 1, "accession" = rel_data_mm_filt$accession[1], "entryName" = rel_data_mm_filt$entryName[1],
                                                         "taxid" = rel_data_mm_filt$taxid[1], order = rel_data_mm_filt$order[1]))
}

colVector <- ifelse(resTblMane_filt_mm$check[match(names(hotspotCountTable_mm), resTblMane_filt_mm$ConvertedPosition)] == "yes", "darkblue", "darkred")
notchLolli <- draw_lolli(rel_data_mm_filt, point_color = colVector)
notchLolli



### drawing newer plots with conservation score (y-axis) for PIK3CA - 
### need to have distinguish between annotations  somehow
### draw it on the human protein
# df <- rel_data_hg_filt2
draw_lolliV2 <- function(df, outline_chain = "black", fill_chain = "grey", label_chains = TRUE,
                         label_size_chain = 4, point_color = NULL, label_domains = FALSE){
  ### point color is vector with legnth of number of points - based on matching
  
  
  ### this version I should make it so y axis is conservation score of each position
  
  require(ggplot2)
  require(ggfittext)
  require(gridExtra)
  require(datawizard)
  
  tmpDf <- df[df$type == "CHAIN", ]
  p <- ggplot() + ylim(c(0,1))
  p <- p + xlim(-max(df$end, na.rm = TRUE) * 0.2, 
                max(df$end, na.rm = TRUE) + max(df$end, na.rm = TRUE) * 
                  0.1)
  p <- p + geom_rect(aes(xmin = tmpDf$begin, xmax = tmpDf$end, ymin = 0.3, ymax = 0.7),
                     color = outline_chain, fill = fill_chain) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_blank(),
                       axis.ticks = element_blank(), axis.text = element_blank()) + xlab("") + ylab("")
  
  if (label_chains == TRUE) {
    p <- p + annotate("text", x = -1, y = 0.5, label = tmpDf$entryName,
                      hjust = 1, size = label_size_chain)
  }
  
  tmpDf3 <- df[df$type == "LOLLIPOP", ]
  ### adding one, since graph of yaxis of the lollipop technically is 1-based because of the domain figure and not 0-based
  ### for lollipop I'm going to add the count to description as character, then make it as numeric. easiest way without adding another column
  
  if (as.numeric(tmpDf3$description[1]) == 0) {
    tmpCount <- rep(as.numeric(tmpDf3$description[1]), nrow(tmpDf3))
  } else{
    tmpCount <- as.numeric(tmpDf3$description)
  }
  
  ### need so factor so that the protein chain + domains all stay the same size between different graphs
  ### the y-axis changes each depending on the data so it shrinks the protein graph portion
  tmpCount2 <- rescale(tmpCount, to= c(1,11), from = c(0,max(tmpCount)))
  if(!is.null(point_color)){
    tmpCount2_color <- point_color
  } else{
    tmpCount2_color <- rep("black", length(tmpCount2))
  }
  tmpDf3$taxid <- factor(tmpDf3$taxid)
  tmpDf3$taxid <- relevel(tmpDf3$taxid, ref =c("Unknown"))
  
  p <- p + geom_segment(aes(x = tmpDf3$begin, xend = tmpDf3$end, y = rep(0.7, length(tmpCount2)), yend = tmpCount2),
                        color = rep("grey", length(tmpDf3$begin)), alpha = 0.3) +
    geom_point(aes(x = tmpDf3$begin, y = tmpCount2, shape = tmpDf3$taxid), alpha = 0.7,size = 3, color = tmpCount2_color) + 
    scale_y_continuous(limits=c(0, 11))
  
  
  tmpDf2 <- df[df$type == "DOMAIN", ]
  p <- p + geom_rect(mapping = aes(xmin = tmpDf2$begin, xmax = tmpDf2$end, ymin = 0.1, ymax = 0.9, fill = tmpDf2$description),
                     show.legend = FALSE)
  
  if (label_domains == TRUE) {
    p <- p + geom_fit_text(aes(xmin = tmpDf2$begin, xmax = tmpDf2$end, ymin = rep(0.21, length(tmpDf2$end)), ymax = rep(0.79, length(tmpDf2$end)), 
                               label = tmpDf2$description))
  }
  
  
  ### adding custom x and y axes, not using different grobs b/c this is inset in graph not bordering it
  p <- p + geom_segment(aes(x = 0, xend = 0, y = 0.9, yend = max(tmpCount2))) +
    geom_segment(aes(x = -max(df$end, na.rm = TRUE) * 0.01, xend = 0, y = 0.9, yend = 0.9)) + 
    geom_segment(aes(x = -max(df$end, na.rm = TRUE) * 0.01, xend = 0, y = max(tmpCount2), yend = max(tmpCount2))) +
    annotate("text", x = -max(df$end, na.rm = TRUE) * 0.01, y = max(tmpCount2), label = max(tmpCount), vjust = -1, hjust = 1) + 
    annotate("text", x = -max(df$end, na.rm = TRUE) * 0.01, y = 0.9, label = 0, vjust = -1, hjust = 2)
  
  
  p <- p + geom_segment(aes(x = 0, xend = max(tmpDf$end), y = 0, yend = 0)) +
    annotate("text", x = 0, y = 0, label = paste0(0, "aa"), hjust = 1.5) + 
    annotate("text", x = max(tmpDf$end), y = 0, label = paste0(max(tmpDf$end), "aa"), hjust = -0.5)
  
  
  return(p)
}



pik3caPositions <- oncokbPositions[which(oncokbPositions$Hugo_Symbol == "EGFR"),]
hgPik3caEnsId <- unique(tmpGraphTblFilt$originalEns[which(tmpGraphTblFilt$gene == "EGFR")])
hgProtId <- uniprotWsTableHg38$Entry[which(uniprotWsTableHg38$From == hgPik3caEnsId)]
hgDomainMatch <- interproDomaninsHg[grep(hgPik3caEnsId, interproDomaninsHg$V1),]
drawProteins::get_features(hgProtId) -> rel_json_hg
drawProteins::feature_to_dataframe(rel_json_hg) -> rel_data_hg
rel_data_hg_filt <- rel_data_hg[which(rel_data_hg$type %in% "CHAIN"),]
rel_data_hg_filt$begin <- 1
rel_data_hg_filt <- rbind(rel_data_hg_filt, data.frame("type" = rep("DOMAIN", nrow(hgDomainMatch)),
                                                       "description" = hgDomainMatch$V6,
                                                       "begin" = hgDomainMatch$V7, "end" = hgDomainMatch$V8,
                                                       "length" = hgDomainMatch$V8 - hgDomainMatch$V7,
                                                       "accession" = rep(rel_data_hg_filt$accession, nrow(hgDomainMatch)),
                                                       "entryName" = rep(rel_data_hg_filt$entryName, nrow(hgDomainMatch)),
                                                       "taxid" = rep(rel_data_hg_filt$taxid, nrow(hgDomainMatch)),
                                                       "order" = rep(rel_data_hg_filt$order, nrow(hgDomainMatch))))


resTblMane_filt_hg <- tmpGraphTblFilt[which(tmpGraphTblFilt$originalEns == hgPik3caEnsId),
                                      c("conservationScore", "OriginalPosition", "ConvertedPosition", "check2", "MutEff")]
 
rel_data_hg_filt2 <- rel_data_hg_filt
for (i in 1:nrow(resTblMane_filt_hg)) {
  rel_data_hg_filt2 <- rbind(rel_data_hg_filt2, data.frame("type" = "LOLLIPOP", "description" = as.character(resTblMane_filt_hg$conservationScore[i]),
                                                         "begin" = as.numeric(resTblMane_filt_hg$OriginalPosition[i]),
                                                         "end" =  as.numeric(resTblMane_filt_hg$OriginalPosition[i]),
                                                         "length" = 1, "accession" = rel_data_hg_filt2$accession[1],
                                                         "entryName" = rel_data_hg_filt2$entryName[1],
                                                         "taxid" = resTblMane_filt_hg$MutEff[i],
                                                         order = rel_data_hg_filt2$order[1]))
}

colVector <- ifelse(resTblMane_filt_hg$check2 == "yes", "darkred", "darkblue")
pik3caLolli <- draw_lolliV2(rel_data_hg_filt2, point_color = colVector)
pik3caLolli

pdf("/mnt/DATA5/tmp/kev/misc/20231207newLollipop.pdf", useDingbats = FALSE, width = 12, height = 6)
pik3caLolli
dev.off()

### need to redo the proportion of conversion results, could be yes, but spurious
### look at  violin plot to see if there is a cutoff for good vs spurious
### data is extremely negatively skewed; so z-scores would be useless; use modified z-score 3.5 cutoff

quantile(tmpGraphTbl3$conservationScore, seq(0, 1, 0.01))

median(tmpGraphTbl3$conservationScore)
mad(tmpGraphTbl3$conservationScore)

tmpGraphTbl3$modZScore <- (0.6745 * (tmpGraphTbl3$conservationScore - median(tmpGraphTbl3$conservationScore))/mad(tmpGraphTbl3$conservationScore))
max(tmpGraphTbl3$conservationScore[which(tmpGraphTbl3$modZScore < -3.5)])
### 15.8 minimum - with the new minimum, I need to redo forest plot

### number of correct conversions based on the second check and cutoff 
nrow(tmpGraphTbl2)
table(tmpGraphTbl2$check2)

