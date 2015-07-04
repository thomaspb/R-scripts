rm(list=ls())



# NOTE. NOT BOTHERING ABOUT NORMALIZING HERE. JUST WANT TO FIND TAGS, NOT COMPARE SAMPLES!

# work
setwd("~/Dropbox/ZERB 2015/CAGE")

# new: from FANTOM website Feb 2015
# Expression (read counts) of robust phase 1 CAGE peaks for human samples with annotation
#@ hg19.cage_peak_counts_ann_decoded.osc.txt.gz.extract
# 1GB file downloaded from RIKEN. Takes a long time to import, so saved as R object afterwards!
#@ CAGEinput <- read.delim("hg19.cage_peak_counts_ann_decoded.osc.tsv", comment.char = "#", header=T)
#@ head(CAGEinput[1,1:5])
#@ head(CAGEinput[,1:1])
#@ View(CAGEinput[1,1:12])
#@ save(CAGEinput,file = "CAGEinput-RAW.Robj") 
load("CAGEinput-RAW.Robj") 
CAGEinput2 <- CAGEinput # create backup data frame


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What does CAGEinput looks like?
colnames(CAGEinput2[1]) # X00Annotation
colnames(CAGEinput2[2]) # short_description
colnames(CAGEinput2[3]) # description
colnames(CAGEinput2[4]) # association_with_transcript
colnames(CAGEinput2[5]) # entrezgene_id
colnames(CAGEinput2[6]) # hgnc_id
colnames(CAGEinput2[7]) # uniprot_id
# CAGEinput2 <- with(CAGEinput2, CAGEinput2[!(hgnc_id  == "" | is.na(hgnc_id)), ]) # keep all here!
#
colnames(CAGEinput2[8]) # first sample
head(CAGEinput2[1]) # e.g. chr10:100013403..100013414,-
head(CAGEinput2)

# get rid of columns 2 to 7
CAGEinputX <- CAGEinput2[ , -which(names(CAGEinput2) %in% c("short_description", "short_description", "description", "association_with_transcript", "entrezgene_id", "hgnc_id", "uniprot_id" ))] ##

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# keep only samples of interest + column 1 = CAGE tag location/X00Annotation
CAGEinputX <- CAGEinputX[ , which(names(CAGEinputX) %in% c("X00Annotation", "MCF7.breast.cancer.cell.line.response.to.EGF1..00hr00min..biol_rep2.CNhs12475.13097.140D1", "MCF7.breast.cancer.cell.line.response.to.EGF1..00hr00min..biol_rep3.CNhs12703.13163.141B4","breast.carcinoma.cell.line.MCF7.CNhs11943.10482.107A5", "breast.carcinoma.cell.line.MDA.MB.453.CNhs10736.10419.106C5", "breast..adult..donor1.CNhs11792.10080.102A8" ))] ## works as expected
#
# NOTE: still do not have rownames! should assign the desired locations as rownames! 
# 
# works, now only have the 5 breast-derived samples of interest 
rownames(CAGEinputX) <- CAGEinputX[,1]
# CAGEinputX <- CAGEinputX[,2:6] # get rid of first column
head(CAGEinputX)
colnames(CAGEinputX)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOW READY TO EXTRACT MATCHES

# choose on of the sub-sections below (corresponds to gene loci of interest)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 22.03.2015: Any ZNF440 tags?
# ZNF440
# hg19 chr19:11925107-11946016  AND +   
TargetChromosome <- "chr19"
TargetStart <- as.numeric("11925107")
TargetStop <-  as.numeric("11946016")
TargetStrand <- ",+"  # include a , (my regex was not perfect)
CAGEinputX <- CAGEinputX[grep("^chr19",rownames(CAGEinputX)),]
# chr19:11925071..11925113,+  p1@ZNF440	CAGE_peak_1_at_ZNF440_5end	11925*071*	11925*113*	chr19	+
TargetStart <- as.numeric(TargetStart - 500) # usually one will paste gene coordinates above, but promoters are 5'
TargetStop <-  as.numeric(TargetStop + 500)
CAGEinputX["chr19:11925071..11925113,+",]   # ZNF440 is there!
which(grepl("chr19:11925071..11925113,+", rownames(CAGEinputX))) # what row #?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 22.03.2015: internal ZNF440 tags. NOT interested in canonical ZNF440 TSS!
TargetChromosome <- "chr19"
TargetStart <- as.numeric("11925107")
TargetStop <-  as.numeric("11946016")
TargetStrand <- ",+"  # include a , (my regex was not perfect)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GHRLOS
# e.g. chr3:10326928..10326933,+
TargetChromosome <- "chr3"
TargetStart <- as.numeric("10327432")
TargetStop <-  as.numeric("10335133")
TargetStart <- as.numeric(TargetStart - 500) # usually one will paste gene coordinates above, but promoters are 5'
TargetStop <-  as.numeric(TargetStop + 500)
TargetStrand <- ",+"  # include a , (my regex was not perfect)
CAGEinputX <- CAGEinputX[grep("^chr3",rownames(CAGEinputX)),]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GHRL
# e.g. chr3:10331943..10331961,-
#want to include exon -1: chr3:10331943..10331961,-  p5@GHRL	CAGE_peak_5_at_GHRL_5end	10331943	10331961 chr3	-
TargetChromosome <- "chr3"
TargetStart <- as.numeric("10327434")
TargetStop <-  as.numeric("10334631")
TargetStart <- as.numeric(TargetStart - 5000) # usually one will paste gene coordinates above, but promoters are 5'
TargetStop <-  as.numeric(TargetStop + 5000)
TargetStrand <- ",-"  # include a , (my regex was not perfect)
CAGEinputX <- CAGEinputX[grep("^chr3",rownames(CAGEinputX)),]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GHSR
TargetChromosome <- "chr3"
TargetStart <- as.numeric("172161081")
TargetStop <-  as.numeric("172166246")
TargetStart <- as.numeric(TargetStart - 500) # usually one will paste gene coordinates above, but promoters are 5'
TargetStop <-  as.numeric(TargetStop + 500)
TargetStrand <- ",-"  # include a , (my regex was not perfect)
CAGEinputX <- CAGEinputX[grep("^chr3",rownames(CAGEinputX)),]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GHSROS --looking for internal CAGE start site only here!
TargetChromosome <- "chr3"
TargetStart <- as.numeric("172161081")
TargetStop <-  as.numeric("172166246")
TargetStrand <- ",+"  # include a , (my regex was not perfect)
CAGEinputX <- CAGEinputX[grep("^chr3",rownames(CAGEinputX)),]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("./CAGEtagSearch.FUN.R")
CAGE.resultDF <- do.call(rbind.data.frame, CAGE.result)  # convert the list to a data frame!
head(CAGE.resultDF)

# write table 
save(CAGE.resultDF,file = "CAGE.resultDF.Robj") 
write.table(CAGE.resultDF, "CAGEtag-output.csv", row.names=T, col.names=T, sep="\t", quote=F)

# only keep samples that have counts
CAGE.resultDF2 <- CAGE.resultDF
CAGE.results.Counts <- CAGE.resultDF2[, !apply(CAGE.resultDF2 == 0, 2, all)]
head(CAGE.results.Counts)
write.table(CAGE.results.Counts, "CAGEtag-output-CountsOnly.csv", row.names=T, col.names=T, sep="\t", quote=F)
