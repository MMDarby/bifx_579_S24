##################################
# creates GRanges list of all reads in .bam that overlap with genomic region
##################################

library(GenomicRanges)
library(ggbio)

######## things to change ############
# 3 things to change: bamDir, setwd(), target_ranges

#Set directories 
# directory with .bam files
#bamDir <- "/net/fs01/mnt/storage/home/mdarby/indexedMouse1BamFiles"
bamDir <- "~/Dropbox/RepeatAnalysis/bamFiles"
#working directory
#setwd("/net/fs01/mnt/storage/home/mdarby/repeatAnalysis/function")
setwd("~/Dropbox/RepeatAnalysis/plotCounts/function")

#Set list of ranges of interest MUST BE STRAND SPECIFIC and match repeat range FOR COUNTS
#load("~/Dropbox/RepeatAnalysis/rmsk_files/reduced_fDbRepeats.rda")
#assign seqlengths
#library(BSgenome.Mmusculus.UCSC.mm9)
#seqlengths(reduced_fDbRepeats) <- seqlengths(Mmusculus)[names(seqlengths(reduced_fDbRepeats))]
#pull out ones on chr15
#reduced_fDbRepeats_chr15 <- reduced_fDbRepeats[seqnames(reduced_fDbRepeats) == "chr15"]
#pull out first 5 to play with.  All five ranges do not have overlapping repeats i.e. rmsk coverage score of 1
#first5 <- reduced_fDbRepeats_chr15[1:5]
        
target_ranges <- first5

## list of bam files
bamFls <- list.files(bamDir, "bam$", full=TRUE)

#make list of sample names for later reference
names(bamFls) <- sub("\\..*", "", basename(bamFls))
samples <- names(bamFls)

############## 
#read in .bam files as gapped alignments
##############

align <- function(fl)
  
{
  aln <- as(readGappedAlignments(fl), "GRanges")
}

alignments <- lapply(bamFls, align)

#########
# make dataframe
#########



##################
# count number of overlapping reads
############## 

ovrlpCount <- function(aln, range)
  
{
  for (i in 1:length(range)) {
    hits <- queryHits(findOverlaps(aln, range[i], ignore.strand=FALSE))
    numberReads <- length(hits)
  }
  return(numberReads)
}
  
counts <- sapply(alignments, ovrlpCount, target_ranges)  


#######peak finding
  #this will work on only one chromosome at a time, must change to chr#
  cov <- coverage(reads)
  cov <- as(cov[15], "GRanges")
  range <- range(score(cov))
  height <- range[2]
  peak <- ranges(cov[score(cov) == height])
  assign(name_gr, reads, envir = .GlobalEnv)
  assign(name_counts, counts, envir = .GlobalEnv)
  assign(name_height, height, envir = .GlobalEnv)
  assign(name_peak, peak, envir = .GlobalEnv)
  ######plot individual samples
  #plot_reads <- autoplot(reads, facets = strand~seqnames, aes(fill = strand))
  #plot_cov <- autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
  #assign(name_reads, plot_reads, envir = .GlobalEnv)
  #assign(name_cov, plot_cov, envir = .GlobalEnv)
  ###### output to global variable
  
}
ovrlps <- sapply(bamFls, ovrlpReads, target_range)

#######peak finding
#cov <- coverage(reads)
#cov <- as(cov[15], "GRanges")
#range <- range(score(cov))
#height <- range[2]
#peak <- ranges(cov[score(cov) == height])

###### Dataframe
stats_fileName <- paste(rangeName, "stats.rda", sep="_")
samples <- c("DPBS_1","DPBS_2", "DPBS_3", "DPBS_4", "DPBS_5", "DPBS_6", "DPBS_7", "DPBS_8", "Pru_1", "Pru_2", "Pru_3", "Pru_4", "Pru_5")
counts <- c(DPBS_1_chr15_count, DPBS_2_chr15_count, DPBS_3_chr15_count, DPBS_4_chr15_count, DPBS_5_chr15_count, DPBS_6_chr15_count, DPBS_7_chr15_count, DPBS_8_chr15_count, Pru_1_chr15_count, Pru_2_chr15_count, Pru_3_chr15_count, Pru_4_chr15_count, Pru_5_chr15_count)
heights <- c(DPBS_1_chr15_height, DPBS_2_chr15_height, DPBS_3_chr15_height, DPBS_4_chr15_height, DPBS_5_chr15_height, DPBS_6_chr15_height, DPBS_7_chr15_height, DPBS_8_chr15_height, Pru_1_chr15_height, Pru_2_chr15_height, Pru_3_chr15_height, Pru_4_chr15_height, Pru_5_chr15_height)
peaks <- c(DPBS_1_chr15_peak, DPBS_2_chr15_peak, DPBS_3_chr15_peak, DPBS_4_chr15_peak, DPBS_5_chr15_peak, DPBS_6_chr15_peak, DPBS_7_chr15_peak, DPBS_8_chr15_peak, Pru_1_chr15_peak, Pru_2_chr15_peak, Pru_3_chr15_peak, Pru_4_chr15_peak, Pru_5_chr15_peak)
stats <- data.frame(list(sample= samples, count = counts, height = heights, peakStart = start(peaks), peakEnd = end(peaks)))
save(stats, file=stats_fileName)

###COMBINED HITS
DPBS_hits <- c(DPBS_1_chr15_hits, DPBS_2_chr15_hits, DPBS_3_chr15_hits, DPBS_4_chr15_hits, DPBS_5_chr15_hits, DPBS_6_chr15_hits, DPBS_7_chr15_hits, DPBS_8_chr15_hits)
DPBS_hits
#GRanges with 38 ranges and 0 metadata columns:

Pru_hits <- c(Pru_1_chr15_hits, Pru_2_chr15_hits, Pru_3_chr15_hits, Pru_4_chr15_hits, Pru_5_chr15_hits)
Pru_hits
#GRanges with 26 ranges and 0 metadata columns:

######## REPEATMASKER HITS
#find the original repeat coordinates from repeatMasker
load("/Users/Miranda/Dropbox/RepeatAnalysis/rmsk_files/featureDb_mm9Repeats.rda")
rmskReps <- queryHits(findOverlaps(featureDb_mm9Repeats, target_range, ignore.strand=FALSE))
rmsk_hits = featureDb_mm9Repeats[rmskReps]

##### Repeats track
repeats <- autoplot(rmsk_hits, facets = strand~seqnames, aes(fill = strand))

####### Combined tracks

DPBS_reads <- autoplot(DPBS_hits, facets = strand~seqnames, aes(fill = strand))
Pru_reads <- autoplot(Pru_hits, facets = strand~seqnames, aes(fill = strand))

DPBS_cov <- autoplot(DPBS_hits, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
Pru_cov <- autoplot(Pru_hits, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))

####### output files

#Save to a pdf
DPBS_combined_fileName <- paste(rangeName, "DPBS_combined.pdf", sep="_")
pdf(file=DPBS_combined_fileName)
tracks(cntrl_reads = DPBS_reads, control_cov = DPBS_cov, repeats = repeats, heights = c(3,3,1))
dev.off() #devices off  - stop sending to the PDF

Pru_combined_fileName <- paste(rangeName, "Pru_combined.pdf", sep="_")
pdf(file= Pru_combined_fileName)
tracks(toxo_reads = Pru_reads, toxo_cov = Pru_cov, repeats = repeats, heights = c(3,3,1))
dev.off() 

###### more graphs

#pdf(file="control_reads.pdf")
#tracks(DPBS_1 = DPBS_1_chr15_reads, DPBS_2 = DPBS_2_chr15_reads, DPBS_3 = DPBS_3_chr15_reads, DPBS_4 = DPBS_4_chr15_reads, DPBS_5 = DPBS_5_chr15_reads, DPBS_6 = DPBS_6_chr15_reads, DPBS_7 = DPBS_7_chr15_reads, DPBS_8 = DPBS_8_chr15_reads) 
#dev.off()
#above won't work - too many tracks?

#pdf(file="toxo_reads.pdf")
#tracks(Pru_1 = Pru_1_chr15_reads, Pru_2 = Pru_2_chr15_reads, Pru_3 = Pru_3_chr15_reads, Pru_4 = Pru_4_chr15_reads, Pru_5 = Pru_5_chr15_reads)
#dev.off()
#above doesn't look pretty because of rectangle height issue (try to fix on next run w/variable in function)

#Following gives ERROR: Error in res$xlim[1] : object of type 'closure' is not subsettable
#pdf(file="control_cov1.pdf")
#tracks(DPBS_1 = DPBS_1_chr15_cov, DPBS_2 = DPBS_2_chr15_cov, DPBS_3 = DPBS_3_chr15_cov, DPBS_4 = DPBS_4_chr15_cov)  
#dev.off()
#pdf(file="control_cov2.pdf")
#tracks(DPBS_5 = DPBS_5_chr15_cov, DPBS_6 = DPBS_6_chr15_cov, DPBS_7 = DPBS_7_chr15_cov, DPBS_8 = DPBS_8_chr15_cov)

#pdf(file="toxo_cov.pdf")
#tracks(Pru_1 = Pru_1_chr15_cov, Pru_2 = Pru_2_chr15_cov, Pru_3 = Pru_3_chr15_cov, Pru_4 = Pru_4_chr15_cov, Pru_5 = Pru_5_chr15_cov)
#dev.off()
#######################
#clear variables
###################
#rm(list=ls())
#gc()
