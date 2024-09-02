##################################
# creates GRanges list of all reads in .bam that overlap with genomic region
# populates data frame with statistics about overlapping reads
# plots overlapping read alignments over the region, for both individual and combined samples
# plots repeats in the genomic range that are designated by repeatMasker
# created by Miranda Darby
# modified 3-7-13
##################################

library(GenomicRanges)
library(ggbio)

######## things to change ############
# 4 things to change: bamDir, setwd(), target_range, rangeName

#Set directories 
# directory with .bam files
#bamDir <- "/net/fs01/mnt/storage/home/mdarby/indexedMouse1BamFiles"
bamDir <- "~/Dropbox/RepeatAnalysis/bamFiles"
#working directory
#setwd("/net/fs01/mnt/storage/home/mdarby/repeatAnalysis/function")
setwd("~/Dropbox/RepeatAnalysis/plotCounts/function")

#Set range of interest MUST BE STRAND SPECIFIC and match repeat range FOR COUNTS
#target_range <- GRanges(seqnames = "chr15", ranges = IRanges(start = 57859370, end = 57859500), strand = "+")

#### Alternative way to set target to equal to range from reducedRepeats
#Import reduced set of repeat coordinates
#redRep <- load("/Users/Miranda/Dropbox/RepeatAnalysis/rmsk_files/reduced_fDbRepeats.rda")
#find section that overlaps with the "four overlapping repeats" 
#overlapArea <- queryHits(findOverlaps(reduced_fDbRepeats, GRanges(seqnames = "chr15", ranges = IRanges( start = 57859370, end = 57859500))))
#ranges = reduced_fDbRepeats[overlapArea]
#####SET + or - for target_range and rangeName
#target_range <- ranges[strand(ranges) == "+", ]

#### Another alternative - choose target range from a list of ranges
#load("chr15plusThree.rda")
load("chr15minusThree.rda")
target_range <- chr15_minusThree_regions[14]
target_range

#### Check that rangeName variable accurately reflects target range and strand - used for naming files
rangeName <- paste((seqnames(target_range)), (start(target_range)), (end(target_range)), "m", sep="_")



#save files ending in bam to bamFls variable
bamFls <- list.files(bamDir, "bam$", full=TRUE)

# the following is a monster function that begs for optimization
ovrlpReads <- function(fl, range)
{
  #create variable names
  names(fl) <- sub("\\..*", "", basename(fl))
  name <- names(fl)
  name_gr <- paste(name,"hits", sep="_")
  name_reads <- paste(name, "reads", sep="_")
  name_cov <- paste(name, "cov", sep="_")
  name_counts <- paste(name, "count", sep="_")
  name_height <- paste(name, "height", sep="_")
  name_peak <- paste(name, "peak", sep="_")
  print(name)
  # import .bam file to GRanges object
  aln <- as(readGappedAlignments(fl), "GRanges")
  # find overlaps between reads and target range
  hits <- queryHits(findOverlaps(aln, range, ignore.strand=FALSE))
  # make GRanges object of overlapping reads
  reads = aln[hits]
  # count number of reads
  counts <- length(reads)
  
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
  
  ######plot individual sample tracks and output to global variable
  #creates a ggbio track that displays all of the reads as stacked rectangles piled on genomic coordinates
  #reads are divided and color coded based on strand
  plot_reads <- autoplot(reads, facets = strand~seqnames, aes(fill = strand))
  #creates a ggbio track that plots coverage of reads on each strand
  plot_cov <- autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
  # assign tracks to global variables that can be called later by ggbio 
  assign(name_reads, plot_reads, envir = .GlobalEnv)
  assign(name_cov, plot_cov, envir = .GlobalEnv)
  
}
ovrlps <- sapply(bamFls, ovrlpReads, target_range)


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

Pru_hits <- c(Pru_1_chr15_hits, Pru_2_chr15_hits, Pru_3_chr15_hits, Pru_4_chr15_hits, Pru_5_chr15_hits)
Pru_hits

######## REPEATMASKER HITS
#find the original repeat coordinates from repeatMasker
load("/Users/Miranda/Dropbox/RepeatAnalysis/rmsk_files/featureDb_mm9Repeats.rda")
rmskReps <- queryHits(findOverlaps(featureDb_mm9Repeats, target_range, ignore.strand=FALSE))
rmsk_hits = featureDb_mm9Repeats[rmskReps]

##### Repeats track
# plot repeats as rectangles stacked on coordinates, color and separate by strand
repeats <- autoplot(rmsk_hits, facets = strand~seqnames, aes(fill = strand))

####### Combined tracks
# plot combined as rectangles stacked on coordinates, color and separate by strand
DPBS_reads <- autoplot(DPBS_hits, facets = strand~seqnames, aes(fill = strand))
Pru_reads <- autoplot(Pru_hits, facets = strand~seqnames, aes(fill = strand))
# plot coverage of combined reads, color and separate by strand
DPBS_cov <- autoplot(DPBS_hits, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
Pru_cov <- autoplot(Pru_hits, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))

####### output files

#Create graphic from tracks created above and save to a pdf
# create file name for pdf
DPBS_combined_fileName <- paste(rangeName, "DPBS_combined.pdf", sep="_")
pdf(file=DPBS_combined_fileName)
# use GGBIO to plot tracks created above to show reads and coverage of combined control samples, plus regions from repeatMasker
# call is: tracks(title = track, heights = list of relative track heights)
tracks(cntrl_reads = DPBS_reads, control_cov = DPBS_cov, repeats = repeats, heights = c(3,3,1))
dev.off() #devices off  - stop sending to the PDF

Pru_combined_fileName <- paste(rangeName, "Pru_combined.pdf", sep="_")
pdf(file= Pru_combined_fileName)
tracks(toxo_reads = Pru_reads, toxo_cov = Pru_cov, repeats = repeats, heights = c(3,3,1))
dev.off() 

###### more graphs #Creates PDFs with stacked plots showing data from individual samples over the gentic region  
# Will only work if all of the samples have reads that overlap the region
#plots with single reads don't look pretty because of rectangle height issue (try to fix with rectangle height variable when writing original plot?)

DPBS_singles_fileName1 <- paste(rangeName, "DPBS_singles1.pdf", sep="_")
pdf(file=DPBS_singles_fileName1)
tracks(DPBS_1 = DPBS_1_chr15_reads, DPBS_2 = DPBS_2_chr15_reads, DPBS_3 = DPBS_3_chr15_reads, DPBS_4 = DPBS_4_chr15_reads) 
dev.off()

DPBS_singles_fileName2 <- paste(rangeName, "DPBS_singles2.pdf", sep="_")
pdf(file=DPBS_singles_fileName2)
tracks(DPBS_5 = DPBS_5_chr15_reads, DPBS_6 = DPBS_6_chr15_reads, DPBS_7 = DPBS_7_chr15_reads, DPBS_8 = DPBS_8_chr15_reads) 
dev.off()

Pru_singles_fileName <- paste(rangeName, "Pru_singles.pdf", sep="_")
pdf(file=Pru_singles_fileName)
tracks(Pru_1 = Pru_1_chr15_reads, Pru_2 = Pru_2_chr15_reads, Pru_3 = Pru_3_chr15_reads, Pru_4 = Pru_4_chr15_reads, Pru_5 = Pru_5_chr15_reads)
dev.off()

DPBS_sing_cov_fileName1 <- paste(rangeName, "DPBS_sing_cov_1.pdf", sep="_")
pdf(file=DPBS_sing_cov_fileName1)
tracks(DPBS_1 = DPBS_1_chr15_cov, DPBS_2 = DPBS_2_chr15_cov, DPBS_3 = DPBS_3_chr15_cov, DPBS_4 = DPBS_4_chr15_cov)  
dev.off()

DPBS_sing_cov_fileName2 <- paste(rangeName, "DPBS_sing_cov_2.pdf", sep="_")
pdf(file=DPBS_sing_cov_fileName2)
tracks(DPBS_5 = DPBS_5_chr15_cov, DPBS_6 = DPBS_6_chr15_cov, DPBS_7 = DPBS_7_chr15_cov, DPBS_8 = DPBS_8_chr15_cov)

Pru_sing_cov_fileName <- paste(rangeName, "Pru_sing_cov.pdf", sep="_")
pdf(file=Pru_sing_cov_fileName)
tracks(Pru_1 = Pru_1_chr15_cov, Pru_2 = Pru_2_chr15_cov, Pru_3 = Pru_3_chr15_cov, Pru_4 = Pru_4_chr15_cov, Pru_5 = Pru_5_chr15_cov)
dev.off()

#######################
#clear variables
###################
rm(list=ls())
gc()
