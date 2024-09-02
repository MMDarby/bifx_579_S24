######################################
# View reads that overlap with repeats
######################################

setwd("~/Dropbox/RepeatAnalysis/")

library(rtracklayer)

############################
# loads practice data, makes and plays with practice track
###########################
data(targets)
targetRanges = IRanges(targets$start, targets$end)
targetTrack = with (targets, GRangesForUCSCGenome("hg18", chrom, targetRanges, strand, name, target))
genome(targetTrack)
head(seqlengths(targetTrack))
head(seqnames(targetTrack))
export(targetTrack, "targets.bed")
restoredTrack = import("targets.bed", asRangedData = FALSE)
subTargetTrack <- targetTrack[1]
#below is strait from vignette but doesn't work
view <- browserView(session, subTargetTrack * -10, pack = "targets")
#following works
browseGenome(targetTrack, range = subTargetTrack * -10)
################################
# make and view track from our data
################################
# Load GRanges object
load("Ddam_Liver_repeats.rda")
# Save GRanges as .bed file
export(Ddam_Liver_repeats, "Ddam_Liver_repeats.bed")

# choose genome browser
session <- browserSession("UCSC")
# to see the names of other supported browsers
# genomeBrowsers()
track(session, "Ddam_Liver_repeats") <- Ddam_Liver_repeats
# This should work but I gett a big error because of the extra chromosomes. Also it calls up HG19. Why? 

# Here is an attept to force this to work:
Ddam_Liver_repTrack = with (Ddam_Liver_repeats, GRangesForUCSCGenome("mm9", chrom=chrom, ranges=Ddam_Liver_repeats$IRanges, strand=Ddam_Liver_repeats$strand))
#Error in as.vector(x, mode) : 
#cannot coerce type 'closure' to vector of type 'any'
# Seems to give all mm9 chromosomes from beginning to end
genome(Ddam_Liver_repTrack) 
export(Ddam_Liver_repTrack, "Ddam_Liver_repTrack.bed")

####Following code doesn't work right
track(session, "Is_mm9") = Ddam_Liver_repTrack
#just grab first feature
sub_repTrack = Ddam_Liver_repTrack[1]
view = browserView(session, sub_repTrack * -10, pack ="Is_mm9")
#ABOVE CLEARLY DOESN"T WORK RIGHT
#Error in resolveTrackIndex(x, i) : Unknown track(s): Is_mm9
#In addition: Warning messages:
 # 1: In `start<-`(`*tmp*`, value = -887379443L) :
  #trimmed start values to be positive
#2: In `end<-`(`*tmp*`, value = 1084574876L) :
 # trimmed end values to be <= seqlengths

#Also doesn't work because it thinks it should be HG19
browseGenome(Ddam_Liver_repeats)

#worth a try
#browseGenome(Ddam_Liver_repeats, range = "mm9" * -10)
#gives error

#DdamLiverTrack <- trackSet(ranges(Ddam_Liver_repeats), chrom(Ddam_Liver_repeats), genome = "mm9")
#ERROR, could not find function "trackSet"
browseGenome(ranges(Ddam_Liver_repeats), chrom(Ddam_Liver_repeats), genome = "mm9")


########## Downloading tracks
loaded_tracks = trackNames(session)

########### viewving tracks
visible_tracks <- trackNames(view)
trackNames(view) = visible_tracks

###########trying again 2-19-13


> getwd()
[1] "/Users/Miranda/Dropbox/RepeatAnalysis"
> DdamLiver15 = load("Ddam_Liver_chr15.rda")
#NOTE: somehow the chr15 file was saved incorrectly so that only the string "chr15" was saved.

#Working on server:
getwd()
#[1] "/net/fs01/mnt/storage/home/mdarby/repeatAnalysis"
load("Ddam_Liver_chr15.rda")
chr15cov = coverage(chr15)
range(chr15cov)
#chr15
#[1,]     0
#[2,] 62613
> Ddam_Liver15_cov = as(chr15cov, "GRanges")
> Ddam_Liver15_cov
#GRanges with 252062 ranges and 1 metadata column:
range(Ddam_Liver15_cov$score)
#[1]     0 62613
readCount_Ddam_Liver_chr15 = subset(Ddam_Liver15_cov, Ddam_Liver15_cov$score>=1)
> readCount_Ddam_Liver_chr15
#GRanges with 234796 ranges and 1 metadata column:
#NOTE:THIS METHOD ELIMINATES STRAND INFO, would need to split by strand first to maintain strand info
> export(readCount_Ddam_Liver_chr15, "Ddam_Liver_chr15_readCount.bedGraph")
> export(readCount_Ddam_Liver_chr15, "Ddam_Liver_chr15_readCount.bed")

library(rtracklayer)
genome(session)="mm9"
trackNames(session)
#download repeatMasker track
rmskTrack = track(session, "rmsk")
head(rmskTrack)
genome(rmskTrack)
#shows that it is mm9
export(rmskTrack, "mm9rmsk.bed")
restoredRmsk = import("mm9rmsk.bed")
# ONLY SOME OF THE READS FROM CHR12 will make from mm9Repeats.rda on server
