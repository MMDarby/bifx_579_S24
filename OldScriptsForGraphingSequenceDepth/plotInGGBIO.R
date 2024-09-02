library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ggbio)
data(genesymbol, package = "biovizBase")

# Set output directory
outDir <- "/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/"
setwd(outDir)
#load files with GRanges
#all control reads combined (slow):
load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/cntrlGR.rda")
#individual MQ40 files (fast):
load("/scratch/temp/Stanley/miranda/repeatAnalysis/human/SarvenRMSK/orbFrontalF1alignGR.rda")
######################################
#HOW GR WERE GENERATED
#########
# #all .bam files in bamDir
# bamFls <- list.files(bamDir, "bam$", full=TRUE)
# #save GR of all reads in each bam file in the list "fl"
# readBams <- function(fl){
#   counter=1
#   for(f in seq_along(fl)){
#     name <- sub("\\..*", "", basename(fl[f]))
#     print(name)
#     fileName <- paste(name,"GR.rda",sep="")
#     # import .bam file to GRanges object
#     aln <- as(readGAlignments(fl[f]), "GRanges")
#     assign(name, aln, envir=.GlobalEnv)
#     save(list=name,file=fileName)
#   }
#   counter=counter+1
# }
# readBams(bamFls)
#######################################

#Plot genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#Find all reads in GRanges that overlap region
findReads <- function(GR,range){
  ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=TRUE))
  ovrlp <- GR[ovrlp]
}

# GRanges object to plot
GR <- cntrlGR
# gene name as string
gene <- "RGS5"
#Set region to plot - can also be in form of: 
#wh <- GRanges(seqnames = "chr14", ranges = IRanges(start = 51461970, end =51462162), strand = "*")
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
