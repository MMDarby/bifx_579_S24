library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ggbio)
data(genesymbol, package = "biovizBase")


#output directory
outDir <- "/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/"
setwd(outDir)
load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/cntrlGR.rda")
#Plot genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts <- transcripts(txdb)
trim9repeat <- GRanges(seqnames = "chr14", ranges = IRanges(start = 51460969, end =51463160 ), strand = "*")
trim9 <- queryHits(findOverlaps(transcripts, trim9repeat, ignore.strand=TRUE))
trim9 <- transcripts[trim9]
trim9 <- reduce(trim9)
#plot transcripts
#plot directly from txdb - "which" must be gene range, not repeat or other
#p.ideo <- plotIdeogram(genome = "hg19")
trim9Plot <- autoplot(txdb,which=trim9, names.expr="gene_id")
reducedTrim9 <- autoplot(txdb,which=trim9,stat="reduce", color="brown", fill="brown")
# pdf(file="trim9test.pdf")
# tracks(full=trim9Plot, reduce=reducedTrim9, heights=c(5,1)) + ylab("") + theme_tracks_sunset()
# dev.off()
# findReads <- function(GR,range,name){
#   ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=TRUE))
#   ovrlp <- GR[ovrlp]
#   assign(name,ovrlp)
# }
# findReads(GR=cntrlGR,range=trim9,name="trim9reads")

trim9reads <- queryHits(findOverlaps(cntrlGR, trim9, ignore.strand=TRUE))
trim9reads <- cntrlGR[trim9reads]
plot_reads <- autoplot(trim9reads, facets = strand~seqnames, aes(fill = strand))
plot_cov <- autoplot(trim9reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))

# trim9trans <- queryHits(findOverlaps(transcripts, trim9, ignore.strand=TRUE))
# trim9trans <- transcripts[trim9trans]
# plot_trans <- autoplot(trim9trans, facets = strand~seqnames, aes(fill = strand))
# 
# pdf(file="trim9reads.pdf")
# tracks(cnrlReads=plot_reads, transcripts=plot_trans, reduced=reducedTrim9, heights=c(5.5,1.5,1))
# dev.off()

#trim9Plot <- autoplot(txdb,which=trim9, names.expr="gene_id", facets = strand~seqnames, aes(fill = strand))

pdf(file="trim9reads2.pdf")
tracks(cnrlReads=plot_reads, transcripts=trim9Plot, reduced=reducedTrim9, heights=c(5.5,1.5,1))
dev.off()

pdf(file="trim9reads3.pdf")
tracks(controlReads=plot_cov, transcripts=trim9Plot, heights=c(7,1))
dev.off()

p1 <- ggplot() + geom_alignment(txdb, which = genesymbol["TRIM9"]) 
ggplot(txdb) + geom_alignment(which = genesymbol["TRIM9"], names.expr = "tx_id:::gene_id")
pdf(file="trim9reads4.pdf")
tracks(transcripts=p1)
dev.off()
#way too slow:
#reads <- ggplot() + geom_alignment(cntrlGR, which = genesymbol["TRIM9"],facets = strand ~ seqnames, aes(color = strand, fill = strand))
pdf(file="trim9reads5.pdf")
tracks(controlReads=reads, transcripts=p1, heights=c(7,1))
dev.off()

################
# new session
#################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ggbio)
data(genesymbol, package = "biovizBase")

#output directory
outDir <- "/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/"
setwd(outDir)
load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/cntrlGR.rda")
#Plot genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

wh <- genesymbol["DNM1L"]

findReads <- function(GR,range){
  ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=TRUE))
  ovrlp <- GR[ovrlp]
}
reads <- findReads(GR=cntrlGR,range=wh)

trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))

pdf(file="DNM1Lcontrol.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

# plotGene <- function(GR,gene){
#   wh <- genesymbol[gene]
#   reads <- findReads(GR=GR,range=wh)
#   trans <- ggplot() + geom_alignment(txdb, which = wh) 
#   ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
#   readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
#   fileName <- paste(gene,".pdf",sep="")
#   pdf(file=fileName)
#   tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
#   dev.off()
# }
# plotGene(GR=cntrlGR,gene="DNM1L")

GR <- cntrlGR
gene <- "DNM1L"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "GAPDH"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "HIF3A"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "BLOC1S3"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "KIAA1671"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "NRG1"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "CDH1"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "RGS5"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

#############################
#Plot tophat reads
################################
getwd()
#[1] "/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep"
F1topHat <- as(readGAlignments("/dcs01/stanley/work/OrbFTopHat/F1/accepted_hits.bam"), "GRanges")
save(F1topHat, file="F1topHatGR.rda")
F2topHat <- as(readGAlignments("/dcs01/stanley/work/OrbFTopHat/F2/accepted_hits.bam"), "GRanges")
save(F2topHat, file="F2topHatGR.rda")
F3topHat <- as(readGAlignments("/dcs01/stanley/work/OrbFTopHat/F3/accepted_hits.bam"), "GRanges")
save(F3topHat, file="F3topHatGR.rda")
F4topHat <- as(readGAlignments("/dcs01/stanley/work/OrbFTopHat/F4/accepted_hits.bam"), "GRanges")
save(F4topHat, file="F4topHatGR.rda")
F5topHat <- as(readGAlignments("/dcs01/stanley/work/OrbFTopHat/F5/accepted_hits.bam"), "GRanges")
save(F5topHat, file="F5topHatGR.rda")

ggbiosave <- function (filename = default_name(plot), plot = last_plot(),
                       device = default_device(filename), path = NULL,
                       scale = 1,
                       width = par("din")[1], height = par("din")[2],
                       units = c("in",
                                 
                                 "cm", "mm"), dpi = 300, limitsize = TRUE, ...)
{
  ## simply comment out the check part
  ## if (!inherits(plot, "ggplot"))
  ##     stop("plot should be a ggplot2 plot")
  eps <- ps <- function(..., width, height) grDevices::postscript(...,
                                                                  width = width, height = height, onefile = FALSE, horizontal =
                                                                    FALSE,
                                                                  paper = "special")
  tex <- function(..., width, height) grDevices::pictex(...,
                                                        width = width, height = height)
  pdf <- function(..., version = "1.4") grDevices::pdf(...,
                                                       version = version)
  svg <- function(...) grDevices::svg(...)
  wmf <- function(..., width, height) grDevices::win.metafile(...,
                                                              width = width, height = height)
  emf <- function(..., width, height) grDevices::win.metafile(...,
                                                              width = width, height = height)
  png <- function(..., width, height) grDevices::png(..., width = width,
                                                     height = height, res = dpi, units = "in")
  jpg <- jpeg <- function(..., width, height) grDevices::jpeg(...,
                                                              width = width, height = height, res = dpi, units = "in")
  bmp <- function(..., width, height) grDevices::bmp(..., width = width,
                                                     height = height, res = dpi, units = "in")
  tiff <- function(..., width, height) grDevices::tiff(...,
                                                       width = width, height = height, res = dpi, units = "in")
  default_name <- function(plot) {
    paste(digest.ggplot(plot), ".pdf", sep = "")
  }
  default_device <- function(filename) {
    pieces <- strsplit(filename, "\\.")[[1]]
    ext <- tolower(pieces[length(pieces)])
    match.fun(ext)
  }
  units <- match.arg(units)
  convert_to_inches <- function(x, units) {
    x <- switch(units, `in` = x, cm = x/2.54, mm = x/2.54/10)
  }
  convert_from_inches <- function(x, units) {
    x <- switch(units, `in` = x, cm = x * 2.54, mm = x *
                  2.54 * 10)
  }
  if (!missing(width)) {
    width <- convert_to_inches(width, units)
  }
  if (!missing(height)) {
    height <- convert_to_inches(height, units)
  }
  if (missing(width) || missing(height)) {
    message("Saving ", prettyNum(convert_from_inches(width *
                                                       scale, units), digits = 3), " x ",
            prettyNum(convert_from_inches(height *
                                            scale, units), digits = 3), " ", units, " image")
  }
  width <- width * scale
  height <- height * scale
  if (limitsize && (width >= 50 || height >= 50)) {
    stop("Dimensions exceed 50 inches (height and width are specified
         in inches/cm/mm, not pixels).",
         " If you are sure you want these dimensions, use
         'limitsize=FALSE'.")
  }
  if (!is.null(path)) {
    filename <- file.path(path, filename)
  }
  device(file = filename, width = width, height = height, ...)
  on.exit(capture.output(dev.off()))
  print(plot)
  invisible()
}
##SAVING DOESN"T WORK
# plotGene <- function(GR,gene){
#   wh <- genesymbol[gene]
#   reads <- findReads(GR=GR,range=wh)
#   trans <- ggplot() + geom_alignment(txdb, which = wh) 
#   ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
#   readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
#   #fileName <- paste(gene,"topHatF1.pdf",sep="")
#   p1=tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
#   ggbiosave(filename="RGS5topHatF1.pdf", plot=p1)
# }
# plotGene(GR=F1topHat,gene="RGS5")

GR <- F1topHat
gene <- "RGS5"
wh <- genesymbol[gene]
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,"topHatF1.pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()

GR <- cntrlGR
gene <- "TRIM9repeat"
wh <- GRanges(seqnames = "chr14", ranges = IRanges(start = 51461970, end =51462162), strand = "*")
reads <- findReads(GR=GR,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
fileName <- paste(gene,".pdf",sep="")
pdf(file=fileName)
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()


##########################################
findAndPlot <- function(GR,rangeList){
  counter=1
  for(r in seq_along(rangeList)){
    # find overlaps between reads and target range
    ovrlp <- queryHits(findOverlaps(GR, rangeList[r], ignore.strand=TRUE))
    # make GRanges object of overlapping reads
    reads = GR[ovrlp]
    #reads are divided and color coded based on strand
    plot_reads <- autoplot(reads, facets = strand~seqnames, aes(fill = strand))
    #creates a ggbio track that plots coverage of reads on each strand
    plot_cov <- autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand)) 
    chr <- as.character(seqnames(rangeList[r]))
    start <- start(rangeList[r])
    end <- end(rangeList[r])
    name <- paste(chr,start,end,sep="_")
    fileName <-paste(name,".pdf",sep="")
    pdf(file=fileName)
    # use GGBIO to plot tracks created above to show reads and coverage of combined control samples, plus regions from repeatMasker
    # call is: tracks(title = track, heights = list of relative track heights)
    tracks(reads = plot_reads, coverage = plot_cov)
    dev.off() #devices off  - stop sending to the PDF
  }
  counter=counter+1
}
findAndPlot(GR=cntrlGR,rangeList=short[1])
