library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ggbio)
data(genesymbol, package = "biovizBase")

#output directory
outDir <- "/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/"
setwd(outDir)
load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/flippedStrandControlRanges.rda")
#Plot genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
findReads <- function(GR,range){
  ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=TRUE))
  ovrlp <- GR[ovrlp]
}

wh <- genesymbol["TACC2"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="TACC2control.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["DNAJC15"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="DNAJC15control.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["ZNF746"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="ZNF746control.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["KHK"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="KHKcontrol.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

# REL peak is not apparent but there does seem to be an antisense transcript that matches the splicing patterns of KHK

wh <- genesymbol["HKR1"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="HKR1control.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["GAPDH"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="GAPDHcontrol.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["NR4A2"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="NR4A2control.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["CGREF1"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="CGREF1control_wholeGene.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["C12orf70"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="C12orf70control_wholeGene.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()

wh <- genesymbol["RCCD1"]
reads <- findReads(GR=flpCntrlGr,range=wh)
trans <- ggplot() + geom_alignment(txdb, which = wh) 
ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
readsPlot <-  autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
pdf(file="RCCD1control_wholeGene.pdf")
tracks(reads=readsPlot, transcripts=trans, heights=c(7,1))
dev.off()
dev.off()
rm(wh)
gc()


load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/flpBPrin7remGR.rda")
load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/CoveragePlots/rin7rem/rin7remMQ40normalization.rda")
library(ggplot2)
library(grid)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plotCovCvBP <- function(chrN,start,end,strand){
  range <- c(start:end)
  chr <- paste("chr",chrN,sep="")
  chr <- c(rep(chr,length(start)[]))
  str <- c(rep(strand,length(start)[]))
  gr <- GRanges(seqnames=chr,ranges=IRanges(start=range,end=range),strand=str)
  cntlCount <- countOverlaps(gr,flpCntrlGr,type="any",ignore.strand=FALSE)
  bpCount <- countOverlaps(gr,rin7remBP,type="any",ignore.strand=FALSE)
  control <- cntlCount/ctlN
  bipolar <- bpCount/bpN
  covCounts <- as.data.frame(cbind(control,bipolar,range))
  limit <- max(max(control),max(bipolar))
  p1 <- ggplot(data=covCounts, aes(x=range, y=control, fill=range)) + 
    geom_bar(colour="black", fill="black", width=1, stat="identity") + 
    guides(fill=FALSE) + 
    theme(axis.title.y = element_text(size=20,lineheight=1,colour="black"),axis.text.y = element_text(size=16,colour="black"),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(1,1.5,0.1,0.5),"lines")) +
    ylab("Control\nreads per billion") + ylim(0,limit) + xlim(start,end) #+ 
  p2 <- ggplot(data=covCounts, aes(x=range, y=bipolar, fill=range)) + 
    geom_bar(colour="black", fill="black", width=1, stat="identity") + 
    guides(fill=FALSE) +
    theme(axis.title.y = element_text(size=20,lineheight=1,colour="black"),axis.text.y = element_text(size=16,colour="black"),axis.text.x = element_text(size=12,colour="black"),axis.title.x = element_text(size=20,lineheight=1,colour="black"),plot.margin = unit(c(0.1,1.5,0.5,0.5),"lines")) +
    xlab(paste("Coordinates on chromosome",chrN)) + ylab("Bipolar\nreads per billion") + ylim(0,limit) + xlim(start,end) #+ 
  #ggtitle("Combined Reads from Bipolar Samples")
  
  filename <- paste(chr,start,end,"CvBP_SpREL.pdf",sep="_")
  pdf(file=filename)
  multiplot(p1,p2,cols=1)
  dev.off() 
}
#chr2:27321106-27321359_-
plotCovCvBP(chrN=2, start=27321106, end=27321359, strand="-")
#chr12:27579086-27579377_+
plotCovCvBP(chrN=12, start=27579086, end=27579377, strand="+")
#chr5:68410807-68410990_-
plotCovCvBP(chrN=5, start=68410807, end=68410990, strand="-")
#chr10:123906276-123906664_-
plotCovCvBP(chrN=10, start=123906276, end=123906664, strand="-")
#chr15:91507008-91507266_+
plotCovCvBP(chrN=15, start=91507008, end=91507266, strand="+")
#chr13:43631034-43633286_+
plotCovCvBP(chrN=13, start=43631034, end=43633286, strand="+")
#chr7:149176844-149177007_-
plotCovCvBP(chrN=7, start=149176844, end=149177007, strand="-")




