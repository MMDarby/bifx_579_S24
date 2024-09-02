getwd()
[1] "/net/fs01/mnt/storage/home/mdarby/repeatAnalysis"

library(GenomicRanges)
library(ggbio)

load("DPBS_1_repeats.rda")
load("Pru_1_repeats.rda")

DPBS_1_repeats_genome <- autoplot(DPBS_1_repeats, coord = "genome")
Pru_1_repeats_genome <- autoplot(Pru_1_repeats, coord = "genome")

target_range <- GRanges(seqnames = "chr15", ranges = IRanges( start = 57859370, end = 57859500), strand = "+")

pdf(file="testPlot.pdf")
tracks(control = DPBS_1_repeats_genome, toxo = Pru_1_repeats_genome) + xlim(target_range)
# took 8 minutes
dev.off()
#did not get correct plot

#try transformToGenome
DPBS_1_repeats_tG <- transformToGenome(DPBS_1_repeats)
#Error: could not find function "transformToGenome"
# tried updating package but still got same error
# quit out and tried something else