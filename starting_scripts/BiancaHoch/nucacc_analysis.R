###########################################
#
# Read and Analyze Snaptron output files from the nucleus accumbens queries made by makeScriptsSnaptron.R.
#
###########################################
rm(list=ls()) #Clear R memory
gc() #empty garbage collector

###Set working directory:
setwd("/Users/biancahoch/Documents/snaptronwork/snapout/nucacc_gtex/")

#Update file path to point toward appropriate folder on the computer
nucacc_list <- list.files(pattern=".tsv") 

#Read files:
nucacc <- sapply(nucacc_list, read.delim, simplify = FALSE)

#Print only the rows in which the number of samples that have one or more reads covering the intron constitute atleast 5% of the total samples searched:

#First find the number of anterior cingulate cortex samples:
load("~/Documents/splice_project/data/analysis/nucacc_rail.Rdata")
nucacc_rail
splitnucacc <- strsplit(nucacc_rail, ",")
unlistsplitnucacc <- unlist(splitnucacc)
print(paste(length(unlistsplitnucacc)))

#See which junctions are present in atleast 5% of samples:
f <- function(df){
  (df[which(df$samples_count >= length(unlistsplitnucacc) *.05), ])
}

#Apply function to list of dataframes:
dfList <- lapply( X=nucacc, FUN=f ) 
#Bind the lists into one dataframe:
allDF <-do.call(rbind,dfList) 
#Write the dataframe to a .csv:
write.csv(allDF, "~/Documents/snaptronwork/snapout_analysis/nucacc_fivepercent.csv")
#Check output:
#nucacc_fivepercent <- read.csv("~/Documents/snaptronwork/snapout_analysis/nucacc_fivepercent.csv", header = TRUE, sep = ",")

#See which junctions are present in atleast 60% of samples:
f2 <- function(df){
  (df[which(df$samples_count >= length(unlistsplitnucacc) *.60), ])
}
#Apply function to list of dataframes:
dfList2 <- lapply( X=nucacc, FUN=f2 ) 
#Bind the lists into one dataframe:
allDF2 <-do.call(rbind,dfList2) 
#Write the dataframe to a .csv: (NOTE: ADD ROWNAMES)
write.csv(allDF2, "~/Documents/snaptronwork/snapout_analysis/nucacc_sixtypercent.csv")
#Check output:
#nucacc_sixtypercent <- read.csv("~/Documents/snaptronwork/snapout_analysis/nucacc_sixtypercent.csv", header = TRUE, sep = ",")

#############################################################
