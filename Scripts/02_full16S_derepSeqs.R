#################################################
### Get reverse complements, find identical
### sequences, dereplicate into unique ASVs.
### PacBio Sequel II full 16S rRNA gene amplicons
### Lou LaMartina, finalized 22 March 2022
#################################################


setwd("~/Desktop/Lab/Projects/PacBio/FINAL")
library(ggplot2)


#################
### load data ###
#################

# ASV table from 01_full16S_removePrimers.R
counts_uniq <- read.csv("./RData/01_counts.csv")
rownames(counts_uniq) <- counts_uniq$orginalFile
counts_uniq <- counts_uniq[-1]


# taxa table from 01_full16S_removePrimers.R
taxa_uniq <- read.csv("./RData/01_taxa.csv")




###############################
### get reverse complements ###
###############################

# making sure that ASVs are not reverse complements of each other

### function ### - takes ~1min
# https://www.r-bloggers.com/2008/11/r-function-to-reverse-and-complement-a-dna-sequence/

reverseComplement.fct <- function(x, rev = TRUE)
{
  x <- toupper(x)
  y <- rep("N", nchar(x))
  xx <- unlist(strsplit(x, NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb] == "A") y[bbb] <- "T"
    if(xx[bbb] == "C") y[bbb] <- "G"
    if(xx[bbb] == "G") y[bbb] <- "C"
    if(xx[bbb] == "T") y[bbb] <- "A"
  }
  if(rev == FALSE)
  {
    for (ccc in (1:nchar(x)))
    {
      if(ccc == 1) yy <- y[ccc] else yy <- paste(yy, y[ccc], sep = "")
    }
  }
  if(rev == TRUE)
  {
    zz <- rep(NA, nchar(x))
    for (ccc in (1:nchar(x)))
    {
      zz[ccc] <- y[nchar(x) + 1 - ccc]
      if(ccc == 1) yy <- zz[ccc] else yy <- paste(yy, zz[ccc], sep = "")
    }
  }
  return(yy)
}


# apply function
start = Sys.time()
revComps.ls <- list()
for (i in 1:length(taxa_uniq$UNIQ)) {
  revComps.ls[[i]] <- c(taxa_uniq$UNIQ[i], reverseComplement.fct(taxa_uniq$R1[i]))
}
Sys.time() - start
# 48.09699 secs


# extract
revComps <- data.frame(do.call(rbind, revComps.ls))
colnames(revComps) <- c("UNIQ", "R2")
taxa_uniq <- merge(taxa_uniq, revComps, by = "UNIQ")
identical(colnames(counts_uniq), taxa_uniq$UNIQ)


# did that work? (randomly select 10, get rc, check if same)
for (i in taxa_uniq$UNIQ[sample(1:length(taxa_uniq$UNIQ), 10)]) {
  cat("\n", i, " ", taxa_uniq$R2[taxa_uniq$UNIQ == i] == 
        reverseComplement.fct(taxa_uniq$R1[taxa_uniq$UNIQ == i]))
}




#################################
### dereplicate reverse comps ###
#################################

# get duplicates
dup_seqs = unique(taxa_uniq$R1)
dup_seqs.ls <- list()
for (i in dup_seqs) {
  dup_seqs.ls[[i]] <- sort(c(taxa_uniq$UNIQ[taxa_uniq$R1 == i], taxa_uniq$UNIQ[taxa_uniq$R2 == i]))
}


# convert to data frame
compares <- data.frame(do.call(rbind, dup_seqs.ls))
colnames(compares) <- c("UNIQ", "UNIQ2")
compares <- unique(merge(compares, taxa_uniq[c(1,9,10)], by = "UNIQ"))
colnames(compares)[1] <- "UNIQ1"


# variable of identical asvs
compares$compare <- paste0(compares$UNIQ1, "__", compares$UNIQ2)


# subset ASVs that do not need to be summed
noDup_seqs <- compares[compares$UNIQ1 == compares$UNIQ2,]
counts_noDup <- counts_uniq[colnames(counts_uniq) %in% noDup_seqs$UNIQ1]


# subset those that do
dup_seqs <- compares[compares$UNIQ1 != compares$UNIQ2,]


# sum counts of duplicate/rc sequences
dups_counts.ls <- list()
for (i in dup_seqs$R1) {
  dups_counts.ls[[i]] <- counts_uniq[dup_seqs$UNIQ1[dup_seqs$R1 == i]] +
    counts_uniq[dup_seqs$UNIQ2[dup_seqs$R1 == i]]
}


# convert to data frame, add singlet asvs
counts_noDup <- cbind(data.frame(do.call(cbind, dups_counts.ls)), counts_noDup)


# make column names R1 reads
taxa_noDup <- subset(taxa_uniq, UNIQ %in% colnames(counts_noDup))
taxa_noDup <- taxa_noDup[match(colnames(counts_noDup), taxa_noDup$UNIQ),]


# add comparisons
colnames(taxa_noDup)[1] <- "UNIQ1"
taxa_noDup <- merge(taxa_noDup, compares[c(1,5)], by = "UNIQ1")




############
### save ###
############

# match
taxa_noDup <- taxa_noDup[match(colnames(counts_noDup), taxa_noDup$UNIQ1),]
identical(colnames(counts_noDup), taxa_noDup$UNIQ1)


# reorder
taxa_noDup <- taxa_noDup[c("UNIQ1", "compare", "Kingdom", "Phylum", "Class",
                           "Order", "Family", "Genus", "Species", "R1", "R2")]


# add file names
counts_noDup <- data.frame(orginalFile = rownames(counts_noDup), counts_noDup)


# write
write.csv(counts_noDup, "./RData/02_counts.csv", row.names = F)
write.csv(taxa_noDup, "./RData/02_taxa.csv", row.names = F)

