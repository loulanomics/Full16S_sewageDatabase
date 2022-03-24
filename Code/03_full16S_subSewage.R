#################################################
### Get reverse complements, find identical
### sequences, dereplicate into unique ASVs.
### PacBio Sequel II full 16S rRNA gene amplicons
### Lou LaMartina, finalized 22 March 2022
#################################################


setwd("~/Desktop/Lab/Projects/PacBio/FINAL")
library(ggplot2)
library(dada2)


#################
### load data ###
#################

# ASV table from 02_full16S_derepSeqs.R
counts_noDup <- read.csv("./RData/02_counts.csv")
rownames(counts_noDup) <- counts_noDup$orginalFile
counts_noDup <- counts_noDup[-1]


# taxa table from 02_full16S_derepSeqs.R
taxa_noDup <- read.csv("./RData/02_taxa.csv")


# sample info
info <- read.csv("./RData/Wastewater_full16S_sample_metadata.csv")




#########################
### wastewater counts ###
#########################

# another dataset was included in the pacbio run,
# subset by sample/file names


# subset wastewater samples/ASVs
rownames(counts_noDup) <- gsub("__", "--", rownames(counts_noDup))
counts_ww <- counts_noDup[rownames(counts_noDup) %in% info$originalName,]
counts_ww <- counts_ww[colSums(counts_ww) > 0]


# sort by abundance, make asv variable called "WW"
counts_ww <- counts_ww[names(sort(colSums(counts_ww), decreasing = T))]
fastas_ww <- data.frame(WW = paste0("WW", sprintf("%04d", 1:ncol(counts_ww))),
                        UNIQ1 = colnames(counts_ww))


# save fasta and names file for mothur -
# see 04_full16S_clusterASV.R
temptax <- subset(taxa_noDup, taxa_noDup$UNIQ1 %in% colnames(counts_ww))
temptax <- temptax[match(colnames(counts_ww), temptax$UNIQ1),]
tempcount <- counts_ww
identical(colnames(tempcount), temptax$UNIQ1)
colnames(tempcount) <- temptax$R1

uniquesToFasta(as.matrix(tempcount), ids = temptax$UNIQ1,
               fout = "./mothur/Wastewater_full16S_23Mar22.fasta")
write.table(cbind(temptax$UNIQ1, temptax$UNIQ1), "./mothur/names", sep = "\t",
            row.names = F, col.names = F, quote = F)

rm(tempcount, temptax)


# add duplicate asvs
fastas_ww <- merge(fastas_ww, taxa_noDup[c(1,2,10,11)], by = "UNIQ1")


# change column names to ASV variable
fastas_ww <- fastas_ww[match(colnames(counts_ww), fastas_ww$UNIQ1),]
identical(colnames(counts_ww), fastas_ww$UNIQ1)
colnames(counts_ww) <- fastas_ww$WW




#######################
### wastewater taxa ###
#######################

# get UNIQ2
fastas_ww$UNIQ2 <- sapply(strsplit(fastas_ww$compare, "__"), '[', 2)


# load original data to decide best assignments
taxa_uniq <- read.csv("./RData/01_taxa.csv")
taxa_uniq[is.na(taxa_uniq)] <- "unclassified"


# taxa of first comparison
uniq1_tax <- taxa_uniq
uniq1_tax <- subset(uniq1_tax, UNIQ %in% fastas_ww$UNIQ1)
colnames(uniq1_tax)[1] <- "UNIQ1"
uniq1_tax <- merge(uniq1_tax, fastas_ww[-c(4:5)], by = "UNIQ1")


# taxa of first comparison
uniq2_tax <- taxa_uniq
uniq2_tax <- subset(uniq2_tax, UNIQ %in% fastas_ww$UNIQ2)
colnames(uniq2_tax)[1] <- "UNIQ2"
uniq2_tax <- merge(uniq2_tax, fastas_ww[-c(4:5)], by = "UNIQ2")


# with ones that have discrepancies
uniq_tax.ls <- list()
notUniq_tax.ls <- list()
for (i in fastas_ww$compare) {
  
  # if there are no discrepancies in tax assignments
  if (identical(paste(uniq1_tax[uniq1_tax$compare == i, 2:8], collapse = " "), 
                paste(uniq2_tax[uniq2_tax$compare == i, 2:8], collapse = " ")) == T) {
    uniq_tax.ls[[i]] <- uniq1_tax[uniq1_tax$compare == i,] }
  
  # if there is more than 1 tax assignment/they are different
  if (identical(paste(uniq1_tax[uniq1_tax$compare == i, 2:8], collapse = " "), 
                paste(uniq2_tax[uniq2_tax$compare == i, 2:8], collapse = " ")) == F) {
    notUniq_tax.ls[[i]] <- rbind(uniq1_tax[uniq1_tax$compare == i,],
                                 uniq2_tax[uniq2_tax$compare == i,]) }
  
}


# convert to data frame
uniq_tax <- data.frame(do.call(rbind, uniq_tax.ls))


# with ones that have discrepancies
for (i in names(notUniq_tax.ls)) {
  
  # add new column
  n = ncol(notUniq_tax.ls[[i]]) + 1
  
  # for each taxon assigned to a given asv
  for (j in 1:nrow(data.frame(notUniq_tax.ls[[i]]))) {
    
    # count how many instances of "unclassified" there are
    notUniq_tax.ls[[i]][j, n] <- length(grep("unclassified", data.frame(notUniq_tax.ls[[i]])[j,]))
  }
}
lastcol = colnames(notUniq_tax.ls[[i]])[ncol(notUniq_tax.ls[[i]])]


# choose asvs with least "unclassifieds"
uniq_tax.ls <- list()
for (i in names(notUniq_tax.ls)) {
  
  # if there are no discrepancies
  if (nrow(notUniq_tax.ls[[i]][notUniq_tax.ls[[i]][lastcol] == min(notUniq_tax.ls[[i]][lastcol]),]) == 1) {
    uniq_tax.ls[[i]] <- notUniq_tax.ls[[i]][notUniq_tax.ls[[i]][lastcol] == min(notUniq_tax.ls[[i]][lastcol]),] }
  
  # if there is more than 1 tax assignment/they are different
  if (nrow(notUniq_tax.ls[[i]][notUniq_tax.ls[[i]][lastcol] == min(notUniq_tax.ls[[i]][lastcol]),]) > 1) {
    notUniq_tax.ls[[i]] <- notUniq_tax.ls[[i]][notUniq_tax.ls[[i]][lastcol] == min(notUniq_tax.ls[[i]][lastcol]),] }
}
notUniq_tax.ls <- notUniq_tax.ls[! names(notUniq_tax.ls) %in% names(uniq_tax.ls)]


# add new unique asvs
uniq_tax <- rbind(uniq_tax, do.call(rbind, uniq_tax.ls)[-ncol(do.call(rbind, uniq_tax.ls))])


# blastn these to get most likely
notUniq_tax.ls[[1]]$R1[1] # defluvii
notUniq_tax.ls[[2]]$R1[1] # casseliflavus


# subset just those
notUniq_tax <- do.call(rbind, notUniq_tax.ls)
notUniq_tax <- subset(notUniq_tax, Species %in% c("defluvii", "casseliflavus"))


# combine
uniq_tax <- rbind(uniq_tax, notUniq_tax[-(ncol(uniq_tax) + 1)])
taxa_ww <- merge(fastas_ww, uniq_tax[c(2:8,11)], by = "compare")


# match
taxa_ww <- taxa_ww[order(taxa_ww$ASV),]
taxa_ww <- taxa_ww[order(taxa_ww$WW),]
identical(taxa_ww$WW, colnames(counts_ww))




############
### save ###
############

# reorder
taxa_ww <- taxa_ww[c("UNIQ1", "WW", "compare", "Kingdom", "Phylum", "Class",
                     "Order", "Family", "Genus", "Species", "R1", "R2")]


# add file names
counts_ww <- data.frame(orginalFile = rownames(counts_ww), counts_ww)


# write
write.csv(counts_ww, "./RData/03_counts.csv", row.names = F)
write.csv(taxa_ww, "./RData/03_taxa.csv", row.names = F)



#########
### fasta

fs1 <- counts_ww[-1]
identical(taxa_ww$WW, colnames(fs1))
colnames(fs1) <- taxa_ww$R1

fs2 <- counts_ww[-1]
identical(taxa_ww$WW, colnames(fs2))
colnames(fs2) <- taxa_ww$R2

fs <- cbind(fs1,fs2)
fs <- fs[names(sort(colSums(fs), decreasing = T))]

ids1 <- taxa_ww[c(2,4:11)]
colnames(ids1)[c(1,9)] <- c("ASV", "read")
ids1$ASV <- gsub("WW", "ASV", ids1$ASV)
ids1$dir <- "R1"

ids2 <- taxa_ww[c(2,4:10,12)]
colnames(ids2)[c(1,9)] <- c("ASV", "read")
ids2$ASV <- gsub("WW", "ASV", ids2$ASV)
ids2$dir <- "R2"

ids <- rbind(ids1, ids2)
ids <- ids[match(colnames(fs), ids$read),]

identical(colnames(fs), ids$read)

ids <- paste0(ids$ASV, "__k_",
              ids$Kingdom, "__p_",
              ids$Phylum, "__c_",
              ids$Class, "__o_",
              ids$Order, "__f_",
              ids$Family, "__g_",
              ids$Genus, "__g_",
              ids$Species, "__count_", colSums(fs), ":", 
              ids$dir)


uniquesToFasta(as.matrix(fs), ids = ids,
               fout = "./RData/Full16S_sewageDatabase_ASV.fasta")

