#################################################
### Cluster ASVs to OTUs with 99.5% similarity
### using mothur pairwise.seqs() and cluster().
### PacBio Sequel II full 16S rRNA gene amplicons
### Lou LaMartina, finalized 23 March 2022
#################################################


setwd("~/Desktop/Lab/Projects/PacBio/FINAL")
library(ggplot2)


#################
### load data ###
#################

# ASV table from 03_full16S_subSewage.R
counts_ww <- read.csv("./RData/03_counts.csv")
rownames(counts_ww) <- counts_ww$orginalFile
counts_ww <- counts_ww[-1]


# taxa table from 03_full16S_subSewage.R
taxa_ww <- read.csv("./RData/03_taxa.csv")


# sample info
info <- read.csv("./RData/Wastewater_full16S_sample_metadata.csv")




######################
### cluster to OTU ###
######################
# see 03_full16S_subSewage.R

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# in terminal

### execute mothur
# $ ~/mothur/mothur

### calculate pairwise distances -
# mothur > pairwise.seqs(fasta=Wastewater_full16S_23Mar22.fasta)
# It took 1164 secs to find distances for 1041 sequences. 301252 distances below cutoff 1.

### assign sequences to OTUs
# mothur > cluster(column=Wastewater_full16S_23Mar22.dist, name=names, cutoff=0.005)
# It took 0 seconds to cluster

# # # # # # # # # # # # # # # # # # # # # # # # # # #


# open cluster results
clusters <- read.table("./mothur/Wastewater_full16S_23Mar22.opti_mcc.list", sep = "\t")
clusters <- data.frame(t(clusters[-c(1,2)]))
colnames(clusters) <- c("OTU", "UNIQ_list")
clusters$OTU <- gsub("Otu", "OTU", clusters$OTU)


# turn ASVs into list
clusters.ls <- list()
for (i in clusters$OTU) {
  clusters.ls[[i]] <- unlist(strsplit(clusters$UNIQ_list[clusters$OTU == i], ","))
}




#####################
### glom taxonomy ###
#####################

# add OTU variable to corresponding ASVs
taxa_ww$OTU <- NA

for (i in clusters$OTU) {
  for (j in 1:length(clusters.ls[[i]]))
    taxa_ww$OTU[taxa_ww$UNIQ1 == clusters.ls[[i]][j]] <- i
}


# find shared taxonomy among OTUs
taxa_otu <- unique(taxa_ww[c("OTU", "Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus", "Species")])
nrow(taxa_otu) # [1] 693
dups_otu <- names(which(table(taxa_otu$OTU) > 1))
length(dups_otu)
# [1] 17


# what species makes them distinct?
dups_otu.ls <- list()
for (i in dups_otu) {
  dups_otu.ls[[i]] <- unique(taxa_ww$Species[taxa_ww$OTU == i])
}


# get proportions of species in duplicated OTUs
props.ls <- list()
for (otu in dups_otu) {
  species <- unique(taxa_ww$Species[taxa_ww$OTU == otu])
  props.ls[[otu]] <- otu
  for (n in 1:length(species)) {
    otu.asvs <- taxa_ww$WW[taxa_ww$OTU == otu]
    sum.otu.asvs <- sum(counts_ww[colnames(counts_ww) %in% otu.asvs])
    otu.asv.specs <- taxa_ww$WW[taxa_ww$OTU == otu & taxa_ww$Species == species[n]]
    sum.otu.asv.specs <- sum(counts_ww[colnames(counts_ww) %in% otu.asv.specs])
    props.ls[[otu]][n] <- sum.otu.asv.specs / sum.otu.asvs
  }
}


# convert to data frame
props <- data.frame(OTU = substr(names(unlist(props.ls)), 1,6),
                    Species = unlist(dups_otu.ls),
                    prop = as.numeric(unlist(props.ls)))
props$prop <- sapply(props$prop, function(n) signif(n, 3))


# make unique OTU variable
dups_otu.df <- taxa_otu[taxa_otu$OTU %in% dups_otu,]
dups_otu.df$OTU_species <- paste0(dups_otu.df$OTU, "_", dups_otu.df$Species)
props$OTU_species <- paste0(props$OTU, "_", props$Species)
dups_otu.df <- merge(dups_otu.df, props[3:4], by = "OTU_species")
dups_otu.df$OTU_prop <- paste0(dups_otu.df$OTU, "_",
                               gsub("0\\.", "", sub("0+$", "", as.character(sapply(dups_otu.df$prop, function(n) signif(n, 2))))))


# add to OTU taxa table
taxa_noDup <- taxa_otu[! taxa_otu$OTU %in% dups_otu,]
taxa_noDup$prop <- 1
taxa_noDup$OTU_prop <- paste0(taxa_noDup$OTU, "_", taxa_noDup$prop)
taxa_noDup$OTU_species <- paste0(taxa_noDup$OTU, "_", taxa_noDup$Species)
taxa_noDup <- taxa_noDup[colnames(dups_otu.df)]
taxa_noDup <- rbind(taxa_noDup, dups_otu.df)
names(which(table(taxa_noDup$OTU_prop) > 1))


# now regroup ASVs in OTUs, given proportions of species assignments
taxa_otu$OTU_species <- paste0(taxa_otu$OTU, "_", taxa_otu$Species)
taxa_otu <- merge(taxa_noDup[c(1,11)], taxa_otu, by = "OTU_species")

taxa_ww$OTU_species <- paste0(taxa_ww$OTU, "_", taxa_ww$Species)
taxa_asv <- merge(taxa_noDup[c(1,11)], taxa_ww, by = "OTU_species")




###################
### glom counts ###
###################

# get proportions of species in duplicated OTUs
asvs.ls <- list()
for (i in unique(taxa_otu$OTU_prop)) {
  asvs.ls[[i]] <- taxa_asv$WW[taxa_asv$OTU_prop == i]
}


# sum samples within otus
counts_otu.ls <- list()
for (i in unique(taxa_otu$OTU_prop)) {
  counts_otu.ls[[i]] <- rowSums(counts_ww[colnames(counts_ww) %in% as.character(unlist(asvs.ls[i]))])
}
counts_otu <- data.frame(do.call(cbind, counts_otu.ls))




############
### save ###
############

# reorder least to greatest
counts_otu <- counts_otu[names(sort(colSums(counts_otu), decreasing = T))]
taxa_otu <- taxa_otu[match(colnames(counts_otu), taxa_otu$OTU_prop),]
identical(taxa_otu$OTU_prop, colnames(counts_otu))

taxa_asv <- taxa_asv[order(taxa_asv$WW),]
identical(taxa_asv$WW, colnames(counts_ww))


# new sample names
info <- info[order(info$originalName),]
counts_otu <- counts_otu[order(rownames(counts_otu)),]
identical(info$originalName, rownames(counts_otu))
rownames(counts_otu) <- info$Sample

counts_ww <- counts_ww[order(rownames(counts_ww)),]
identical(info$originalName, rownames(counts_ww))
rownames(counts_ww) <- info$Sample


# add samples names
counts_otu <- data.frame(Sample = rownames(counts_otu), counts_otu)
counts_asv <- data.frame(Sample = rownames(counts_ww), counts_ww)


# reorder
taxa_asv <- taxa_asv[c("WW", "OTU", "OTU_prop", "Kingdom", "Phylum", "Class",
                       "Order", "Family", "Genus", "Species", "R1", "R2")]
taxa_otu <- taxa_otu[c("OTU", "OTU_prop", "Kingdom", "Phylum", "Class",
                       "Order", "Family", "Genus", "Species")]

colnames(counts_asv) <- gsub("WW", "ASV", colnames(counts_asv))
taxa_asv$WW <- gsub("WW", "ASV", taxa_asv$WW)
colnames(taxa_asv)[1] <- "ASV"


# save
write.csv(counts_otu, "./RData/Wastewater_full16S_OTU_counts.csv", row.names = F, na = "")
write.csv(counts_asv, "./RData/Wastewater_full16S_ASV_counts.csv", row.names = F, na = "")
write.csv(taxa_otu, "./RData/Wastewater_full16S_OTU_taxonomy.csv", row.names = F, na = "")
write.csv(taxa_asv, "./RData/Wastewater_full16S_ASV_taxonomy.csv", row.names = F, na = "")
