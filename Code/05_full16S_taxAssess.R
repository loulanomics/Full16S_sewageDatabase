#################################################
### Analysis of ASVs and OTUs.
### PacBio Sequel II full 16S rRNA gene amplicons
### Lou LaMartina, finalized 23 March 2022
#################################################


setwd("~/Desktop/Lab/Projects/PacBio/FINAL")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(OTUtable)
library(vegan)
library(scales)


#################
### load data ###
#################

# ASV counts table from 04_full16S_clusterASV.R
counts_asv <- read.csv("./RData/Wastewater_full16S_ASV_counts.csv")
rownames(counts_asv) <- counts_asv$Sample
counts_asv <- counts_asv[-1]


# ASV taxa table from 04_full16S_clusterASV.R
taxa_asv <- read.csv("./RData/Wastewater_full16S_ASV_taxonomy.csv")


# OTU counts table from 04_full16S_clusterASV.R
counts_otu <- read.csv("./RData/Wastewater_full16S_OTU_counts.csv")
rownames(counts_otu) <- counts_otu$Sample
counts_otu <- counts_otu[-1]


# ASV taxa table from 04_full16S_clusterASV.R
taxa_otu <- read.csv("./RData/Wastewater_full16S_OTU_taxonomy.csv")


# sample info
info <- read.csv("./RData/Wastewater_full16S_sample_metadata.csv")


# format date
info$Date <- as.Date(info$Date, format = "%m/%d/%y")




#############################
### lowest classification ###
#############################

# get lowest level of taxonomic classifications for each OTU

# create vector of classifications
classes <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# make sure columns are in this order:
taxa_otu <- taxa_otu[c("OTU_prop", classes)]


# change unclassified to blank
taxa_otu[taxa_otu == "unclassified"] <- NA



###################################
### OTUs unclassified at each level

# empty list
unknownOTUs.ls <- list()

# for each taxa class/column in taxa data frame,
for (i in classes){

  # in the "kingdom" column,
  if (i == "Kingdom") {

    # if it's blank, obtain its OTU names -> have list of all OTUs
    # that are not classified to kingdom level.
    unknownOTUs.ls[[i]] <- taxa_otu[is.na(taxa_otu[[i]]), "OTU_prop"]

    # in the other columns,
  } else if (i != "Kingdom") {

    # if they are blank, obtain their OTU names, but don't include ones found previously.
    unknownOTUs.ls[[i]] <- setdiff(taxa_otu[is.na(taxa_otu[[i]]), "OTU_prop"],
                                   taxa_otu[is.na(taxa_otu[[match(i, classes)]]), "OTU_prop"])
  }
}


taxa_otu$OTU_prop[which(is.na(taxa_otu[[i]]))] %in% classes.ls[[j]]$OTU_prop



#################################
### make names their lowest class

# empty list
classes.ls <- list()

# for each column in tax data frame,
for (i in classes) {

  j = classes[match(i, classes) - 1]
  props = subset(taxa_otu, OTU_prop %in% unknownOTUs.ls[[i]])$OTU_prop

  # in kingdom column, if it has unclassified OTUs,
  if (i == "Kingdom" & length(props) > 0) {

    # if they are blank, call them "unclassified"
    classes.ls[[i]] <- data.frame(OTU_prop = props,
                                  Name = "unclassified") }

  # for all the others, if it has unclassified OTUs,
  if (i != "Kingdom" & length(props) > 0) {

    # call them by their lowest level of classification = the column before it.
    classes.ls[[i]] <- data.frame(OTU_prop = props, Name = paste0(j, "__",
                                                                  subset(taxa_otu, OTU_prop %in% unknownOTUs.ls[[i]])[[match(i, classes)]]))
  }
}



##########################
### add names to otu table

# combine into single data frame
classes.df <- dplyr::bind_rows(classes.ls)



# for those classified down to lowest level, give them those names
low = classes[length(classes)]
classes.df <- rbind(classes.df,
                    data.frame(OTU_prop = taxa_otu[is.na(taxa_otu[low]) == F, "OTU_prop"],
                               Name = paste0(low, "__", taxa_otu[is.na(taxa_otu[low]) == F, low])))


# remove characters R doesn't like
classes.df$Name <- gsub("-", "_", classes.df$Name)
classes.df$Name <- gsub("\\(", "_", classes.df$Name)
classes.df$Name <- gsub("\\)", "", classes.df$Name)
classes.df$Name <- gsub(" ", "_", classes.df$Name)


# add to taxa_otu table
taxa_otu <- merge(classes.df, taxa_otu, by = "OTU_prop")




################
### sum OTUs ###
################

# need variable for those that are the same species but different higher classes
dup_names <- unique(taxa_otu[-1])
dup_names <- dup_names$Name[dup_names$Name %in% names(which(table(dup_names$Name) > 1))]
dup_names$Name <- paste0("Genus__", dup_names$Genus, "__", dup_names$Species)

taxa_otu$Name[taxa_otu$Name %in% dup_names] <- 
  paste0("Genus__", taxa_otu$Genus[taxa_otu$Name %in% dup_names], "__",
         taxa_otu$Species[taxa_otu$Name %in% dup_names])


# get OTUs for each class
nameOTUs.ls <- list()
for (i in unique(taxa_otu$Name)) {
  nameOTUs.ls[[i]] <- as.character(taxa_otu$OTU_prop[taxa_otu$Name == i])
}


# sum those OTUs in each sample
counts_name.ls <- list()
for( i in unique(taxa_otu$Name)) {
  counts_name.ls[[i]] <- rowSums(counts_otu[colnames(counts_otu) %in% unlist(nameOTUs.ls[i])])
}


# make new counts table - "Name" as column instead of OTUs
counts_name <- data.frame(do.call(cbind, counts_name.ls))


# order by abundance
counts_name <- counts_name[names(sort(colSums(counts_name), decreasing = T))]


# convert to relative abundance
relabun_name <- counts_name / rowSums(counts_name)




#######################
### subset datasets ###
#######################

# samples diverse spatially & temporally were purposefully chosen.
# want to compare those

# keep relevant info
set_info <- info[1:7]


# get spring & fall samples, and those in "transition period" (from lamartina 2021)
fall_smps <- subset(set_info, City == "Milwaukee")
fall_smps <- fall_smps$Sample[format(fall_smps$Date, "%b") %in% c("Aug", "Sep", "Oct", "Nov")]
fall_smps

spring_smps <- subset(set_info, City == "Milwaukee")
spring_smps <- spring_smps$Sample[format(spring_smps$Date, "%b") %in% c("Feb", "Mar", "Apr", "May")]
spring_smps

trans_smps <- subset(set_info, City == "Milwaukee")
trans_smps <- trans_smps$Sample[format(trans_smps$Date, "%b") %in% c("Dec", "Jan", "Jun", "Jul")]
trans_smps


# get samples from north & south US
south_smps <- subset(set_info, City != "Milwaukee")
south_smps <- south_smps[order(south_smps$Temp_profile_metric, decreasing = T),]
south_smps <- south_smps$Sample[1:11]
south_smps

north_smps <- subset(set_info, City != "Milwaukee")
north_smps <- north_smps[order(north_smps$Temp_profile_metric, decreasing = F),]
north_smps <- north_smps$Sample[1:11]
north_smps


# add variables to set_info
set_info$Set[set_info$Sample %in% fall_smps] <- "Fall"
set_info$Set[set_info$Sample %in% spring_smps] <- "Spring"
set_info$Set[set_info$Sample %in% trans_smps] <- "Transition"

set_info$Set[set_info$Sample %in% north_smps] <- "North"
set_info$Set[set_info$Sample %in% south_smps] <- "South"




################
### top OTUs ###
################

# get most abundant groups in each category
n = 6
fall_top <- counts_name[rownames(counts_name) %in% set_info$Sample[set_info$Set == "Fall"],]
fall_top <- fall_top[names(sort(colSums(fall_top), decreasing = T))][1:n]
fall_top
fall_top <- colnames(fall_top)

spring_top <- counts_name[rownames(counts_name) %in% set_info$Sample[set_info$Set == "Spring"],]
spring_top <- spring_top[names(sort(colSums(spring_top), decreasing = T))][1:n]
spring_top
spring_top <- colnames(spring_top)

trans_top <- counts_name[rownames(counts_name) %in% set_info$Sample[set_info$Set == "Transition"],]
trans_top <- trans_top[names(sort(colSums(trans_top), decreasing = T))][1:n]
trans_top
trans_top <- colnames(trans_top)

north_top <- counts_name[rownames(counts_name) %in% set_info$Sample[set_info$Set == "North"],]
north_top <- north_top[names(sort(colSums(north_top), decreasing = T))][1:n]
north_top
north_top <- colnames(north_top)

south_top <- counts_name[rownames(counts_name) %in% set_info$Sample[set_info$Set == "South"],]
south_top <- south_top[names(sort(colSums(south_top), decreasing = T))][1:n]
south_top
south_top <- colnames(south_top)


# get uniques of each
tops <- unique(c(fall_top, spring_top, trans_top, north_top, south_top))
tops


# add top names, call everything else "other"
taxa_tops <- unique(data.frame(Name = taxa_otu$Name, Top = taxa_otu$Name))
taxa_tops$Top[! taxa_tops$Top %in% tops] <- "Other"


# sum counts in sets
counts_tops <- data.frame(Other = rowSums(counts_name[! colnames(counts_name) %in% tops]), counts_name)
counts_tops <- counts_tops[colnames(counts_tops) %in% c(tops, "Other")]
counts_tops <- data.frame(Sample = rownames(counts_tops), counts_tops)
counts_tops <- merge(set_info[c(1,8)], counts_tops, by = "Sample")[-1]
counts_tops <- aggregate(. ~ Set, sum, data = counts_tops)




#########################
### stacked bar chart ###
#########################

# order taxa, make labels
tops <- names(sort(colSums(counts_tops[colnames(counts_tops) %in% tops]), decreasing = T))
taxlabels <- data.frame(Name = tops, Order.f = paste(sprintf("%02d", 1:length(tops)), tops))
taxlabels <- rbind(taxlabels, c("Other", paste(length(tops) + 1, "Other")))


# order axis by cloacibacterium abundance ("warm" indicator)
setlabels <- counts_tops$Set[order(counts_tops$Genus__Cloacibacterium, decreasing = T)]
setlabels <- data.frame(Set = setlabels, Order.x = paste(1:length(setlabels), setlabels))
setlabels$Lab <- c("South\nUSA", "North\nUSA", "Fall\nperiod", 
                   "Transition\nperiod", "Spring\nperiod")
labs.x <- setlabels$Lab
names(labs.x) <- setlabels$Order.x


# melt, add labels
tops.m <- melt(counts_tops, variable.name = "Name", value.name = "Count")
tops.m <- merge(tops.m, taxlabels, by = "Name")
tops.m <- merge(tops.m, setlabels, by = "Set")


# legend labels
taxlabels <- unique(merge(taxlabels, taxa_otu[-1], by = "Name", all.x = T))
taxlabels$Lab <- taxlabels$Genus
taxlabels$Lab[is.na(taxlabels$Species) == F] <- paste0(taxlabels$Genus[is.na(taxlabels$Species) == F], " (",
                                                       taxlabels$Species[is.na(taxlabels$Species) == F], ")")
taxlabels$Lab[taxlabels$Name == "Other"] <- "Other"
labs.f <- taxlabels$Lab
names(labs.f) <- taxlabels$Order.f


# plot
bar.plot <-
  ggplot(tops.m, aes(x = Order.x, y = Count, fill = Order.f)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(labels = labs.f,
                    values = c(brewer.pal(11, "Paired"), 
                               "gold3", "#ceede3", "darkslategray4", "grey90")) +
  scale_x_discrete(labels = labs.x) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black", vjust = 8),
        axis.text.y = element_text(size = 10, color = "black", margin = unit(c(0,-0.25,0,0), "cm")),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold", margin = unit(c(-0.4,0,0,0), "cm")),
        legend.title = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.background = element_blank(),
        legend.box.margin = margin(c(0, 0, 0, -15)),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  guides(fill = guide_legend(keyheight = 1.25, keywidth = 0.75, units = "in", ncol = 1)) +
  labs(y = "Taxa proportions", x = "Sample set", fill = "Common taxa:\nGenus (species)")
bar.plot

ggsave("./Plots/barplot.pdf", plot = bar.plot, device = "pdf", width = 6, height = 6, units = "in")




##########################
### dendrogram heatmap ###
##########################

# taxa found in at least 6 samples
commons <- list()
for (i in colnames(counts_name)) {
  n = length(which(counts_name[[i]] > 0))
  if (n > 5) {
    commons[[i]] <- i 
  }
}
commons <- unlist(commons)


# relative abundance samples = rows & cols equal 1
dendro_abun <- relabun_name[colnames(relabun_name) %in% commons]


# convert to z score  = (x-μ)/σ
dendro_abun <- t(zscore(t(dendro_abun)))


# euclidian distances of OTUs & samples
otu.dist <- vegdist(t(dendro_abun), method = "euclidian")
smp.dist <- vegdist(dendro_abun, method = "euclidian")


# cluster based on euclidian distances
otu.clus <- hclust(otu.dist, method = "average")
smp.clus <- hclust(smp.dist, method = "average")


# extract order of dendro leaves
smp_order = smp.clus$labels[c(smp.clus$order)]
otu_order = otu.clus$labels[c(otu.clus$order)]


# make same heatmap in ggplot2, so i can actually modify it
dendro.df <- data.frame(Sample = rownames(dendro_abun), dendro_abun)


# order samples
dendro.df <- data.frame(smp_order = paste(sprintf("%02d", 1:length(smp_order)), 
                                       smp_order), dendro.df)


# melt for plotting
dendro.m <- melt(dendro.df, variable.name = "OTU", value.name = "Abundance")


# sample labels
smplabels <- data.frame(Sample = smp_order)
smplabels <- merge(smplabels, set_info, by = "Sample")
smplabels$Lab <- paste0(smplabels$Set, " (", smplabels$State, ", ", format(smplabels$Date, "%b"), ")")


# reorder
smplabels <- smplabels[match(smp_order, smplabels$Sample),]
labs.y <- smplabels$Lab
names(labs.y) <- smp_order


# legend labels
taxlabels <- data.frame(Name = otu_order)
taxlabels <- unique(merge(taxlabels, taxa_otu[-1], by = "Name"))

taxlabels$Lab[is.na(taxlabels$Species) == F] <- 
  paste0(taxlabels$Phylum[is.na(taxlabels$Species) == F], " (",
         taxlabels$Genus[is.na(taxlabels$Species) == F], " ",
         taxlabels$Species[is.na(taxlabels$Species) == F], ")")

taxlabels$Lab[is.na(taxlabels$Genus) == F & is.na(taxlabels$Lab)] <- 
  paste0(taxlabels$Phylum[is.na(taxlabels$Genus) == F & is.na(taxlabels$Lab)], " (",
         taxlabels$Genus[is.na(taxlabels$Genus) == F & is.na(taxlabels$Lab)], ")")

taxlabels$Lab[is.na(taxlabels$Family) == F & is.na(taxlabels$Lab)] <- 
  paste0(taxlabels$Phylum[is.na(taxlabels$Family) == F & is.na(taxlabels$Lab)], " (",
         taxlabels$Family[is.na(taxlabels$Family) == F & is.na(taxlabels$Lab)], ")")

taxlabels$Lab[is.na(taxlabels$Order) == F & is.na(taxlabels$Lab)] <- 
  paste0(taxlabels$Phylum[is.na(taxlabels$Order) == F & is.na(taxlabels$Lab)], " (",
         taxlabels$Order[is.na(taxlabels$Order) == F & is.na(taxlabels$Lab)], ")")


# change this one
taxlabels$Lab <- gsub("T34", "Burkholderiales", taxlabels$Lab)


# reorder
taxlabels <- taxlabels[match(otu_order, taxlabels$Name),]
labs.x <- taxlabels$Lab
names(labs.x) <- otu_order


# distribute colors in heat map
heatvals <- sort(c(0, as.numeric(quantile(dendro.m$Abundance, probs = seq(0, 1, 0.1)))))
heatvals <- rescale(heatvals)
heatcols <- c(colorRampPalette(brewer.pal(9, "Blues")[1:3])(10), brewer.pal(9, "Blues")[c(5,9)])


# plot
heat.plot <-
  ggplot(dendro.m, aes(y = Sample, x = OTU, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradientn(colors = heatcols, values = heatvals) +
  scale_y_discrete(limits = smp_order, labels = labs.y, position = "right", ) +
  scale_x_discrete(limits = otu_order, labels = labs.x) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = "black", angle = 90, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "bottom",
        legend.title.align = 0.5,
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(direction = "horizontal", title.position = "top",
                               barwidth = 10, barheight = 0.5, units = "in")) +
  labs(x = "Common unique taxa (OTU):\nPhylum (lowest classification)", fill = "Z score",
       y = "Sample set: (state, month)")
heat.plot

ggsave("./Plots/heat.pdf", plot = heat.plot, device = "pdf", width = 10, height = 10, units = "in")


# save for dendrograms
heatmap(as.matrix(dendro_abun),
        Rowv = as.dendrogram(smp.clus),
        Colv = as.dendrogram(otu.clus),
        margins = c(12, 5))
