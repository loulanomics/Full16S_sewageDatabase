#################################################
### Remove short reads and trim off residual
### primer sequences from unique reads.
### PacBio Sequel II full 16S rRNA gene amplicons
### Lou LaMartina, finalized 22 March 2022
#################################################


setwd("~/Desktop/Lab/Projects/PacBio/FINAL")
library(ggplot2)


#################
### load data ###
#################

# ASV table from galaxy
counts_orig <- read.table("./RData/Full16S_dada2_sequencetable.txt", sep = "", header = T)


# simplify names
colnames(counts_orig) <- gsub("\\..", "__", gsub(".fastq.gz", "", gsub("trim.", "", colnames(counts_orig))))


# taxonomy table
taxa_orig <- read.table("./RData/Full16S_assignTaxonomy_addSpecies.txt", sep = "\t", header = T)
colnames(taxa_orig)[1] <- "fasta"




##########################
### remove short reads ###
##########################

# inspect sequence length distribution -
# given Ecoli position 1492 is the position of the reverse primer,
# and 27 is the position of the forward primer, we expect lengths to
# be around 1400 bp, with room for  variability among taxonomic groups.
# we also are not expecting/setting a maximum length.

# get bp lengths of ASVs
seq_distribution <- data.frame(table(nchar(rownames(counts_orig)))) 
colnames(seq_distribution) <- c("SeqLength", "Frequency")
seq_distribution$SeqLength <- as.numeric(as.character(seq_distribution$SeqLength))


# minimum sequence length = 5% below the median length
# maximum length = the longest seq in the dataset (or none, in other words)
medianLength <- median(nchar(rownames(counts_orig)))
minLength <- floor(medianLength * 0.95)
maxLength <- max(seq_distribution$SeqLength)


# remove those, save into new data frame
counts_untrim <- data.frame(t(counts_orig[nchar(rownames(counts_orig)) %in% seq(minLength, maxLength), ]))
counts_untrim <- counts_untrim[,names(sort(colSums(counts_untrim), decreasing = T))]


# save fasta sequences after short reads removed, make new ASV variable called "ORIG"
fastas_untrim <- data.frame(ORIG = paste0("UNTRIM", sprintf("%04d", 1:ncol(counts_untrim))),
                          fasta = colnames(counts_untrim))


# match taxa table to new ASVs
taxa_untrim <- subset(taxa_orig, fasta %in% fastas_untrim$fasta)


# calculate loss when short reads removed
nrow(counts_orig)     # [1] 5974
sum(counts_orig)      # [1] 2,676,766 total reads
ncol(counts_untrim)          # [1] 5064 ASVs
sum(counts_untrim)           # [1] 1,048,998 reads after trimming
(sum(counts_untrim) - sum(counts_orig)) / sum(counts_orig) * 100 
# [1] -60.81099 % lost


# visualize length distribution
dist.plot <-
  ggplot(seq_distribution, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black") +
  geom_text(data = data.frame(x = medianLength,
                              y = max(seq_distribution$Frequency) - 10),
            aes(x, y, label = paste0("median\n", medianLength)), color = "red") +
  geom_text(data = data.frame(x = minLength, y = 10),
            aes(x, y, label = paste0("min\n(median - 5%)\n", minLength)), color = "red") +
  geom_text(data = data.frame(x = maxLength, y = 10),
            aes(x, y, label = paste0("max\n(longest read)\n", maxLength)), color = "red") +
  scale_x_continuous(breaks = ceiling(quantile(seq_distribution$SeqLength)),
                     labels = ceiling(quantile(seq_distribution$SeqLength)))  +
  labs(x = "Sequence length (bp; quantiles labelled)", y = "Frequency")
dist.plot

ggsave("./Plots/seq_distribution.pdf", plot = dist.plot, device = "pdf", width = 16, height = 4, units = "in")




#######################
### convert primers ###
#######################

# we need to find and remove primers that cutadapt did not get.
# here we are converting ambiguous IUPAC bases to all their possibilities,
# since ambiguous bases are not in real reads.
# https://droog.gs.washington.edu/parc/images/iupac.html


# data frame of IUPAC nucleotides -
# for grep pattern matching later, the pipe | considers all of those
# possibilities at that position when inside brackets
iupac <- data.frame(Code = c("A", "T", "C", "G", 
                             "M", "R", "W", "S", "Y", "K", 
                             "V", "H", "D", "B", "N"),
                    Meaning = c("A", "T", "C", "G", 
                                "A|C", "A|G", "A|T", "C|G" , "C|T", "G|T", 
                                "A|C|G", "A|C|T", "A|G|T", "C|G|T", "G|A|T|C"),
                    Complement = c("T", "A", "G", "C", 
                                   "K", "Y", "W", "S", "R", "M", 
                                   "B", "D", "H", "V", "N"))


# subset only ambiguous ones (ie. not ATCG)
ambs <- iupac$Code[5:15]
names(ambs) <- iupac$Meaning[5:15]


# primers that need to be removed
F27 = "AGRGTTYGATYMTGGCTCAG"
R1492 = "RGYTACCTTGTTACGACTT"


# forward primer reverse complements
F1492 <- list()
for (i in 1 : length(unlist(strsplit(R1492, "")))) {
  primer.split <- unlist(strsplit(R1492, ""))
  F1492[[i]] <- iupac$Complement[iupac$Code == primer.split[[i]]]
}
F1492 = paste(rev(F1492), collapse = "")


# reverse primer reverse complements
R27 <- list()
for (i in 1 : length(unlist(strsplit(F27, "")))) {
  primer.split <- unlist(strsplit(F27, ""))
  R27[[i]] <- iupac$Complement[iupac$Code == primer.split[[i]]]
}
R27 = paste(rev(R27), collapse = "")


# combine into data frame
primers = data.frame(Primer = c("F27", "R27", "R1492", "F1492"),
                     Seq = c(F27, R27, R1492, F1492))




#######################
### partial primers ###
#######################

# i think cutadapt left many primers because there are many
# that are only partial. cutadapt can usually handle this
# but there are also many deletions. in future, increase
# allowed number of errors in the 3' primer.


# get partial primers, down to 9 bp length -
# this creates a vector of new primers, all with one less bp.
# as in, the first 9 bases of F primer, first 10, first 11....
primers.ls <- list()
for (i in primers$Seq) {
  for (j in 0 : (nchar(i) - 11)) {
    primers.ls[[i]] [j + 1] <- substr(i, 1, nchar(i) - j)
  }
}


# convert list to data frame - remove number added into names when listing
primers.df <- data.frame(Seq = gsub("[[:digit:]]+", "", names(unlist(primers.ls))),
                         Sub = unlist(primers.ls))


# add to original primer list
primers.df <- merge(primers, primers.df, by = "Seq")


# create a variable showing which primer it is and how many bp removed
for (i in primers$Primer) {
  primers.df$Primer_sub[primers.df$Primer == i] <- 
    paste0(i, "_", 1:length(unique(primers.df$Sub[primers$Primer == i])) - 1)
}
   

# expand ambiguous bases
primers.df$Ambig <- primers.df$Sub
for (i in 1:length(primers.df$Ambig)) {
  primer.split <- unlist(strsplit(primers.df$Ambig[i], ""))
  ambs.match <- primer.split[primer.split %in% ambs]
  ambs[ambs %in% ambs.match]
  
  # for primers with ambiguous bases,
  if(length(ambs.match) > 0) {
    for (j in 1:length(ambs.match)) {
      cat("\n", primers.df$Primer_sub[i], ":", primers.df$Ambig[i], "->", ambs.match[j], "-> ")
      
      # in brackets, expand them to their possibilties
      primers.df$Ambig[i] <- gsub(as.character(ambs[ambs %in% ambs.match[j]]),
                         paste0("[", names(ambs[ambs %in% ambs.match[j]]), "]"),
                         primers.df$Ambig[i])
      cat(primers.df$Ambig[i]) }
    
    # if no ambiguities, skip
  } else {
    cat("\n", primers.df$Primer[i], ":", primers.df$Ambig[i], "-> no ambiguities")
  }
}




#####################
### match primers ###
#####################

# find primers in sequences using pattern matching grep
matches.ls <- list()
for (i in fastas_untrim$fasta) {
  for (j in 1:length(primers.df$Ambig)) {
    
    # for sequences with primers, 
    if(grepl(primers.df$Ambig[j], i) == TRUE) {
      
      # which primer is it, what position in the read is it, and what was the sequence
      matches.ls[[i]] <- c(primers.df$Primer[j], primers.df$Primer_sub[j], 
                         unlist(gregexpr(primers.df$Ambig[j], i)), i)
    }
  }
}


# convert list to data frame
matches <- data.frame(do.call(rbind, matches.ls))
colnames(matches) <- c("Primer", "Primer_sub", "Position", "fasta")
matches$Position <- as.numeric(matches$Position)


# add sequences, including those without primers found
matches <- merge(fastas_untrim, matches, by = "fasta", all = T)


# which partial primers were found most?
table(matches$Primer_sub)
# F1492_9  F27_10 R1492_9  R27_10 
#    2837       2      14    2176 


# what position are they found at most?
aggregate(Position ~ Primer_sub, summary, data = matches)
match.stats <- aggregate(Position ~ Primer_sub, mean, data = matches)
# Primer_sub        Min.        1st Qu.         Median        Mean        3rd Qu.        Max.
#    F1492_9       1379           1440            1459        1500           1479        3288
#     F27_10        759            759             759         759            759         759
#    R1492_9          1           1117            1122        1367           2027        2027
#     R27_10        708           1440           1460         1487           1479        3093


# visualize the frequency of primers at what positions -
# i don't want to remove any "real" sequences
match.plot <-
  ggplot(matches[complete.cases(matches),], aes(x = Primer_sub, y = Position)) +
  geom_boxplot() + geom_jitter(width = 0.25) +
  geom_text(data = match.stats, color = "red",
            label = paste0("mean\n", ceiling(match.stats$Position)))
match.plot

ggsave("./Plots/primer_positions.pdf", plot = match.plot, device = "pdf", width = 6, height = 4, units = "in")


# choosing only the revcomps and only primers at 1300+
# ... the few that are before 1300 are uncommon and could be real sequences
matches_trim <- matches[matches$Primer %in% c("F1492", "R27") & matches$Position > 1300,]


# trim at those positions
matches_trim$R1 <- NA
for (i in 1:nrow(matches_trim)) {
  matches_trim$R1[i] <- substr(matches_trim$fasta[i], 1, matches_trim$Position[i] - 1)
}
unique(nchar(matches_trim$R1) == matches_trim$Position - 1)


# replace sequences with trimmed ones in new column "R1"
matches$R1 <- matches$fasta
matches <- matches[! matches$fasta %in% matches_trim$fasta,]
matches <- rbind(matches_trim, matches)


# just checking that worked - any that are "false" were trimmed
# because the original fasta sequence does not match the new R1 column
matches$notTrimmed <- matches$fasta == matches$R1




#################################
### dereplicate trimmed reads ###
#################################

################
### counts table

# find ASVs with same unique sequence
uniqs.ls <- list()
for (i in unique(matches$R1)) {
  uniqs.ls[[i]] <- matches$ORIG[matches$R1 == i]
}


# change counts column names to ASV
counts_trim <- counts_untrim
matches <- matches[order(matches$ORIG),]
identical(matches$fasta, colnames(counts_trim))
colnames(counts_trim) <- matches$ORIG


# sum counts within unique sequences
counts_uniq.ls <- list()
for (i in unique(matches$R1)) {
  counts_uniq.ls[[i]] <- rowSums(counts_trim[colnames(counts_trim) %in% as.character(unlist(uniqs.ls[i]))])
}
counts_uniq <- data.frame(do.call(cbind, counts_uniq.ls))


# sort by abundance
counts_uniq <- counts_uniq[,names(sort(colSums(counts_uniq), decreasing = T))]


# make new asv variable called "UNIQ"
fastas_uniq <- data.frame(UNIQ = paste0("UNIQ", sprintf("%04d", 1:ncol(counts_uniq))),
                          R1 = colnames(counts_uniq))


# change column names to new ASV variable
identical(fastas_uniq$R1, colnames(counts_uniq))
colnames(counts_uniq) <- fastas_uniq$UNIQ



##############
### taxa table

# change taxa fastas to trimmed versions
taxa_dups <- taxa_orig[match(matches$fasta, taxa_orig$fasta),]
identical(matches$fasta, taxa_dups$fasta)
taxa_dups <- cbind(matches[c(2,6)], taxa_dups[-1])
taxa_dups[is.na(taxa_dups)] <- "unclassified"


# condense to unique ASVs
noDups_tax.ls <- list()
dups_tax.ls <- list()
for (i in names(uniqs.ls)) {
  
  # if there are no discrepancies in tax assignments
  if(nrow(unique(subset(taxa_dups[-1], R1 %in% i))) == 1) {
    noDups_tax.ls[[i]] <- unique(subset(taxa_dups[-1], R1 %in% i)) }
  
  # if there is more than 1 tax assignment/they are different
  if(nrow(unique(subset(taxa_dups[-1], R1 %in% i))) > 1) {
    dups_tax.ls[[i]] <- unique(subset(taxa_dups[-1], R1 %in% i)) }
}


# convert to data frame
noDups_tax <- data.frame(do.call(rbind, noDups_tax.ls))


# with ones that have discrepancies
for (i in names(dups_tax.ls)) {
  
  # add new column
  n = ncol(dups_tax.ls[[i]]) + 1
  
  # for each taxon assigned to a given asv
  for (j in 1:nrow(data.frame(dups_tax.ls[[i]]))) {
    
    # count how many instances of "unclassified" there are
    dups_tax.ls[[i]][j, n] <- length(grep("unclassified", data.frame(dups_tax.ls[[i]])[j,]))
  }
}


# choose asvs with least "unclassifieds"
noDups_tax.ls <- list()
for (i in names(dups_tax.ls)) {
  
  # if there are no discrepancies
  if (nrow(dups_tax.ls[[i]][dups_tax.ls[[i]]$V9 == min(dups_tax.ls[[i]]$V9),]) == 1) {
    noDups_tax.ls[[i]] <- dups_tax.ls[[i]][dups_tax.ls[[i]]$V9 == min(dups_tax.ls[[i]]$V9),] }
  
  # if there is more than 1 tax assignment/they are different
  if (nrow(dups_tax.ls[[i]][dups_tax.ls[[i]]$V9 == min(dups_tax.ls[[i]]$V9),]) > 1) {
    dups_tax.ls[[i]] <- dups_tax.ls[[i]][dups_tax.ls[[i]]$V9 == min(dups_tax.ls[[i]]$V9),] }
}
dups_tax.ls <- dups_tax.ls[! names(dups_tax.ls) %in% names(noDups_tax.ls)]


# convert to data frame
noDups_tax <- rbind(noDups_tax, data.frame(do.call(rbind, noDups_tax.ls)[-9]))


# one can't decide if it's acidovorax defluvii or carolinensis.
# blastn search said most likely defluvii
dups_tax.ls[[1]] <- dups_tax.ls[[1]][dups_tax.ls[[1]]$Species == "defluvii",]


# combine with uniques
dups_tax <- dups_tax.ls[[1]]
noDups_tax <- rbind(noDups_tax, dups_tax[-9])
unique(noDups_tax$R1 %in% fastas_uniq$R1)


# add ASV
taxa_uniq <- merge(noDups_tax, fastas_uniq, by = "R1")




############
### save ###
############

# match
taxa_uniq <- taxa_uniq[order(taxa_uniq$UNIQ),]
identical(colnames(counts_uniq), taxa_uniq$UNIQ)


# reorder
taxa_uniq <- taxa_uniq[c("UNIQ", "Kingdom", "Phylum", "Class",
                         "Order", "Family", "Genus", "Species", "R1")]


# add file names
counts_uniq <- data.frame(orginalFile = rownames(counts_uniq), counts_uniq)


# write
write.csv(counts_uniq, "./RData/01_counts.csv", row.names = F)
write.csv(taxa_uniq, "./RData/01_taxa.csv", row.names = F)


