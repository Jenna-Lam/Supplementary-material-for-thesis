library('plyr')
library('dplyr')

setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/RNAtagSeq (differential gene expression)/Total_DESeq2_output_all_comparisons_cappable/')

files <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/RNAtagSeq (differential gene expression)/Total_DESeq2_output_all_comparisons_cappable/')
files = files[grep(".csv", files)] 
gene_names = (read.csv("stress_vs_37_exponentialcontrol.csv", sep=",", header=F))$V1
gene_names = gsub('"','',gene_names)


master_counts <- data.frame(matrix(" ", nrow=3381, ncol =1))
master_counts[1] <- gene_names
colnames(master_counts) = "X"
for(index in 1:length(files)) {
  #fileName = paste0(myfolder,file)
  counts = read.delim(files[index], sep=",")#separate by comma
  #condition = gsub('"','',counts$V1)
  counts = select(counts, X,log2FoldChange,pvalue,padj)
  
  counts$X = gsub('"','',counts$X)
  master_counts = join(master_counts, counts, by ="X")
  master_counts[1,(3*index-1)] = files[index] 
}
#rownames(master_counts) = master_counts$V1  #use the locus tag
names(master_counts)[1] = "gene_names"


#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts, file="master_counts.txt", sep = "\t", row.names=FALSE)

###################################################################################################################################################
###########  Significant genes only
###################################################################################################################################################
master_Sigcounts <- data.frame(matrix(" ", nrow=3381, ncol =1))
master_Sigcounts[1] <- gene_names
colnames(master_Sigcounts) = "X"
for(index in 1:length(files)) {
  counts = read.csv(files[index], sep=",")#separate by comma
  #condition = gsub('"','',counts$V1)
  countsSig <- counts[which(counts$padj <=0.05),]
  countsSig = select(countsSig, X,log2FoldChange)
  #names(counts)[2] = file
  countsSig$X = gsub('"','',countsSig$X)
  
  master_Sigcounts = join(master_Sigcounts, countsSig, by ="X")
  #master_Sigcounts[1,(3*index-1)] = files[index] 
  master_Sigcounts[1,(index+1)] = files[index]
}
#rownames(master_counts) = master_counts$V1  #use the locus tag
names(master_Sigcounts)[1] = "gene_names"
master_Sigcounts[is.na(master_Sigcounts)] <- 0


#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_Sigcounts, file="master_Sigcounts.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

#############################################################################################################################
#everything above log2fold change of >1 and <1 and p <= 0.05

master_Sig1counts <- data.frame(matrix(" ", nrow=3381, ncol =1))
master_Sig1counts[1] <- gene_names
colnames(master_Sig1counts) = "X"
for(index in 1:length(files)) {
  counts = read.csv(files[index], sep=",")#separate by comma
  #condition = gsub('"','',counts$V1)
  counts_p <- counts[which(counts$padj <=0.05),]
  counts_over_1 <- counts_p[which(counts_p$log2FoldChange >=1),]
  counts_below_1 <- counts_p[which(counts_p$log2FoldChange <= -1),]
  countsSig = rbind(counts_over_1, counts_below_1)
  countsSig = select(countsSig, X,log2FoldChange)
  #names(counts)[2] = file
  countsSig$X = gsub('"','',countsSig$X)
  
  master_Sig1counts = join(master_Sig1counts, countsSig, by ="X")
  #master_Sigcounts[1,(3*index-1)] = files[index] 
  master_Sig1counts[1,(index+1)] = files[index]
}
#rownames(master_counts) = master_counts$V1  #use the locus tag
names(master_Sig1counts)[1] = "gene_names"
master_Sig1counts[is.na(master_Sig1counts)] <- 0


#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_Sig1counts, file="master_Sigcutoff_1_counts.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

binary = read.delim("master_Sigcutoff_1_counts.txt", header = T)
for(i in 1:nrow(binary)){
  for(j in 2:ncol(binary)){
    if(binary[i,j] > 0){
      binary[i,j] = 1
    }
    if(binary[i,j] < 0){
      binary[i,j] = -1
    }
  }
}

write.table(binary, file="differentially_expressed_comparisons_binary_cappable_cutoff1.txt", sep = "\t")
#############################################################################################################################
#######Growth phase temperature effects
#############################################################################################################################
setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_temperature_effects/')
files_t <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_temperature_effects/')
files_t = files_t[grep(".csv", files_t)] 
gene_names_t = (read.csv("37_early_stationary_vs_37_exponentialcontrol.csv", sep=",", header=F))$V1
gene_names_t = gsub('"','',gene_names_t)

master_counts_t <- data.frame(matrix(" ", nrow=3122, ncol =1))
master_counts_t[1] <- gene_names_t
colnames(master_counts_t) = "X"
for(index in 1:length(files_t)) {
  #fileName = paste0(myfolder,file)
  counts = read.csv(files_t[index], sep=",")#separate by comma
  countsSig <- counts[which(counts$padj <0.05),]
  #condition = gsub('"','',counts$V1)
  countsSig = select(countsSig, X, log2FoldChange)
  
  countsSig$X = gsub('"','',countsSig$X)
  master_counts_t = join(master_counts_t, countsSig, by ="X")
  master_counts_t[1,(index+1)] = files_t[index]
  #names(counts)[2] = file
}
master_counts_t[is.na(master_counts_t)] <- 0

#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts_t, file="temperature_effectsSig.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

#############################################################################################################################
####### bile effects
#############################################################################################################################

setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_bile_effects/')
files_b <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_bile_effects/')
files_b = files_b[grep(".csv", files_b)] 
gene_names_b = (read.csv("sodium_deoxycholate_exponential_vs_37_exponentialcontrol.csv", sep=",", header=F))$V1
gene_names_b = gsub('"','',gene_names_b)

master_counts_b <- data.frame(matrix(" ", nrow=3122, ncol =1))#change number of rows accordingly
master_counts_b[1] <- gene_names_b
colnames(master_counts_b) = "X"
for(index in 1:length(files_b)) {
  #fileName = paste0(myfolder,file)
  counts = read.csv(files_b[index], sep=",")#separate by comma
  countsSig <- counts[which(counts$padj <0.05),]
  #condition = gsub('"','',counts$V1)
  countsSig = select(countsSig, X, log2FoldChange)
  
  countsSig$X = gsub('"','',countsSig$X)
  master_counts_b = join(master_counts_b, countsSig, by ="X")
  master_counts_b[1,(index+1)] = files_b[index]
  #names(counts)[2] = file
}
master_counts_b[is.na(master_counts_b)] <- 0

#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts_b, file="bile_effectsSig.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

#############################################################################################################################

setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_with_ncRNAs/results_bile_effects_version2/')
files_b <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_with_ncRNAs/results_bile_effects_version2/')
files_b = files_b[grep(".csv", files_b)] 
gene_names_b = (read.csv("sodium_deoxycholate_exponential_vs_37_exponentialcontrol.csv", sep=",", header=F))$V1
gene_names_b = gsub('"','',gene_names_b)

master_counts_b <- data.frame(matrix(" ", nrow=2087, ncol =1))
master_counts_b[1] <- gene_names_b
colnames(master_counts_b) = "X"
for(index in 1:length(files_b)) {
  #fileName = paste0(myfolder,file)
  counts = read.csv(files_b[index], sep=",")#separate by comma
  countsSig <- counts[which(counts$padj <0.05),]
  #condition = gsub('"','',counts$V1)
  countsSig = select(countsSig, X, log2FoldChange)
  
  countsSig$X = gsub('"','',countsSig$X)
  master_counts_b = join(master_counts_b, countsSig, by ="X")
  master_counts_b[1,(index+1)] = files_b[index]
  #names(counts)[2] = file
}
master_counts_b[is.na(master_counts_b)] <- 0

#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts_b, file="bile_effectsSig_v2.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

#############################################################################################################################
####### iron effects
#############################################################################################################################

setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_iron_effects/')
files_iron <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_iron_effects/')
files_iron = files_iron[grep(".csv", files_iron)] 
gene_names_iron = (read.csv("37_degrees_vs_iron_limitationcontrol.csv", sep=",", header=F))$V1
gene_names_iron = gsub('"','',gene_names_iron)

master_counts_iron <- data.frame(matrix(" ", nrow=3122, ncol =1))
master_counts_iron[1] <- gene_names_iron
colnames(master_counts_iron) = "X"
for(index in 1:length(files_iron)) {
  #fileName = paste0(myfolder,file)
  counts = read.csv(files_iron[index], sep=",")#separate by comma
  countsSig <- counts[which(counts$padj <0.05),]
  #condition = gsub('"','',counts$V1)
  countsSig = select(countsSig, X, log2FoldChange)
  
  countsSig$X = gsub('"','',countsSig$X)
  master_counts_iron = join(master_counts_iron, countsSig, by ="X")
  master_counts_iron[1,(index+1)] = files_iron[index]
  #names(counts)[2] = file
}
master_counts_iron[is.na(master_counts_iron)] <- 0

#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts_iron, file="iron_effectsSig.txt", sep = "\t", row.names=FALSE, col.names = FALSE)


#############################################################################################################################
####### chicken juice and cold effects
#############################################################################################################################

setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_cold_effects/')
files_cold <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_cold_effects/')
files_cold = files_cold[grep(".csv", files_cold)] 
gene_names_cold = (read.csv("cold_vs_37_exponentialcontrol.csv", sep=",", header=F))$V1
gene_names_cold = gsub('"','',gene_names_cold)

master_counts_cold <- data.frame(matrix(" ", nrow=3122, ncol =1))
master_counts_cold[1] <- gene_names_cold
colnames(master_counts_cold) = "X"
for(index in 1:length(files_cold)) {
  #fileName = paste0(myfolder,file)
  counts = read.csv(files_cold[index], sep=",")#separate by comma
  countsSig <- counts[which(counts$padj <0.05),]
  #condition = gsub('"','',counts$V1)
  countsSig = select(countsSig, X, log2FoldChange)
  
  countsSig$X = gsub('"','',countsSig$X)
  master_counts_cold = join(master_counts_cold, countsSig, by ="X")
  master_counts_cold[1,(index+1)] = files_cold[index]
  #names(counts)[2] = file
}
master_counts_cold[is.na(master_counts_cold)] <- 0

#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts_cold, file="cold_effectsSig.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

#############################################################################################################################
####### stresses
#############################################################################################################################

setwd('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_stress_effects/')
files_s <- list.files('/Users/u1560837/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_23_06_19/results_stress_effects/')
files_s = files_s[grep(".csv", files_s)] 
gene_names_s = (read.csv("acid_vs_37_exponentialcontrol.csv", sep=",", header=F))$V1
gene_names_s = gsub('"','',gene_names_s)

master_counts_s <- data.frame(matrix(" ", nrow=3122, ncol =1))
master_counts_s[1] <- gene_names_s
colnames(master_counts_s) = "X"
for(index in 1:length(files_s)) {
  #fileName = paste0(myfolder,file)
  counts = read.csv(files_s[index], sep=",")#separate by comma
  countsSig <- counts[which(counts$padj <0.05),]
  #condition = gsub('"','',counts$V1)
  countsSig = select(countsSig, X, log2FoldChange)
  
  countsSig$X = gsub('"','',countsSig$X)
  master_counts_s = join(master_counts_s, countsSig, by ="X")
  master_counts_s[1,(index+1)] = files_s[index]
  #names(counts)[2] = file
}
master_counts_s[is.na(master_counts_s)] <- 0

#colnames(counts.all) = gsub("_"only_mapped_filt.txt","",files)
write.table(master_counts_s, file="stress_effectsSig.txt", sep = "\t", row.names=FALSE, col.names = FALSE)

