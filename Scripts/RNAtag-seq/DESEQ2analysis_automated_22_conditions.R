library('DESeq2')
library('ggplot2')
library('grid')
library('gridExtra')
library('ggfortify')

setwd('~/OneDrive - University of Warwick/OneDrive/PhD/filtered_text_files/gff_ncRNAs_cappable_13_08_19/')
#read the count.data file with all the read counts
count.all = read.csv("./counts.all.csv")#remember to change this if you repeat, this if for no ncRNAs PCA plot
rownames(count.all) = count.all$X
count.all = count.all[,-1]
col.data = read.csv("./master_metadata.csv")
rownames(col.data) = col.data$Sample_name
#check the column names between matrix and counts.all
col.data$Sample_name[!col.data$Sample_name %in% colnames(count.all)]
#check which names are not matching
colnames(count.all)[!colnames(count.all) %in% col.data$Sample_name]

ddscond <- DESeqDataSetFromMatrix(countData = count.all, colData = col.data, design=~Condition)
colData(ddscond)
ddscond<-DESeq(ddscond)
rldTcond <- rlog(ddscond, blind=TRUE)
jennapalette <- c("#FFFFFF", "#FFA8BB", "#0075DC", "#993F00", "#191919", "#005C31", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#FF0010", "#5EF1F2", "#00998F", "#740AFF", "#FFFF00", "#990000","#F0A3FF", "#669999", "#FFA405", "#426600")
plotPCA(rldTcond, intgroup=c('Condition'))+scale_color_manual(values=jennapalette)


sampleDists <- dist(t(assay(rldTcond)))
sampleDistMatrix <- as.matrix(sampleDists)
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,cex = 0.85, angle_col = 315)

#pairwise comparisons
condition_names = unique(col.data[,7])
#con2 = unique(ddscond$Condition)
for (index1 in  condition_names){
  for(index2 in condition_names){
    if(index1 == index2){
      next
    }else{
      print(index2)
      print(index1)
      cond1 <- toString(index1)
      cond2 <- toString(index2)
      res <- results(ddscond, contrast=c("Condition",cond2,cond1),alpha =.05)
      mcols(res, use.names = TRUE)
      write.csv(as.data.frame(res),file=paste("results_pairwise_combination/",cond2,"_vs_",cond1,"control.csv",sep=""))

    }
  }
}

################################################################################
###### comp 1
################################################################################
#group comparisons
col.data2 = read.csv("./metadata_for_groups.csv")
rownames(col.data2) = col.data$Sample_name
#check the column names between matrix and counts.all
col.data2$Sample_name[!col.data2$Sample_name %in% colnames(count.all)]
#check which names are not matching
colnames(count.all)[!colnames(count.all) %in% col.data2$Sample_name]

ddsgroup <- DESeqDataSetFromMatrix(countData = count.all, colData = col.data2, design=~Temperature_growth_phase+Stress)



