setwd("~/OneDrive - University of Warwick/OneDrive/PhD/Cappable-Seq/TSS_TPM_conditions/") #working directory
TPM_og = read.csv("TPM_TSS.csv") #List of TPM values for TSS transcripts
rownames(TPM_og) = TPM_og$X
TPM_og = TPM_og[,-1]
TPM_og = TPM_og %>% select(-cj100_1, -cj100_2, -cj100_3)

TSS_iron = TPM_og %>% select(ILE_1, ILE_2, ILE_3, ILS_1, ILS_2, ILS_3, IRS_1, IRS_2, IRS_3, IRE_1, IRE_2, IRE_3)

TSS_iron_filter = TSS_iron[which(((apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
  unlist(strsplit(x,"_"))[1])) == "ILE"], 1, function(y) sum(y>10)) >1) &
    (apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
      unlist(strsplit(x,"_"))[1])) == "ILS"], 1, function(y) sum(y>10)) >1))&
    ((apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
      unlist(strsplit(x,"_"))[1])) == "IRE"], 1, function(y) sum(y<=10)) >1) &
       (apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
         unlist(strsplit(x,"_"))[1])) == "IRS"], 1, function(y) sum(y<=10)) >1))),]

write.csv(TSS_iron_filter, file="TSS_iron_all.csv")

TSS_iron_sep = TSS_iron[which(((apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x) #if up in ILE and down in IRE
  unlist(strsplit(x,"_"))[1])) == "ILE"], 1, function(y) sum(y>10)) >1) & (apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
    unlist(strsplit(x,"_"))[1])) == "IRE"], 1, function(y) sum(y<=10)) >1)| #or
    ((apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x) #if up in ILS and down in IRS
      unlist(strsplit(x,"_"))[1])) == "ILS"], 1, function(y) sum(y>10)) >1)) &
       (apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
         unlist(strsplit(x,"_"))[1])) == "IRS"], 1, function(y) sum(y<=10)) >1))),]

TSS_iron_og = TPM_og[which(((apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
  unlist(strsplit(x,"_"))[1])) == "ILE"], 1, function(y) sum(y>10)) >1) &
    (apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
      unlist(strsplit(x,"_"))[1])) == "ILS"], 1, function(y) sum(y>10)) >1))&
    ((apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
      unlist(strsplit(x,"_"))[1])) == "IRE"], 1, function(y) sum(y<=10)) >1) &
       (apply(TSS_iron[,unlist(lapply(colnames(TSS_iron),function(x)
         unlist(strsplit(x,"_"))[1])) == "IRS"], 1, function(y) sum(y<=10)) >1))),]

write.csv(TSS_iron_og, file="TSS_iron_37M.csv")
