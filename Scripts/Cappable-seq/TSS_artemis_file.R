############################################################################
######   Creating an artemis file for visualisation of TSS 

library("biofiles")
setwd("~/OneDrive - University of Warwick/OneDrive/PhD/Cappable-Seq/R_scripts_artemis_file/")

file = read.csv("categorised_TSS_correct_rank.csv")
Arnaud = read.delim("12864_2013_5327_MOESM23_ESM.art", header = F) #just to see what the artemis file looks like
TSS_list = file$V4
TSS_strand = file$V7
locus_tag = file$Locus_tag
TSS_type = file$TSS_type
TSS_identify = file$TSS_identifier

output_vec = rep(NA, (length(TSS_list)*4)) # <- starting point
for(i in 1:length(output_vec)){
  FT = "FT"
  output_vec[i] = FT #add feature at the beginning
  
}
colour_function <- function(TSStype){ #to determine the colour of the cell based on TSS type
  if(TSStype == "Primary"){
    return(12)
  }
  if(TSStype == "Secondary"){
    return(7)
  }
  if(TSStype == "Internal"){
    return(10)
  }
  if(TSStype == "Antisense"){
    return(8)
  }
  if(TSStype == "Orphan"){
    return(2)
  }
}
sigma_70 = as.character(file$sigma70)
sigma_70[which(is.na(sigma_70))] = "N"
sigma_54 = as.character(file$sigma54)
sigma_54[which(is.na(sigma_54))] = "N"
sigma_28 = as.character(file$sigma28)
sigma_28[which(is.na(sigma_28))] = "N"
sigma_factor <- function(sigma70, sigma54, sigma28){ #for writing a string for sigma factor
  if(sigma70 == "Y" & sigma54 == "Y"){
    return("sigma 70 and 54 promoter")
  }
  else if(sigma70 == "Y" & sigma28 == "Y"){
    return("sigma 70 and 28 promoter")
  }
  else if(sigma70 == "Y(weak)" & sigma54 == "Y"){
    return("sigma 70 (weak) and 54 promoter")
  }
  else if(sigma70 == "Y(weak)" & sigma28 == "Y"){
    return("sigma 70 (weak) and 28 promoter")
  }
  else if(sigma70 == "Y(very weak)" & sigma54 == "Y"){
    return("sigma 70 (very weak) and 54 promoter")
  }
  else if(sigma70 == "Y(very weak)" & sigma28 == "Y"){
    return("sigma 70 (very weak) and 28 promoter")
  }
  else if(sigma70 == "Y" & sigma54 =="N" & sigma28 == "N"){
    return("sigma 70 promoter")
  }
  else if(sigma70 == "Y(weak)" & sigma54 =="N" & sigma28 == "N"){
    return("sigma 70 (weak) promoter")
  }
  else if( sigma70 == "Y(very weak)" & sigma54 == "N" & sigma28 == "N"){
    return("sigma 70 (very weak) promoter")
  }
  else if(sigma54 == "Y" & sigma70 == "N"){
    return("sigma 54 promoter")
  }
  else if(sigma28 == "Y" & sigma70 == "N"){
    return("sigma 28 promoter")
  }
  else{
    return("no promoter motif")
  }
}

sharma_TSS = as.numeric(file$sharma_published)
sharma_TSS[which(is.na(sharma_TSS))] = 0
arnaud_TSS = as.numeric(file$arnaud_published)
arnaud_TSS[which(is.na(arnaud_TSS))] = 0
published <- function(sharma, arnaud){
  if(sharma > 0  & arnaud >0){
    return("found in Dugar et al, 2013 and Porcelli et al, 2013")
  }
  else if(sharma > 0 & arnaud == 0){
    return("found in Dugar et al, 2013")
  }
  else if(sharma == 0 & arnaud > 0){
    return("found in Porcelli et al, 2013")
  }
  else if(sharma == 0 & arnaud == 0){
    return("novel putative promoter")
  }
}

ToNER = as.character(file$ToNER_enriched)
ToNER[which(is.na(ToNER))] = "N"
toner_enriched <- function(toner){
  if(toner == "N"){
    
  }
  else{
    return("enriched with ToNER")
  }
}
for(j in 1:length(TSS_strand)){#no 4 times the size
  #print(j)
  #print((4*j)+1) #to access the first index of every 4 blocks
#but you need to make it stop otherwise it will do 5197-1 * 4 +1 which is like 20000 something but you do have that times 4 ...
  output_vec[(4*j)-3] = paste0(output_vec[(4*j)-3],"   ", "promoter","        ")#add 3 spaces, then promoter then 8 spaces
  output_vec[(4*j)-2] = paste0(output_vec[(4*j)-2], "                    ", '/locus_tag="', locus_tag[j], '"')#add 18 spaces, then \locus_tag= also " is really annoying so use ' and also / is annoying
  output_vec[(4*j)-1] = paste0(output_vec[(4*j)-1], "                    ", '/note="Promoter ', locus_tag[j],", ", TSS_type[j], " TSS, ", sigma_factor(sigma_70[j], sigma_54[j], sigma_28[j]), ", ", published(sharma_TSS[j], arnaud_TSS[j]), ", ", toner_enriched(ToNER[j]),", ",TSS_identify[j],'"' ) #the backslash isn't actually there
  output_vec[(4*j)] = paste0(output_vec[(4*j)], "                    ", "/colour=", colour_function(TSS_type[j]))
  if(TSS_strand[j] == "+"){
      #print("Yes")
     output_vec[(4*j)-3] = paste0(output_vec[(4*j)-3],as.character(TSS_list[j]-49), "..", as.character(as.numeric(TSS_list[j])))#minus 50 cos upstream
  }
  else{
     output_vec[(4*j)-3] = paste0(output_vec[(4*j)-3],"complement(",as.character(TSS_list[j]), "..", as.character(as.numeric(TSS_list[j])+49),")")
    }
}
output_vec
length(output_vec)
test = as.data.frame(output_vec)
rownames(test) <- c()
write.csv(test, file="TSS_artemis_file.csv")
write.FeatureTable(output_vec, file="test")


###### promoters only

output_vec = rep(NA, length(TSS_list)) # <- starting point
for(i in 1:length(output_vec)){
  FT = "FT"
  output_vec[i] = FT #add feature at the beginning
  
}
for(j in 1:length(TSS_strand)){#no 4 times the size
  #print(j)
  #print((4*j)+1) #to access the first index of every 4 blocks
  #but you need to make it stop otherwise it will do 5197-1 * 4 +1 which is like 20000 something but you do have that times 4 ...
  output_vec[j] = paste0(output_vec[j],"   ", "promoter","        ")#add 3 spaces, then promoter then 8 spaces
#  output_vec[(4*j)-2] = paste0(output_vec[(4*j)-2], "                    ", '/locus_tag="', locus_tag[j], '"')#add 18 spaces, then \locus_tag= also " is really annoying so use ' and also / is annoying
#  output_vec[(4*j)-1] = paste0(output_vec[(4*j)-1], "                    ", '/note="Promoter ', locus_tag[j],", ", TSS_type[j], " TSS, ", sigma_factor(sigma_70[j], sigma_54[j], sigma_28[j]), ", ", published(sharma_TSS[j], arnaud_TSS[j]), ", ", toner_enriched(ToNER[j]),'"' ) #the backslash isn't actually there
#  output_vec[(4*j)] = paste0(output_vec[(4*j)], "                    ", "/colour=", colour_function(TSS_type[j]))
  if(TSS_strand[j] == "+"){
    #print("Yes")
    output_vec[j] = paste0(output_vec[j],as.character(TSS_list[j]), "..", as.character(as.numeric(TSS_list[j])+50))
  }
  else{
    output_vec[j] = paste0(output_vec[j],"complement(",as.character(TSS_list[j]), "..", as.character(as.numeric(TSS_list[j])+50),")")
  }
}
output_vec
length(output_vec)

##########3 tying to figure out what the hell is going on
