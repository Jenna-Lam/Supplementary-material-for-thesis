####need a script which determines whether the TSS is primary, secondary, antisense, internal etc ... advice from Stephen
setwd("~/OneDrive - University of Warwick/OneDrive/PhD/Cappable-Seq/R_scripts_categorise_TSS/")
###find the matches from ToNER

TSS_file = read.csv("TSS_sigma_factors.csv")
TSS_pos = TSS_file$V4
TSS_strand = TSS_file$V7
TSS_enriched = TSS_file$ToNER_enriched
ToNER_file=read.csv("ToNER_sorted.csv")
ToNER_TSS = ToNER_file$position
ToNER_strand=ToNER_file$strand
ToNER_region = ToNER_file$region
ToNER_feat = ToNER_file$feature
ToNER_dist = ToNER_file$distance

split_TSS_enriched <- function(input){#this function splits the numbers from ToNER enriched column ALSO MADE REDUNDANT COS WE USE GFF FILE
  input = as.vector(input)
  output = as.numeric(unlist(strsplit(input, ";")))
  return(output)
}
threshold=200
TSS_type = rep(NA, length(TSS_pos))
for(i in 1: length(TSS_pos)){
  print(paste0(i,"/",length(TSS_pos)))
  for(j in 1:length(ToNER_TSS)){
    #if TSS from cappable-seq is also enriched in ToNER
    if(is.na(TSS_enriched[i]) == "FALSE" &ToNER_TSS[j] %in% split_TSS_enriched(TSS_enriched[i]) &TSS_strand[i] == ToNER_strand[j]){
      if(as.character(ToNER_feat[j]) == "gene,CDS"){
        if(as.character(ToNER_region[j]) == "Upstream"){
          if(ToNER_dist[j] < threshold){
            TSS_type[i] = "Primary"
            print("primary")
          }
          else if(ToNER_dist[j] >= threshold & is.na(TSS_type[i])){
            TSS_type[i] = "Secondary"
            print("secondary")

          }
        }
        if(substr(as.character(ToNER_region[j]), 1,5) == "Intra"){
          TSS_type[i] = "Internal"
          print("internal")
        }
      }
      if(as.character(ToNER_feat[j]) == "ncRNA"){
    }
  }
}

######################################################  <- START HERE
######do this with gff file

TSS_file = read.csv("TSS_sigma_factors.csv")
TSS_pos = TSS_file$V4
TSS_strand = TSS_file$V7
TSS_rrs = TSS_file$V6
GFF_file = read.delim("11168_genome_genes_only_jenna_published.gff3.txt")#modified to have gene features, pseudogenes, and rRNA only don't include tRNAs ...
GFF_file = GFF_file[-1,]
GFF_start = GFF_file$X.2
GFF_feature = GFF_file$X.1
GFF_strand = GFF_file$X.5
GFF_end = GFF_file$X.3
GFF_locus_tag = GFF_file$X.7


#retrieve locus tag, different rules for each, twas a pain
new_locus_tag = rep(NA,length(GFF_locus_tag))
for(i in 1:length(GFF_locus_tag)){
  ID = as.character(GFF_locus_tag[i])#make it a string
  ID_split = strsplit(ID, "locus_tag=")[[1]]#need double bracket 1 to make first element become a vector
  if(GFF_feature[i] == "pseudogene"){
    tag = strsplit(ID_split[length(ID_split)], ";")[[1]]
    pseudo_tag = tag[1]
    part = strsplit(tag[2], "=")[[1]]
    if(part[1] == "part"){
      part_tag = paste0("p", part[2])
    }
    new_locus_tag[i] = paste0(pseudo_tag, part_tag)
    print(new_locus_tag[i])
  }
  if(GFF_feature[i] == "rRNA"){
    tag = strsplit(ID_split[length(ID_split)], "=")[[1]]
    rna_tag = strsplit(tag[length(tag)], ' ')[[1]] #the last element is the 16S ribosomal RNA and split by space
    rna_tag_beg = rna_tag[1]# get the first element so 16S
    rna_num = strsplit(tag[2],";")[[1]] #split the second element from before so the /name=rRNA1 fpr e.g.
    rna_num_1 = rna_num[1]
    rna = paste0(rna_tag_beg,"_", rna_num_1)
    new_locus_tag[i] = rna
  }
  if(GFF_feature[i] =="gene"){
    tag = ID_split[length(ID_split)] #length of the vector to get the last entry
    new_tag = (strsplit(tag,";")[[1]])[1]
    if(startsWith(new_tag, "I")){#some tRNAs are classed as genes
      alt_tag = (strsplit(tag, ";")[[1]])[2]#Name=tRNA_Gly
      name = strsplit(alt_tag, "=")[[1]]#Name tRNA_Gly
      gene_number = strsplit(new_tag, "=")[[1]]#ID gene844
      together = paste0(gene_number[2], "_", name[2])
      new_locus_tag[i] = together
    }
    else{
      new_locus_tag[i] = new_tag
    }

  }

}

length(new_locus_tag)#should be same as GFF_locus_tag which is 1723

threshold=500 #300 is what sharma used, use 500 which is what arnoud used and redefine secondary promoters
TSS_type = rep(NA, length(TSS_pos))
locus_tag = rep(NA, length(TSS_pos))

primary_TSS <- function(TSSpos, GFFstart, GFFend, GFFstrand,num){#function to determine the primary TSS
  if(GFFstrand == "+"){
    if((GFFstart - TSSpos)<=num & ((GFFstart - TSSpos) >0)){#start of gene minus the TSS is less than 300bp but also over 0bp otherwise will count negative numbers
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  }
  else{
    if((TSSpos - GFFend)<= num & ((TSSpos - GFFend)>0)){#start of TSS minus the gene end since it is in the forward direction
      return(TRUE)
    }
    else{#need this cos of stupid argument is of length zero
      return(FALSE)
    }
  }


}

internal_TSS <- function(TSSpos, GFFstart, GFFend){
  if(TSSpos %in% (GFFstart:GFFend)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

for(j in 1:length(GFF_start)){#for every position in GFF file
  print(paste0(j,"/",length(GFF_start)))
  num = 0
  for(i in 1:length(TSS_pos)){#loop through the TSS positions
    
    if(GFF_strand[j] == TSS_strand[i]){#if the strand position is the same
      if(primary_TSS(TSS_pos[i], GFF_start[j], GFF_end[j], GFF_strand[j], threshold)){#first because otherwise will overwrite internal
        if(is.na(TSS_type[i])){
          TSS_type[i] = "Primary"
          print("primary")
          locus_tag[i] = new_locus_tag[j]
        }
        
        
      }
      
      else if(internal_TSS(TSS_pos[i], GFF_start[j], GFF_end[j]) == TRUE){
        num = num+1
        if(num > 1){
          TSS_type[i] = "Internal"
          print("internal")
          locus_tag[i] = paste0(new_locus_tag[j], "_int_", num)
        }
        else{
          TSS_type[i] = "Internal"
          print("internal")
          locus_tag[i] = paste0(new_locus_tag[j], "_int")
        }
        
        
      }



    }
     else{
       next
     } 
  }
}




TSS_file$TSS_type = TSS_type
TSS_file$Locus_tag = locus_tag


secondary_TSS<- function(firstTSS, secondTSS, strand, firstTSStype, secondTSStype, GFFstart, GFFend){
  if(strand == "+"){
    if((GFFstart - firstTSS) <= 500 & (GFFstart - firstTSS) > 0 & is.na(firstTSStype) & secondTSS %in% ((GFFstart-300):GFFstart)){#within 500bp of gene start, entry is NA, and primary downstream of secondary TSS(i.e. the first one if positive strand)
      return(TRUE)
      #& secondTSS %in% ((GFFstart-300):GFFstart) insert back in if it doesn't work, yes need it because if value = NA it doesn't work
    }
    else{
      return(FALSE)
    }
  }#always need a FALSE for every TRUE/FALSE statement in all functions -_-
  else if((secondTSS - GFFend) <= 500 & (secondTSS - GFFend) > 0 & is.na(secondTSStype) & firstTSS %in% (GFFend:GFFend+300)){#otherwise if negative
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  
}


antisense_TSS <- function(TSSpos, GFFstart,GFFend,TSSstrand){
  if(TSSstrand == "-"){
    if(TSSpos %in% (GFFstart:GFFend) | ((TSSpos - GFFend) <= 100 & (TSSpos - GFFend)>0)){#if opposite strand (in for loop) and within gene, or within 100nt of start or end of gene feature
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  }
  else{
    if(TSSpos %in% (GFFstart:GFFend) |((GFFstart-TSSpos) <= 100 & (GFFstart-TSSpos) >0)){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  }

}
#########loop after for secondary, and antisense and orphan

################################# DON'T RUN THIS, MADE REDUNDANT NOW  
TSS_type = TSS_file$TSS_type
locus_tag = TSS_file$Locus_tag
#newTSS_type = rep(NA, length(TSS_type))
for(a in 1:(length(TSS_type))){#two loops so can compare the TSS in pairs, first TSS
    #print(a)
#second TSS
    for(c in 1:length(GFF_start)){
      #print(paste0("c=", c, "_a=", a))
      if(TSS_strand[a] == GFF_strand[c]){#if matches with GFF strand
          #print("yes")
        if(TSS_strand[a] == TSS_strand[a+1] & a < (length(TSS_type))){#if first strand matches with second TSS strand, and doesn't count the last one
            if(secondary_TSS(TSS_pos[a], TSS_pos[a+1], TSS_strand[a], TSS_type[a], TSS_type[a+1], GFF_start[c], GFF_end[c])){ #firstTSS, secondTSS, strand, firstTSStype, secondTSStype, GFFstart
            
            TSS_type[a] = "Secondary"
            print(paste0("secondary_c_", c, "_a_", a))
            locus_tag[a] = new_locus_tag[c]
          }
          else{
            next
              }
            }

             }
      else{#need an else for everything ...
        next
      }
      
  }

}

###################################################################################################### -< start here
############antisense
TSS_type = TSS_file$TSS_type
locus_tag = TSS_file$Locus_tag
for(j in 1:length(GFF_start)){#for every TSS
  num=0
  for(i in 1:length(TSS_type)){#for every GFF entry
    
    if(GFF_strand[j] != TSS_strand[i] & is.na(TSS_type[i])){#if strand info DOES NOT match for antisense and also if not catergorised already
      if(antisense_TSS(TSS_pos[i], GFF_start[j], GFF_end[j], TSS_strand[i])){#TSSpos, GFFstart,GFFend,TSSstrand
        num = num+1
        if(num >1){
          TSS_type[i] = "Antisense"
          print(paste0("antisense_i=",i,"_j=",j))
          locus_tag[i] = paste0("as_",new_locus_tag[j], "_", num)
        }
        else{
          TSS_type[i] = "Antisense"
          print(paste0("antisense_i=",i,"_j=",j))
          locus_tag[i] = paste0("as_",new_locus_tag[j])
        }

      }
    }
  }
}

####orphan
num =0
for(i in 1:length(TSS_type)){
  if(is.na(TSS_type[i])){
    num=num+1
    TSS_type[i] = "Orphan"
    locus_tag[i] = paste0("intergenic_", num)
  }
}

TSS_file$TSS_type = TSS_type
TSS_file$Locus_tag = locus_tag





######## alternative secondary promoters


##### if doesn't work then clear everything
write.csv(TSS_file, file="for_Stephen.csv")

alt_TSS = TSS_file
alt_TSS = read.csv("for_Stephen.csv")
primary = alt_TSS[which(alt_TSS$TSS_type == "Primary"),]
primary = primary[order(primary$Locus_tag),]#order by locus tag so names are together
non_primary = alt_TSS[which(alt_TSS$TSS_type != "Primary"),]
old_locus = primary$Locus_tag
new_type = c()
new_locus_vec = c()
TSS_rrs = primary$V9
length(TSS_rrs)
rrs = c()
for(a in 1:length(TSS_rrs)){
  new_rrs = strsplit(as.character(TSS_rrs[a]), "_")[[1]][4]
  rrs = c(rrs, new_rrs)
}
length(rrs)

i=1

while(i <= length(old_locus)){
#while(i < 10){
  index = c(i)
  j=1
  first = as.character(old_locus[i])
  print(paste0(i,"_is_i"))
  while(j >= 1){
    second = as.character(old_locus[j+i]) #i+j moves along each successive index
    if (is.na(second)){
      
      break
    }
    #print(c(first,second))
    if(first != second){ #no need for strand info as assume that positive and negative locus tag names are the same, if opposite strand won't have c
      #new_locus_vec = c(new_locus_vec, as.character(old_locus[i]))
      
      print ("not the same")
     
      
      break
    }
    else{
      index = c(index, j+i)#row number
      j = j+1 #move onto the next position
      print("same")
      #print(j)
    }
  }
  #
  rss_score=rrs[index]
  print(paste0("rss_score_is_", rss_score))
  for(number in rank(-(as.numeric(rss_score)))){#rank the vector with all the indexes
    
    print(paste0("the numbers are: ",number))
    print(paste0("the rank are: ",rank(-(as.numeric(rss_score)))))# needs to be negative so the highest number is ranked 1st
    if(number == 1 & length(rank(-(as.numeric(rss_score))))>1){
      print(paste0("this is ranked 1, i =", i, " j =", j ))
      new_type = c(new_type, "Primary")
      lt= paste0(as.character(old_locus[i]), "_1")
      new_locus_vec = c(new_locus_vec, lt)
    }
    else if (number == 1 & length(rank(-(as.numeric(rss_score))))==1){
      new_type = c(new_type, "Primary")
      lt= as.character(old_locus[i])
      new_locus_vec = c(new_locus_vec, lt)
    }
    else if (number>1){ #because some cases have the same rrs score so it is name 1.5 instead. Need to somehow sort those by position instead
      print(paste0("this is alternative, i =", i, " j =", j ))
      new_type = c(new_type, "Secondary")
      lt= paste0(as.character(old_locus[i]), "_", number)
      new_locus_vec = c(new_locus_vec, lt)
      print(lt)
    }
    

  }
  print(paste0((i+j),"_i_and_j):(length_of_new_locus_vec)",length(new_locus_vec)))
  if ((i+j)-(length(new_locus_vec)) >1){
    print(i+j)
    break
  }
  i=i+j #need this so it skips the window that i was comparing to i.e. j
  j=1
  #if(i == length(old_locus )){
   # new_type=c(new_type,"Primary")
   # new_locus_vec = c(new_locus_vec, old_locus[i])
   # break
  #}
}


length(new_locus_vec)
length(new_type)

primary$TSS_type = new_type
primary$Locus_tag = new_locus_vec

##### naming the rest
##### manually change the 1.5, change the script later
test = rbind(primary, non_primary)
test = test[order(test$V4),]

write.csv(test, file="categorised_TSS_correct_rank2.csv")
