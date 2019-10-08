##########script for predicting promoters 
##########sigma28 = CGATwt (w means weak A or T) 6-8 nt upstream
##########sigma 54 = GGaa-N6-tTGCTt 8-13 nt upstream
##########sigma 70 = gnTAnaAT 4-8 nt upstream
library("Biostrings") #package reads the fastasequence
setwd("~/OneDrive - University of Warwick/OneDrive/PhD/Cappable-Seq/R_scripts_promoter_prediction/")
fasta_file = readDNAStringSet("11168_genome_updated.fasta")#function to read the fasta sequence
fasta_sequence_forward = fasta_file$`AL111168.1 `
fasta_sequence_reverse = reverseComplement(fasta_sequence_forward)
fasta_complement = complement(fasta_sequence_forward)
TSS_input = read.delim("TSS_enriched_cluster_10.gtf", header = F)
TSS_positions = TSS_input[,4]#get the 4th column for the TSS
TSS_strand = TSS_input[,7]#7th column for the strand

#sigma28 = "CGATwt"
#sigma54 = "GGaa-N6-tTGCTt"
#sigma70 = "gnTAnaAT"
promoter<- function(string, fasta_sequence){#function to find the promoter motif in the genome
  #vector_of_positions = c()
  #output_vector = rep(NA, length(TSS_positions))
  sigma_motif = DNAString(string)
  result = matchPattern(sigma_motif, fasta_sequence, fixed=FALSE)#finds a match but also fixed = false makes it more lenient so it can account for w
  return(result@ranges@start)
}
promoter70<- function(string, fasta_sequence){#function to find the promoter motif in the genome
  #vector_of_positions = c()
  #output_vector = rep(NA, length(TSS_positions))
  sigma_motif = DNAString(string)
  result = matchPattern(sigma_motif, fasta_sequence, fixed=FALSE, max.mismatch = 1)#finds a match but also fixed = false makes it more lenient so it can account for w
  return(result@ranges@start)
}
promoter_reverse<- function(string, fasta_sequence){#function to find the promoter motif in the genome

  sigma_motif = reverseComplement(DNAString(string))#essential to find the antisense promoters
  result = matchPattern(sigma_motif, fasta_sequence, fixed=FALSE)#finds a match but also fixed = false makes it more lenient so it can account for w
  return(result@ranges@start)#+(result@ranges@width -1)) #add this to get the end but then change the range if this is the case
}
promoter_reverse70<- function(string, fasta_sequence){#function to find the promoter motif in the genome
  
  sigma_motif = reverseComplement(DNAString(string))#essential to find the antisense promoters
  result = matchPattern(sigma_motif, fasta_sequence, fixed=FALSE, max.mismatch = 1)#finds a match but also fixed = false makes it more lenient so it can account for w
  return(result@ranges@start)#+(result@ranges@width -1)) #add this to get the end but then change the range if this is the case
}

veryweakF<- function(string, fasta_sequence){#function to find the promoter motif in the genome
  #vector_of_positions = c()
  #output_vector = rep(NA, length(TSS_positions))
  sigma_motif = DNAString(string)
  result = matchPattern(sigma_motif, fasta_sequence, fixed=FALSE, max.mismatch = 2)#finds a match but also fixed = false makes it more lenient so it can account for w
  return(result@ranges@start)
}

veryweakR<- function(string, fasta_sequence){#function to find the promoter motif in the genome
  
  sigma_motif = reverseComplement(DNAString(string))#essential to find the antisense promoters
  result = matchPattern(sigma_motif, fasta_sequence, fixed=FALSE, max.mismatch = 2)#finds a match but also fixed = false makes it more lenient so it can account for w
  return(result@ranges@start)#+(result@ranges@width -1)) #add this to get the end but then change the range if this is the case
}
#search for promoter motifs in the whole
sigma28F = promoter("CGATwt", fasta_sequence_forward) #the motifs from porcelli paper
sigma28R = promoter_reverse("CGATwt", fasta_sequence_forward) 
sigma54F = promoter("GGnnNNNNNNnTGCT", fasta_sequence_forward)
sigma54R = promoter_reverse("GGnnNNNNNNnTGCT", fasta_sequence_forward)
sigma70F = promoter("gnTAnaaT", fasta_sequence_forward)
sigma70R = promoter_reverse("gnTAnaaT", fasta_sequence_forward)
extra70F = promoter70("gnTAnaaT", fasta_sequence_forward)
extra70R = promoter_reverse70("gnTAnaaT", fasta_sequence_forward)
vw70F = veryweakF("gnTAnaaT", fasta_sequence_forward)
vw70R = veryweakR("gnTAnaaT", fasta_sequence_forward)

output28 = rep(NA, length(TSS_positions))#empty vectors 
output54 = rep(NA, length(TSS_positions))
output70 = rep(NA, length(TSS_positions))

for(i in 1:length(TSS_positions)){#do we need strand info ... yes we do

 # if(TSS_positions[i] <6){
 #   next
  #}
  if(TSS_strand[i] == "+"){
    range28F = ((TSS_positions[i]-8):(TSS_positions[i]-6))#6-8 nt upstream reaches the end of the motif so add the length of the sequence too so add 6
    range54F = ((TSS_positions[i]-13):(TSS_positions[i]-8))#add 15 because account for the TSS at the beginning so 8+15 -1
    range70F = ((TSS_positions[i]-8): (TSS_positions[i]-4))#add 8
    for(pos28F in 1:length(sigma28F)){
      if(sigma28F[pos28F+6] %in% range28F){
        #print("works28F")
        output28[i] = "Y"
        break
      }
    
    }
   for(pos54F in 1:length(sigma54F+14)){
      if(sigma54F[pos54F+15] %in% range54F){
        #print("works54F")
        output54[i] = "Y"
        break
      }   
    }
    for(pos70F in 1:length(sigma70F+7)){
      
      if(sigma70F[pos70F+8] %in% range70F){
       # print("works70F")
        output70[i] = "Y"
        break
      }
      
    }
      for(eF in 1:length(extra70F)){
        if(extra70F[eF] %in% range70F){
          if(is.na(output70[i])){
            output70[i] = "Y(weak)"
            break
          }
          
        }
      }
    if(is.na(output70[i])){
      for(vF in 1:length(vw70F)){
        if(vw70F[vF] %in% range70F){
          output70[i] = "Y(very weak)"
          break
        }
      }
    }

    
  }
  if(TSS_strand[i] == "-"){
    range28R = ((TSS_positions[i]+6):(TSS_positions[i]+8))#this is because when it is on the negative strand the start of the motif is the closest to the TSS, TSS is already 1 but need to catch the end of the motif ...
   
    range54R = ((TSS_positions[i]+8):(TSS_positions[i]+13))
    range70R = ((TSS_positions[i]+4): (TSS_positions[i]+8))
    for(pos28R in 1:length(sigma28R)){
     # print(range28R)
     # print(sigma28R[pos28R])
      #print(" ")
      if(sigma28R[pos28R] %in% range28R){
        
        #print("works28R")
        output28[i] = "Y"
        break
      }
    }
    for(pos54R in 1:length(sigma54R)){
      if(sigma54R[pos54R] %in% range54R){
        #print("works54R")
        output54[i] = "Y"
        break
      }   
    }
    for(pos70R in 1:length(sigma70R)){
      
      if(sigma70R[pos70R] %in% range70R){
        #print("works70R")
        output70[i] = "Y"
        break
      }

    }
    for(eR in 1:length(extra70R)){
      if(extra70R[eR] %in% range70R){
        if(is.na(output70[i])){
          output70[i] = "Y(weak)"
          break
        }
      }
    }
    if(is.na(output70[i])){
      for(vR in 1:length(vw70R)){
        if(vw70R[vR] %in% range70R){
          output70[i] = "Y(very weak)"
          break
        }
      }
    }
  }
}

length(which(output28 == "Y"))
length(which(output54 == "Y"))
length(which(output70 == "Y"))
length(which(output70 == "Y(weak)"))
length(which(output70 == "Y(very weak)"))

cap_TSS = read.csv("novel_TSS_sig_enriched.csv")
cap_TSS$sigma70 = output70
cap_TSS$sigma54 = output54
cap_TSS$sigma28 = output28



#function 1 . Find position of sigma 28, return a vector of position. Use readDNAStringSet to generate the fasta sequence, then use a for loop to find position of GCAT
#function 2 and 3. Do the same thing for sigma 54, sigma 70
#function 4: compare 2 vector (output from the previous 3 functions), also with the maximum distance as the third input. Scan through both vectors and find their distance. If smaller than the 3rd input, then that means they should be together. 

novel_TSS = read.csv("novel_TSS_both.csv")

promoter_regions = c()
for(position in 1:length(TSS_positions)){
  for(index in 1:length(fasta_sequence_forward)){#the position of the nucleotide in fasta
    if(TSS_positions[position] == index){#if TSS equals the position of the nucleotide
      if(TSS_strand[position] == "+"){
        promoter = (fasta_sequence_forward[(index-49):index])
        promoter_regions = c(promoter_regions, promoter)
      }
      if(TSS_strand[position] == "-"){
      
        promoter = reverseComplement((fasta_sequence_forward[index:(index+49)]))
        promoter_regions = c(promoter_regions, promoter)
      }
      

    }
  }
}

promoter_string = as(promoter_regions, "DNAStringSet")#needs to be converted to DNAstringset first
promoter_vec = rep(NA, length(TSS_positions))
for(i in 1:length(TSS_positions)){
  promoter_vec[i] = as.character(promoter_string[i])
}
cap_TSS$promoters = promoter_vec
write.csv(cap_TSS, file="TSS_sigma_factors.csv")

writeXStringSet(promoter_string, file="promoter_regions_rc.fa")
promoter = read.delim("promoter_regions_rc.fa", header = F)
promoters = as.character(promoter$V1)
output.df = data.frame(stringsAsFactors = FALSE)
id =1
for(line in 1:length(promoters)){
  
  #print(promoters[line])
    if(line %% 2 != 0){
      t = promoters[line]
      t = paste0(t, "promoter",id)
      print(id)
      id=id+1
      
    }
  else{
    t = promoters[line]
  }
  output.df=rbind(output.df, t, stringsAsFactors=FALSE)
}


write.csv(output.df, file="promoter_regions_rc_id.fa")

################################################################################################
###### fasta file for common TSS
common_TSS = read.csv("common_TSS.csv", header = F)
common_positions = common_TSS$V1
common_strand = common_TSS$V2
promoter_regions = c()
for(position in 1:length(common_positions)){
  for(index in 1:length(fasta_sequence_forward)){#the position of the nucleotide in fasta
    if(common_positions[position] == index){#if TSS equals the position of the nucleotide
      if(common_strand[position] == "+"){
        promoter = (fasta_sequence_forward[(index-49):index])
        promoter_regions = c(promoter_regions, promoter)
      }
      if(common_strand[position] == "-"){
       
        promoter = rev((fasta_sequence_forward[index:(index+49)]))
        promoter_regions = c(promoter_regions, promoter)
      }
      
      
    }
  }
}

promoter_string = as(promoter_regions, "DNAStringSet")
writeXStringSet(promoter_string, file="promoter_regions_common.fa")
promoter = read.delim("promoter_regions_common.fa", header = F)
promoters = as.character(promoter$V1)
output.df = data.frame(stringsAsFactors = FALSE)
id =1
for(line in 1:length(promoters)){
  
  #print(promoters[line])
  if(line %% 2 != 0){
    t = promoters[line]
    t = paste0(t, "promoter",id)
    print(id)
    id=id+1
    
  }
  else{
    t = promoters[line]
  }
  output.df=rbind(output.df, t, stringsAsFactors=FALSE)
}


write.csv(output.df, file="promoter_regions_common_id.fa")

################################################################################################
###### fasta file for Arnaud supplementary file
supp_TSS = read.csv("Arnaud_supplementary_promoter.csv", header = T)
pub_positions = supp_TSS$published_TSS
dRNA_positions = supp_TSS$dRNA_TSS
supp_strand = supp_TSS$strand

promoter_regions_pub = c()
promoter_regions_dRNA = c()
for(position in 1:length(pub_positions)){
  for(index in 1:length(fasta_sequence_forward)){#the position of the nucleotide in fasta
    if(pub_positions[position] == index){#if TSS equals the position of the nucleotide
      if(supp_strand[position] == "+"){
        promoter_pub = (fasta_sequence_forward[(index-49):index])
        promoter_regions_pub = c(promoter_regions_pub, promoter_pub)
      }
      if(supp_strand[position] == "-"){
        promoter_pub = rev((fasta_sequence_forward[index:(index+49)]))
        promoter_regions_pub = c(promoter_regions_pub, promoter_pub)
      }
      
      
    }
    if(dRNA_positions[position] == index){#if TSS equals the position of the nucleotide
      if(supp_strand[position] == "+"){
        promoter_dRNA = (fasta_sequence_forward[(index-49):index])
        promoter_regions_dRNA = c(promoter_regions_dRNA, promoter_dRNA)
      }
      if(supp_strand[position] == "-"){
        promoter_dRNA = rev((fasta_sequence_forward[index:(index+49)]))
        promoter_regions_dRNA = c(promoter_regions_dRNA, promoter_dRNA)
      }
      
      
    }
  }
}

promoter_string_pub = as(promoter_regions_pub, "DNAStringSet")
promoter_string_dRNA = as(promoter_regions_dRNA, "DNAStringSet")
writeXStringSet(promoter_string_pub, file="promoter_regions_pub.fa")
writeXStringSet(promoter_string_dRNA, file="promoter_regions_dRNA.fa")
promoter_pub = read.delim("promoter_regions_pub.fa", header = F)
promoter_dRNA = read.delim("promoter_regions_dRNA.fa", header = F)
promoters_pub = as.character(promoter_pub$V1)
promoters_dRNA = as.character(promoter_dRNA$V1)
output.df_pub = data.frame(stringsAsFactors = FALSE)
output.df_dRNA = data.frame(stringsAsFactors = FALSE)
id =1
for(line in 1:length(promoters_pub)){
  
  #print(promoters[line])
  if(line %% 2 != 0){
    t = promoters_pub[line]
    t = paste0(t, "promoter",id)
    s = promoters_dRNA[line]
    s = paste0(s, "promoter",id)
    print(id)
    id=id+1
    
  }
  else{
    t = promoters_pub[line]
    s = promoters_dRNA[line]
  }
  output.df_pub=rbind(output.df_pub, t, stringsAsFactors=FALSE)
  output.df_dRNA=rbind(output.df_dRNA, s, stringsAsFactors=FALSE)
}


write.csv(output.df_pub, file="promoter_regions_pub_id.fa")
write.csv(output.df_dRNA, file="promoter_regions_dRNA_id.fa")

####################################################################################################################
sigma70pwm = read.delim("sigma70_pssm.txt", header = F)
sigma70matrix <- rbind(A=(sigma70pwm$V1)[24:29], C=(sigma70pwm$V2)[24:29], G=(sigma70pwm$V3)[24:29], T=(sigma70pwm$V4)[24:29]) #motif is the last 6 TAnnAT from MEME downloaded the count matrix and delimited by tab in excel, columns in alphabetical order ACGT
forward70<- matchPWM(sigma70matrix, fasta_sequence_forward, min.score="90%")
reverse70 <- matchPWM(sigma70matrix, fasta_sequence_reverse,min.score="90%")

sigma54pwm = read.delim("sigma54_pssm.txt", header = F)
sigma54matrix <- rbind(A=(sigma54pwm$V1)[6:22], C=(sigma54pwm$V2)[6:22], G=(sigma54pwm$V3)[6:22], T=(sigma54pwm$V4)[6:22])
forward54<- matchPWM(sigma54matrix, fasta_sequence_forward, min.score="90%")
reverse54 <-matchPWM(sigma54matrix, fasta_sequence_reverse, min.score="90%")
