
# This script creates an exclusion criteria file for MIST analysis. 
# In file should have two columns: "Sample" and "Batch"
# Batches should number samples 1-N, N being the number of batches of things to apply exclusion criteria to. 0 will be ignored.

# Example in_file:

# Sample	Batch
# nsp3_WT	0
# S_WT	1
# S_ManausP1	1
# S_B.1.351	1
# N_WT	2
# N_P80R	2
# N_T205I	2
# N_H145Y	2
# N_R203K	2
# N_A220V	2
# E_WT	3
# E_P71L 	3
# nsp2_WT	4
# nsp2_I121F	4
# nsp13_WT	5
# nsp13_E342D	5
# orf8_WT	6
# orf8_E92K	6
# orf8_nature	6
# N_nature	2
# S_nature	1
# nsp13_nature	5
# nsp2_nature	4
# E_nature	3

# Written by M. Bouhaddou 2021-04-12

make_exclusion_criteria_MIST = function(in_file){
  
  #Libraries
  library(data.table)
  library(gtools)
  
  # Load data
  D = fread(in_file)
  nums = unique(D$Batch)[unique(D$Batch)>0]
  D_final <- data.frame()
  
  for (i in nums){
    samples = D$Sample[D$Batch==nums[i]]
    
    a <- data.frame(V1 = samples)
    
    b <- sapply(samples, function(x){
      x1 <- samples[-which(samples == x)]
      return(x1)
    })
    
    if(is.null(nrow(b))){
      a <- as.matrix(cbind(a, b))
    }else{
      a <- as.matrix(cbind(a, t(b)))
    }
    
    # a=permutations(n = length(samples), r = length(samples), v = samples, repeats.allowed = FALSE)
    
    # Remove duplications
    # a = a[!duplicated(a[,1]),]
    
    vec = rep(NA,nrow=nrow(a))
    for (n in 1:nrow(a)){
      vec[n]=paste(a[n,2:ncol(a)],collapse="|")
    }
    
    # Create final data.frame
    if (which(nums == i)==1){
      D_final = data.frame(V1 = a[,1],V2 = vec)
    } else {
      D_final = rbind(D_final,data.frame(V1 = a[,1],V2 = vec))
    }
  }
  return(D_final)
}

# fwrite(D_final,file="exclusion_criteria_Ready.txt",sep='\t',col.names=F)


