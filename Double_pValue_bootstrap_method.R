# Official Double p-value script
# Jorge Bustamante

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(edgeR))

#Resample genes, requires gene_signatyures as list, and data matrix
resampling_method <- function(gene_signatures, ge_data){
  #Establish variables
  rs_gs <- vector(mode = "list", length = length(gene_signatures))
  names(rs_gs) <- names(gene_signatures)
  
  rs_gs_gw <- vector(mode = "list", length = length(gene_signatures))
  names(rs_gs_gw) <- names(gene_signatures)
  
  #Unlist genes from all gene signatures
  ##gene_wide
  genes <- unlist(gene_signatures) %>% as.vector()
  
  ##genome_wide
  g_w <- rownames(ge_data)
  
  #Resampling  - From genes and g_w
  for (i in 1:length(gene_signatures)){
    rs_gs[i][[1]] <- sample(genes, size = length(gene_signatures[i][[1]] ))
    
    #Adjust for duplicated values in dataset
    #Note: This part fails if there are not enough unique values for every index
    while (any (duplicated (rs_gs[i]) ) == T ){
      dup <- which( rs_gs[[i]] == rs_gs[[i]][duplicated( rs_gs[[i]] )] [1] )
      rs_gs[[i]][dup[2]] <- sample(genes[!genes %in% rs_gs[[i]]], size = 1, replace = F)
    }
    
    rs_gs_gw[i][[1]] <- sample(g_w, size = length(gene_signatures[i][[1]] ))
  }
  
  
  #return both in a single list
  return(list(rs_gs, rs_gs_gw))
}


#Set seed
set.seed(123)

#Insert matrix and list of gene signatures

#####################################

# gc_mat <- #Insert matrix you wish to use here
##Note: gc_mat should have unique rownames since it is a matrix

# ls <- #Insert list of gene signatures that will be used


#####################################

#Establish normal gsva scores
og_gsva <- gsva(gc_mat, ls, method = "gsva")

#Establish number of resamplings wanting to be performed
n_rs <- 500

#Initialize counters for pvalues

p_v_counter_g <- matrix(0, nrow = nrow(og_gsva), ncol = ncol(og_gsva), dimnames = list(rownames(og_gsva), colnames(og_gsva))) %>% as.data.frame()

p_v_counter_genome <- matrix(0, nrow = nrow(og_gsva), ncol = ncol(og_gsva), dimnames = list(rownames(og_gsva), colnames(og_gsva))) %>% as.data.frame()

#Resamples
for (i in 1:n_rs){
  cat(i, "\n")
  
  #Get new gene signatures
  rs_gene_sig <- resampling_method(ls, gc_mat)
  
  
  #Get new scores, both gene_sig and genomewide
  new_gsva <- gsva(gc_mat, rs_gene_sig[[1]])
  
  new_gsva_genome <- gsva(gc_mat, rs_gene_sig[[2]])
  
  #Comparison counter
  df1 <- 1*(abs(new_gsva) >= abs(og_gsva))
  
  p_v_counter_g <- p_v_counter_g + df1
  
  df2 <- 1*(abs(new_gsva_genome) >= abs(og_gsva))
  
  p_v_counter_genome <- p_v_counter_genome + df2
  
  
  
}

#P-value Calculation

pvalues_g <- (p_v_counter_g + 1) / (n_rs + 1)

pvalues_genome <- (p_v_counter_genome + 1) / (n_rs + 1)

pvalues_g <- pvalues_g %>% round(digits = 3)

pvalues_genome <- pvalues_genome %>% round(digits = 3)

#Save files

pvalues_g %>% write.csv(file = "p_value_gene_sig_500bootstrap.csv")
pvalues_genome %>% write.csv(file = "p_value_genomewide_500bootstrap.csv")
