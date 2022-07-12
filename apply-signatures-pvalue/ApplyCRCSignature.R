library(tidyverse)
library(patchwork)
library(ggplot2)
library(forcats)

#set working directory
setwd()
# FUNCTIONS


#Helper function - Resample genes - Takes dataframe with 2 columns - gene and groups
resampling_method <- function(signatures){
  #Turn signatures into list
  l <- length(unique(signatures$group ) )
  gene_signatures <- vector("list", l)
  for (i in 1:l){
    names(gene_signatures)[i] <- unique(signatures$group)[i]
    gene_signatures[[i]] <- signatures$gene[signatures$group == unique(signatures$group)[i] ]
  }
  
  #Establish variables
  rs_gs <- vector(mode = "list", length = length(gene_signatures))
  names(rs_gs) <- names(gene_signatures)
  
  #Unlist genes from all gene signatures
  ##gene_wide
  genes <- signatures$gene %>% as.vector()
  
  #Resampling  - From genes
  for (i in 1:length(gene_signatures)){
    rs_gs[i][[1]] <- sample(genes, size = length(gene_signatures[i][[1]]), replace = T )
    
    #Adjust for duplicated values in dataset
    #Note: This part fails if there are not enough unique values for every index
    while (any (duplicated (rs_gs[i]) ) == T ){
      dup <- which( rs_gs[[i]] == rs_gs[[i]][duplicated( rs_gs[[i]] )] [1] )
      rs_gs[[i]][dup[2]] <- sample(genes[!genes %in% rs_gs[[i]]], size = 1, replace = T)
    }
    
    
  }
  
  new_df <- data.frame(gene = unlist(rs_gs, use.names = F), 
                       group = names(setNames(unlist(rs_gs, use.names = F), 
                                              rep(names(rs_gs), lengths(rs_gs)))))
  
  #return 
  return(new_df)
}




## main function, sig: signatures to apply; data: copy number data
ApplySignature <-function(sig, data){
  sig_n <- count(sig, group)
  tib_wide <- inner_join(sig, data)
  tib <- tib_wide %>% gather(sample_id, value, -group, -gene)
  
  # Row normalization 
  CalZscore <- function(value){(value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)}
  tib_norm <- tib %>%
    ## scale each gene  across all samples
    group_by(gene) %>%
    mutate(value = CalZscore(value)) %>%
    ungroup()
  
  # Calculate scores 
  tib_norm_mean <- tib_norm %>%
    group_by(group, sample_id) %>%
    summarise(value = mean(value, na.rm = TRUE)) %>%
    spread(sample_id, value)
  mat <- column_to_rownames(tib_norm_mean, "group") 
  v <- column_to_rownames(sig_n, "group") %>% sqrt() %>% pull(n)
  score <- mat * v
  score_tib <- as_tibble(score, rownames = "group") %>%
    gather(sample, score, -group) %>%
    mutate(sample = factor(sample, levels = str_sort(unique(sample), numeric = TRUE))) %>%
    mutate(group = factor(group, levels = c(str_sort(unique(group), numeric = TRUE), "unknown"))) %>%
    arrange(sample) 
  
  score_tib$group[score_tib$score<0] <- "unknown"
  
  assignment <- 
    score_tib %>% group_by(sample) %>% top_n(1, score) %>% 
    ungroup()  %>% rename(case_id = sample) 
  
  return(list(score = score_tib, max_score = assignment))
}

PlotDots <- function(scores){
  p <- scores %>% 
    ggplot(aes(x=sample, y = score, color = group, size = score, alpha = score)) +
    geom_point() +
    geom_point(color = "black", shape = 1) +
    facet_wrap("sample", nrow = 1, scales = "free_x")+
    theme_bw()
  return(p)
}

PlotClinical <- function(assignment, feature_column_name){
  assignment_count <- assignment %>% count(group, {{feature_column_name}})
  p <- ggplot(
    assignment_count, 
    aes(x = fct_rev(group), y = n, fill = {{feature_column_name}},  label = n)
  ) +
    geom_bar(stat = "identity") +
    geom_text(position = position_stack(vjust = 0.5) , fontface = "bold", size = 5)+
    coord_flip()+
    theme_bw()+
    theme(axis.text  = element_text(size = 12, face = "bold"),
          legend.position = "right")+
    xlab("Group")
  
  
  return(p)
}


p_value <- function(sig, data, num_resample){
  set.seed(123)
  
  #Obtain original scores
  og_s <- ApplySignature(sig, data)
  score <- as.data.frame(og_s$max_score$score, row.names = colnames(data)[2:ncol(data)])
  
  #Initialize Counter
  p_v_counter <- data.frame(row.names = rownames(score))
  p_v_counter[1:nrow(p_v_counter), 1] <- 0
  
  #Initialize resampling list
  LoR <- vector("list", num_resample)
  #Make list
  for (i in 1:num_resample){
    LoR[[i]] <- resampling_method(sig)
  }
  
  #initialize scores
  scores_rs <- vector("list", num_resample)
  
  #Make resampled scores and acquire p_value counter
  for(i in 1:num_resample){
    scores_rs[[i]] <- ApplySignature(LoR[[i]], data )
    
    #Comparison - add to p_v_counter
    s_of_i <- as.data.frame(scores_rs[[i]]$max_score$score)
    
    df1 <- 1*(s_of_i >= score)
    
    p_v_counter <- p_v_counter + df1
  }
  
  #Calculation of p-value
  p_value <- (p_v_counter/num_resample)
  
  #Obtain all max scores across all resamples
  max_scores_rs <- as.data.frame(matrix(nrow = num_resample, ncol = (ncol(data)-1)))
  colnames(max_scores_rs) <- colnames(data)[2:ncol(data)]
  
  for (i in 1:num_resample){
    max_scores_rs[i, 1:ncol(max_scores_rs)] <- sapply(scores_rs[[i]][[2]][[3]], "[")
  }
  
  #returns p-value, max scores of resamples, and max scores of original
  return(list(p_value, max_scores_rs, og_s[[2]][[3]]))
  
}

#Visualize all scores, including resample. Takes in p_value list from above

v_score <- function(all_scores){
  #turn scores into dataframe that works with ggplot
  df <- pivot_longer(all_scores[[2]], cols = 1:ncol(all_scores[[2]]), names_to = "Sample", values_to = "Scores")
  
  #Ensures order of samples
  df <- df %>% mutate(Sample = fct_relevel(Sample, unique(df$Sample)))
  
  #get number of samples
  l_sam <- length(unique(df$Sample))
  
  #min maxing
  md <- min(unlist(c(all_scores[[2]][1:100,], all_scores[[3]])))
  mu <- max(unlist(c(all_scores[[2]][1:100,], all_scores[[3]])))
  
  #acquire segments based on number of samples
  seg <- data.frame(stx = seq(0.75, (l_sam - 0.25) ), sty = all_scores[[3]], enx = seq(1.25, (l_sam + 0.25)), eny = all_scores[[3]] )
  
  
  #make the plot
  v <- ggplot(df, aes(x = Sample, y = Scores)) +
    geom_point() +
    geom_segment(data = seg, x = seg$stx, y = seg$sty, xend = seg$enx, yend = seg$eny, color = "red") +
    ylim(md, mu)
  
  
  #return the plot
  return(v)
  
}

# 1 Read CRC signatures ----
sig_data_p <- read_tsv("data/primary_signature.tsv")
sig_data_m <- read_tsv("data/metastatic_signature.tsv")
sig_data_pm <- bind_rows(sig_data_p,sig_data_m) %>% 
  mutate(group = str_remove_all(group, '[0-9]+')) # 


# 2 Read patient info and copy number data ----
## patient info 
clin <- read_tsv("data/GSE100787 Primary Matched Metastasis/From Dennis Wylie/Features_gse100787_enhanced.tsv" )
clin_list <- split(clin, clin$sample_type) 
## copy number
cn_path <- "data/GSE100787 Primary Matched Metastasis/From Dennis Wylie/gse100787_genewise_most_extreme_segmented_cn_values.tsv.gz"
cn_data <- vroom::vroom(cn_path) %>% rename(gene = gene_name) %>% 
  filter(!duplicated(gene)) # remove duplicated gene from raw data
## primary & metastatic
cn_data_pm <- cn_data  %>%  
  select(gene, clin$sample) %>% # select sample columns
  purrr::set_names(c("gene", clin$case_id)) # rename "GSM" with shorter case id
cn_data_p <- cn_data_pm %>% select(gene, ends_with("p")) # primary
cn_data_m <- cn_data_pm %>% select(gene, ends_with("m")) # Metastatic

# 3 Running ----
scores_p <- ApplySignature(sig = sig_data_p, data = cn_data_p)
scores_m <- ApplySignature(sig_data_m, cn_data_m)

clin_p <- scores_p[["max_score"]] %>% rename(cluster = group) %>% left_join(clin)
clin_m <- scores_m[["max_score"]] %>% rename(cluster = group) %>% left_join(clin)
scores_p[[2]]
# Extact the signature with the max score ----
assigned_sig <- scores_p[[2]] %>% 
  left_join(clin)  # combine signature assignment and clinical info for each sample
## Viusalize all the scores 
PlotDots(scores_p[[1]]) 
## Viusalize clinial information
PlotClinical(assigned_sig, mss) + 
  PlotClinical(assigned_sig, kras_wt) +
  PlotClinical(assigned_sig, braf_wt) +
  PlotClinical(assigned_sig, left_sided) +
  PlotClinical(assigned_sig, outcome) +
  PlotClinical(assigned_sig, liver_metastasis_synchronous_metachronous_with_primary) 


# 2.1 read copy number data for pre/post
gse_m <- read.csv("gse152178/GSE152178_CN_gene_level.csv")
feat_gse <- readxl::read_xlsx("gse152178/sample_info_CTM2-11-e401-s003_45 Pre-Post Samples with CN data.xlsx")
feat_gse$`Sample name` <- gsub("\\-", ".", feat_gse$`Sample name`)

#Get rid of duplicates (one of the duplicates had no copy number)
dp <- feat_gse$Patient[duplicated(feat_gse$Patient)  ]

feat_gse <- feat_gse[feat_gse$Patient %in% dp, ]
feat_gse <- feat_gse[-13,]

#separate by pre and post treatment
feat_PRE <- feat_gse[feat_gse$`Pre or Post-treatment for LM` == "Pre", ]
feat_Post <- feat_gse[feat_gse$`Pre or Post-treatment for LM` != "Pre", ]

feat_PRE$case_id <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
feat_Post$case_id <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")

#finding samples
gse_pre <- gse_m[, substr(colnames(gse_m), 12, 21) %in% substr(feat_PRE$`Sample name`, 1, 10) ]
gse_post <- gse_m[, substr(colnames(gse_m), 12, 21) %in% substr(feat_Post$`Sample name`, 1, 10) ]
gse_post<- gse_post[, c(1, 2, 3, 4, 5, 8, 6, 7)]

#add column (sample) and rows (gene name)
colnames(gse_pre) <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
colnames(gse_post) <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")

gse_pre_cn <- add_column(gse_pre, gene = gse_m$gene, .before = 1)
gse_post_cn <- add_column(gse_post, gene = gse_m$gene, .before = 1)

#make p_value for samples
p_pre <- p_value(sig_data_m, gse_pre_cn, 100)

p_post <- p_value(sig_data_m, gse_pre_cn, 100)

#get visual of scores
g <- v_score(p_pre)
g_po <- v_score(p_post)

g
g_po


#original scores
scores_pre <- ApplySignature(sig_data_m, gse_pre_cn)
scores_post <- ApplySignature(sig_data_m, gse_post_cn)

PlotDots(scores_pre[[1]])
PlotDots(scores_post[[1]])

a_sig_pre <- scores_pre[[2]] %>% left_join(feat_PRE)
a_sig_post <- scores_post[[2]] %>% left_join(feat_Post)

PlotClinical(a_sig_pre, group) +
  PlotClinical(a_sig_post, group)
















