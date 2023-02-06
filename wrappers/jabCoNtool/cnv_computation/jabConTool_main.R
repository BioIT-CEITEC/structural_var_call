suppressMessages(library(data.table))
suppressMessages(library(stringi))
suppressMessages(library(zoo))
suppressPackageStartupMessages(library(proxy))
library(fitdistrplus)
library(ggfortify)
library(parallel)
library(pheatmap)
library(stats)
library(dplyr)
library(ClusterR)
library(mclust)

# setwd("/mnt/ssd/ssd_1/workspace/vojta/cfDNA_PCa_Kluge/")

source(paste0(script_dir,"/jabConTool_func_load_inputs.R"))
source(paste0(script_dir,"/estimate_tumorload.R"))

#TODO sex chromosome estimation

#set_constants
sample_regex <<- ".*/(.*)/jabCoNtool.*"


min_corelation_threshold <<- 0.9
dist_transition_treshold <<- 50000


###########################
###########################
######               ######
######   SNP part    ######
######               ######
###########################
###########################


compute_snp_based_nloglike <- function(snp_tab){
  
  # if(all(snp_tab$TL == 1)){
  #   cn_het_var_count_table <- cn_het_var_count_table[-1]
  # }
  
  snp_based_cn_nloglike <- data.table(id = rep(seq_along(snp_tab$chr),each = nrow(cn_het_var_count_table)),
                                      TL = rep(snp_tab$TL,each = nrow(cn_het_var_count_table)),
                                      alt_count = rep(snp_tab$alt_count,each = nrow(cn_het_var_count_table)),
                                      ref_count = rep(snp_tab$ref_count,each = nrow(cn_het_var_count_table)),
                                      cn_category = rep(cn_het_var_count_table$cn_category,nrow(snp_tab)),
                                      het_var_count = rep(cn_het_var_count_table$het_var_count,nrow(snp_tab)),
                                      cn_count = rep(cn_het_var_count_table$cn_count,nrow(snp_tab)))
  snp_based_cn_nloglike[,expected_het_ratio := ((1 - TL) + TL * het_var_count) / (2 * (1 - TL) + TL * cn_count)]
  snp_based_cn_nloglike[,snp_prob := dbeta(expected_het_ratio,alt_count + 1,ref_count + 1)]
  snp_based_cn_nloglike[,snp_nloglike := -log(snp_prob)]
  snp_based_cn_nloglike <- snp_based_cn_nloglike[,.(snp_nloglike = min(snp_nloglike)),by = .(id,cn_category)]
  
  if(any(snp_based_cn_nloglike$TL == 1)){
    snp_based_cn_nloglike[TL == 1,mean_snp_nloglike := mean(snp_nloglike,na.rm = T),by = id]
    snp_based_cn_nloglike[TL == 1 & cn_category == "0",snp_nloglike := mean_snp_nloglike]
  }  
  
  snp_tab <- cbind(snp_tab,dcast.data.table(snp_based_cn_nloglike,formula = id ~ cn_category,value.var = "snp_nloglike"))
  snp_tab[,id := NULL]
  snp_tab[ , (cn_categories_vec) := lapply(.SD, "+", -log(pop_HET_probability)), .SDcols = cn_categories_vec]
  
  
  return(snp_tab)
}



###########################
###########################
######               ######
###### COVRAGE part  ######
######               ######
###########################
###########################


vector.weighted.mean <- function(value_vec,weight_matrix){
  return(colSums(value_vec * weight_matrix) / colSums(weight_matrix))
}

vector.weighted.sd <- function(value_vec,weight_matrix){
  mean_dif_mat <- matrix(value_vec - rep(vector.weighted.mean(value_vec,weight_matrix),each = length(value_vec)),ncol = length(value_vec))^2
  non_zero_weight_count <- rowSums(matrix(as.logical(weight_matrix),ncol = length(value_vec)))
  
  return(sqrt(colSums(mean_dif_mat * weight_matrix) / (((non_zero_weight_count - 1) / non_zero_weight_count)*colSums(weight_matrix))))
}

compute_distribution_pramaters <- function(cov_tab,library_type = "panel"){
  
  
  if(any(names(cov_tab) == "cov_norm_factor")){
    cov_tab[,norm_vec := TL * cov_norm_factor + (1 - TL)]
  } else {
    cov_tab[,norm_vec := 1]
  }
  
  if(library_type != "wgs"){
    mat <- dcast.data.table(cov_tab,formula = region_id ~ sample,value.var = "cov",fill = 0)
    mat_row_name <- mat$region_id
    mat <- as.matrix(mat[,-1,with = F])
    rownames(mat) <- mat_row_name
    cor_mat <- as.matrix(simil(mat,method = cor, by_rows = F,use = "complete.obs"))  
    cor_mat <- cor_mat - min_corelation_threshold
    cor_mat[cor_mat < 0] <- 0
    cor_mat <- cor_mat / (1 - min_corelation_threshold)
    diag(cor_mat) <- 1
    
    cov_tab[,norm_dist_mean := vector.weighted.mean(cov,cor_mat),by = "region_id"]
    cov_tab[,norm_dist_sd := vector.weighted.sd(cov,cor_mat),by = "region_id"]
    
  } else {
    normalization_sample_type <- tail(sort(unique(cov_tab$type)),1)
    norm_nbinom_dist <- fitdist(as.integer(cov_tab[type == normalization_sample_type]$cov / cov_tab[type == normalization_sample_type]$norm_vec),distr = "nbinom")
    
    mu <- norm_nbinom_dist$estimate["mu"]
    size <- norm_nbinom_dist$estimate["size"]
    cov_tab[,nbinom_mean := mu]
    cov_tab[,nbinom_var := mu + mu^2/size]
    
  }
  
  cov_tab[,norm_vec := NULL]
  
  return(cov_tab)
  
}

#TODO implement differently with binom and byas 
get_normal_negative_log_likelihoods <- function(cn_count,cov_tab,library_type){
  cn_count_to_normal <- cn_count / 2
  
  mean <- cov_tab$TL * (cn_count_to_normal * cov_tab$norm_dist_mean) + cov_tab$norm_dist_mean * (1 - cov_tab$TL)
  sd <- cov_tab$TL * (cn_count_to_normal * cov_tab$norm_dist_sd) + cov_tab$norm_dist_sd * (1 - cov_tab$TL)  
  if(cn_count_to_normal == 0){
    mean[cov_tab$TL == 1] <- 1
    sd[cov_tab$TL == 1] <- cov_tab$norm_dist_sd / 4
  }
  nloglike <- -log(1 - abs(pnorm(cov_tab$cov,mean = mean,sd = sd) - 0.5))
  return(nloglike)
}

#TODO implement differently with binom and byas 
get_nbinom_negative_log_likelihoods <- function(cn_count,cov_tab){
  cn_count_to_normal <- cn_count / 2
  
  mu <- cov_tab$TL * (cn_count_to_normal * cov_tab$nbinom_mean) + cov_tab$nbinom_mean * (1 - cov_tab$TL)
  var <- cov_tab$TL * (cn_count_to_normal * cov_tab$nbinom_var) + cov_tab$nbinom_var * (1 - cov_tab$TL)  
  if(cn_count_to_normal == 0){
    mu[cov_tab$TL == 1] <- 1
    var[cov_tab$TL == 1] <- cov_tab$nbinom_var / 4
  }
  
  size <- mu^2 / (var - mu) 
  nloglike <- -log(dnbinom(round(cov_tab$cov),mu = mu,size = size))
  return(nloglike)
}

compute_coverage_based_nloglike <- function(cov_tab,library_type,cn_categories_tab){
  
  cov_tab <- compute_distribution_pramaters(cov_tab,library_type)
  
  if(library_type != "wgs"){
    res_nloglike <- sapply(cn_categories_tab$cn_count,get_normal_negative_log_likelihoods,cov_tab = cov_tab)
  } else {
    res_nloglike <- sapply(cn_categories_tab$cn_count,get_nbinom_negative_log_likelihoods,cov_tab = cov_tab)
  }
  
  colnames(res_nloglike) <- cn_categories_tab$cn_category
  cov_tab <- cbind(cov_tab,res_nloglike)
  cov_tab <- cov_tab[type == "call"]
  
  
  setkey(cov_tab,chr,start,end,sample)
  
  return(cov_tab)
}


###########################
###########################
######               ######
###### bayes part   ######
######               ######
###########################
###########################


compute_byes_net_likelihood <-function(trans_mat_list,per_region_state_nloglike_matrix,only_model_nloglike = T,cn_categories_vec){
  
  # per_region_state_nloglike_matrix = combined_nloglike_matrix
  # only_model_nloglike = estimate_only
  
  n_regions <- length(trans_mat_list)
  prev = matrix(0, n_regions-1, length(cn_categories_vec))
  omega = matrix(0, n_regions, length(cn_categories_vec))
  
  
  omega[1,] = trans_mat_list[[1]][1,] + per_region_state_nloglike_matrix[1,]
  
  index_vec <- seq_along(trans_mat_list)[-1]
  
  trans_index <- 1:ncol(trans_mat_list[[1]]) - 1
  trans_mat_cols <- ncol(trans_mat_list[[1]])
  
  for(i in index_vec){
    probs_mat <- omega[i - 1, ] + trans_mat_list[[i]]
    prev[i - 1, ] <- apply(probs_mat,2,which.min)
    omega[i, ] <- probs_mat[prev[i - 1, ] + trans_index * trans_mat_cols] + per_region_state_nloglike_matrix[i,]
  }
  
  if(only_model_nloglike){
    return(min(omega[nrow(omega),]))
  }
  
  predict_state_vector <- integer(length = n_regions)
  predict_state_vector[length(predict_state_vector)] <- which.min(omega[nrow(omega),])
  
  
  for(i in (length(predict_state_vector)-1):1){
    predict_state_vector[i] <- prev[i,predict_state_vector[i + 1]]
  }
  
  return(list(min(omega[nrow(omega),]),predict_state_vector))
}

###########################
###########################
######               ######
###### wrapper part  ######
######               ######
###########################
###########################

predict_CNV_model <- function(sample_tab,trans_mat_list,cov_tab,snp_tab,library_type,categories_default_tabs,TL_select = NULL,estimate_only = F){
  
  # sample_tab <- process_sample_tab
  sample_tab <- copy(sample_tab)
  
  if(!is.null(TL_select)){
    sample_tab[,TL := sample_tab[[TL_select]]]
  } else {
    if(!any(names(sample_tab) == "TL")){
      sample_tab[,TL := 1]
    }
  }
  
  tictoc::tic()
  print(paste0("cov + snp_prepar: "))
  
  if(any(names(cov_tab) == "cn_pred")){
    cov_tab <- merge(sample_tab[,.(sample,type,TL_new)],cov_tab,by = "sample")
    cov_tab[cn_pred == "2",TL := TL_new]
    cov_tab[,TL_new := NULL]
    cov_tab <- merge(cov_tab,categories_default_tabs$cn_categories_tab[,.(cn_pred = cn_category,cov_norm_factor)],by = "cn_pred")
  } else {
    cov_tab <- merge(sample_tab[,.(sample,type,TL)],cov_tab,by = "sample")
  }
  cov_nloglike_tab <- compute_coverage_based_nloglike(cov_tab,library_type,categories_default_tabs$cn_categories_tab)
  
  if(!is.null(snp_tab)){
    if(any(names(snp_tab) == "cn_pred")){
      snp_tab <- merge(sample_tab[,.(sample,type,TL_new)],snp_tab,by = "sample")
      snp_tab[cn_pred == "2",TL := TL_new]
      snp_tab[,TL_new := NULL]
      
    } else {
      snp_tab <- merge(sample_tab[,.(sample,TL)],snp_tab,by = "sample")
    }
  }
  
  tictoc::toc()
  
  cn_categories_vec <- categories_default_tabs$cn_categories_tab$cn_category
  
  per_sample_res <- lapply(sample_tab[type == "call"]$sample,function(sel_sample){
    
    #sel_sample = "141"
    
    if(!is.null(snp_tab)){
      tictoc::tic()
      print(paste0("sample: ",sel_sample))
      snp_nloglike_tab <- compute_snp_based_nloglike(snp_tab[sample == sel_sample])
      tictoc::toc()
      
      tictoc::tic()
      combined_nloglike_tab <- rbind(cov_nloglike_tab[sample == sel_sample,c("region_id",cn_categories_vec),with = F],snp_nloglike_tab[,c("region_id",cn_categories_vec),with = F])
      setkey(combined_nloglike_tab,region_id) 
      combined_nloglike_tab <- melt.data.table(combined_nloglike_tab,id.vars = "region_id")
      combined_nloglike_tab <- combined_nloglike_tab[,.(value = sum(value)),by = .(region_id,variable)]
      combined_nloglike_matrix <- matrix(data = combined_nloglike_tab$value,nrow = length(trans_mat_list))
      tictoc::toc()
    } else {
      combined_nloglike_matrix <- as.matrix(cov_nloglike_tab[sample == sel_sample,cn_categories_vec,with = F])
    }
    tictoc::tic()
    res <- compute_byes_net_likelihood(trans_mat_list = trans_mat_list,
                                       per_region_state_nloglike_matrix = combined_nloglike_matrix,
                                       only_model_nloglike = estimate_only,
                                       cn_categories_vec = cn_categories_vec)
    
    tictoc::toc()
    return(res)
  })
  
  if(estimate_only){
    TL_estimate_tab <- sample_tab[type == "call",c("sample",TL_select),with = F]
    TL_estimate_tab[,(paste0(TL_select,"_nloglike")) := unlist(per_sample_res)]
    
    return(TL_estimate_tab)
  } else {
    model_estimate_tab <- sample_tab[type == "call",.(sample,TL,nloglike = sapply(per_sample_res,function(x) x[[1]]))]
    cn_predict_tab <- lapply(per_sample_res,function(x) data.table(region_id = unique(cov_tab$region_id),cn_pred = x[[2]]))
    names(cn_predict_tab) <- sample_tab[type == "call"]$sample
    cn_predict_tab <- rbindlist(cn_predict_tab,use.names = T,idcol = "sample")
    
    return(list(model_estimate_tab,cn_predict_tab))
  }
}

estimate_per_sample_tumor_load <- function(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list){
  
  tumor_load_precision <- 0.02
  iterations <- 1:floor(log2(1 / tumor_load_precision))
  
  sample_tab[,high_TL := 0.75]
  sample_tab[,low_TL := 0.25]
  res_sample_tab_list <- list()
  
  for(i in iterations){
    tictoc::tic()
    print(paste0("iter: ",i))
    res_list <- lapply(c("high_TL","low_TL"),function(x) predict_CNV_model(sample_tab,trans_mat_list,cov_tab,snp_tab,library_type,TL_select = x,estimate_only = T))
    res_sample_tab <- merge(res_list[[1]],res_list[[2]],by = "sample")
    res_sample_tab[,best_TL := high_TL]
    res_sample_tab[low_TL_nloglike < high_TL_nloglike,best_TL := low_TL]
    setorder(res_sample_tab,sample)
    setorder(sample_tab,sample)
    sample_tab[type == "call",high_TL := res_sample_tab$best_TL + 1/2^(i+2)]
    sample_tab[type == "call",low_TL := res_sample_tab$best_TL - 1/2^(i+2)]
    res_sample_tab_list[[i]] <- res_sample_tab
    tictoc::toc()
  }
  
  
  save(res_sample_tab_list,file ="res_sample_tab_list.Rdata")
  # load("res_sample_tab_list.Rdata")
  # final_TL_estimates <- rbind(final_TL_estimates,data.table(sample = "HC",TL = 1))
  
  final_TL_estimates <- tail(res_sample_tab_list,1)[[1]][,.(sample,TL = best_TL)]
  names(res_sample_tab_list) <- iterations
  TL_estimates_table <- rbindlist(res_sample_tab_list,use.names = T,idcol = "iteration")
  
  return(list(final_estimates,TL_estimates_table))
}

test_TL_from_N_longest_CNVs <- function(final_estimates,cn_categories_tab,N_longest_CNVs = 10,min_CNV_len = mean(final_estimates[,.(end - start)]$V1) * 3){
  CNV_tab <- final_estimates[cn_pred_norm != 3]
  CNV_tab <- CNV_tab[,.(start = min(start),end = max(end),cov = median(cov),cn_pred_id = cn_pred_norm[1]),by = .(sample,chr,cn_id)]
  CNV_tab[,cnv_length := end - start]
  # CNV_tab <- CNV_tab[cnv_length > min_CNV_len]
  setorder(CNV_tab,-cnv_length)
  
  CNV_tab <- CNV_tab[,.SD[1:min(N_longest_CNVs,nrow(.SD))],by = sample]
  CNV_tab <- merge(CNV_tab,final_estimates[,.(norm_cov = median(cov / cn_categories_tab$cov_norm_factor[cn_pred_norm],na.rm = T)),by = sample],by = "sample")
  CNV_tab[,CNV_pred_TL := (cov / norm_cov - 1) / (cn_categories_tab$cn_count[cn_pred_id] / 2 - 1)]
  CNV_tab <- CNV_tab[!is.na(CNV_pred_TL)]
  TL_estimates_tab <- CNV_tab[,.(pred_TL = weighted.mean(x = CNV_pred_TL,w = cnv_length),pred_TL_sd = sqrt(weighted.mean( (CNV_pred_TL - weighted.mean(x = CNV_pred_TL,w = cnv_length))^2, cnv_length )),CNV_size_Mbp = sum(cnv_length) / 10^6),by =sample]
  TL_estimates_tab[,rel_CNV_size := round(CNV_size_Mbp / (sum(final_estimates[sample == sample[1],max(end),by = chr]$V1) / 10 ^ 6) * 100,2) ]
  TL_estimates_tab[pred_TL > 0.99,pred_TL := 0.99]
  return(TL_estimates_tab)
}

predict_CNVs <- function(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list,categories_default_tabs,initial_TL = 0.5,iterations = 1){
  
  CNV_pred_list <- list()
  TL_estimates_tab_list <- list()
  if(any(names(sample_tab) == "TL")){
    TL_estimates_tab_list[[1]] <- data.table(sample = sample_tab$sample,
                                             pred_TL = sample_tab$TL,
                                             pred_TL_sd = 0,
                                             CNV_size_kb = 0,
                                             rel_CNV_size = 0)
    sample_tab[,TL := NULL]
  } else {
    TL_estimates_tab_list[[1]] <- data.table(sample = sample_tab$sample,
                                             pred_TL = initial_TL,
                                             pred_TL_sd = 0,
                                             CNV_size_kb = 0,
                                             rel_CNV_size = 0)
  }
 
  
  for(i in seq_along(vector(length = iterations))){
    # i <-1
    print(TL_estimates_tab_list[[i]][,.(sample,TL = pred_TL)])
    print(sample_tab)
    
    tictoc::tic()
    print(paste0("iter: ",i))
    process_sample_tab <- merge.data.table(sample_tab,TL_estimates_tab_list[[i]][,.(sample,TL = pred_TL)],by = "sample",all.x = T)
    print(process_sample_tab)
    process_sample_tab[is.na(TL),TL := 0.1]
    
    res <- predict_CNV_model(process_sample_tab,trans_mat_list,cov_tab,snp_tab,library_type,categories_default_tabs)
    
    # save(res,file = "test_prediction.Rdata")
    # load("test_prediction.Rdata")
    TL_estimates_table <- res[[1]]
    final_estimates <- res[[2]]
    
    final_estimates[,cn_pred_norm := cn_pred]
    #in WGS low coverage the LOH prediction based only on snps is off
    final_estimates[cn_pred == 7,cn_pred_norm := 3]
    final_estimates <- merge.data.table(cov_tab[,.(sample,region_id,chr,start,end,cov)],final_estimates,by = c("sample","region_id"))
    
    rle_res <- final_estimates[,rle(cn_pred_norm),by = .(sample,chr)]
    rle_res[,cn_id := seq_along(values),by = .(sample,chr)]
    final_estimates[,cn_id := rep(rle_res$cn_id,rle_res$lengths)]
    TL_estimates_tab_list[[i + 1]] <- test_TL_from_N_longest_CNVs(final_estimates,categories_default_tabs$cn_categories_tab,3)
    TL_estimates_tab_list[[i + 1]] <- merge.data.table(TL_estimates_tab_list[[i]][,.(sample)],TL_estimates_tab_list[[i + 1]],by = "sample",all.x = T)
    TL_estimates_tab_list[[i + 1]][is.na(CNV_size_Mbp),CNV_size_Mbp := 0]
    TL_estimates_tab_list[[i + 1]][is.na(rel_CNV_size),rel_CNV_size := 0]
    
    
    CNV_pred_list[[i]] <- final_estimates
    
    tictoc::toc()
  }
  
  return(list(CNV_pred_list,TL_estimates_tab_list))
  
}


###########################
###########################
######               ######
######  postprocess  ######
######               ######
###########################
###########################



run_all <- function(args){  
  out_filename <- args[1]
  panel_intervals_filename <- args[2]  
  panel_snps_filename <- args[3] #filename or "no_use_snps"
  calling_type <- args[4] #tumor_only, tumor_normal, germline
  library_type <- args[5] #wgs, panel
  GC_normalization_file <- args[6] #filename or "no_GC_norm"
  cytoband_file <- args[7] #filename or "no_cytoband"
  prior_est_tumor_ratio <- as.logical(args[8])
  cov_tab_filenames <- args[(which(args == "cov") + 1):length(args)] 
  
  
  #create sample table
  if(calling_type == "tumor_normal"){
    normal_cov_tab_filenames <- args[(which(args == "normal_cov") + 1):(which(args == "cov") - 1)]
    sample_tab <- data.table(cov_tab_filenames = c(cov_tab_filenames,normal_cov_tab_filenames),
                             type = c(rep("call",length(cov_tab_filenames)),rep("normal",length(normal_cov_tab_filenames))))
  } else {
    sample_tab <- data.table(cov_tab_filenames = cov_tab_filenames,
                             type = c(rep("call",length(cov_tab_filenames))))
  }
  sample_tab[,sample := gsub(sample_regex,"\\1",cov_tab_filenames)] 
  
  if(panel_snps_filename != "no_use_snps"){
    sample_tab[,snp_tab_filenames := gsub(".region_coverage.tsv",".snpAF.tsv",cov_tab_filenames)]  
  }
  setcolorder(sample_tab,c("sample","type"))
  
  #load all data
  res <- load_and_prefilter_sample_data(sample_tab,
                                        panel_intervals_filename,
                                        panel_snps_filename,
                                        library_type,
                                        GC_normalization_file,
                                        cytoband_file)
  cov_tab <- res[[1]]
  snp_tab <- res[[2]]
  sample_tab[,c("cov_tab_filenames","snp_tab_filenames") := NULL]
  
  #prepare default transition matrices for selected cn
  categories_default_tabs <- create_copy_number_categories_default_tabs(predict_LOH = !is.null(snp_tab))
  transition_matrix <- create_transition_matrix(categories_default_tabs$cn_categories_tab)
  trans_mat_list <- prepare_transition_matrix_list(transition_matrix,categories_default_tabs$cn_categories_tab,cov_tab,dist_transition_treshold,calling_type)
  
  if(calling_type != "germline" & library_type == "wgs" & prior_est_tumor_ratio == T){
    est_ratio_score_tab <- estimate_tumor_ratio_from_hist(cov_tab)
    est_tumor_ratio_tab <- est_ratio_score_tab[,.(TL = TL[1]),by = sample]
    sample_tab <- merge(sample_tab,est_tumor_ratio_tab,by = "sample")
  }

  save(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list,categories_default_tabs,file = "test_run_objects_260123.Rdata")
  # load("test_run_objects_260123.Rdata")
  # iterations = 1
  
  res <- predict_CNVs(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list,categories_default_tabs,initial_TL = 0.99,iterations = 3)
  save(res,file = "test_CNV_call_all_samples_260123.Rdata")
  # load("test_optim_1_init_0_05.Rdata")
  # TL_estimates_tab_list <- res[[2]]
  # final_estimates <- res[[1]]
  # 
  # final_estimates[,cn_pred_vals := cn_categories_tab$cn_category[cn_pred_norm]]
  
  # per_sample_cov_means <- final_estimates[,.(mean_cov = mean(cov / cn_categories_tab$cov_norm_factor[cn_pred],na.rm = T)),by = sample]
  
  
  # load("test_optim_1_init_0_05.Rdata")
  # TL_est <- rbindlist(TL_estimates_tab_list,use.names = T,idcol = "iter")
  # TL_est[,min_iter := which.min(pred_TL),by = sample]
  # TL_est <- TL_est[,.SD[min_iter[1]],by = sample]
  # TL_est[,iter := min_iter - 1]
  # TL_est[iter == 0,iter := 1]
  # 
  # final_estimates
  # CNV_pred_list <- rbindlist(CNV_pred_list,use.names = T,idcol = "iter")
  # CNV_pred_list <- merge(CNV_pred_list,TL_est[,.(sample,iter)],by = c("sample","iter"))
  # CNV_vars <- CNV_pred_list[,.(chr = chr[1],start = min(start),end = max(end),cov = mean(cov),cn_pred_norm = cn_pred_norm[1]),by = .(sample,cn_id)]
  # 
  # sample_tab[,TL := TL_est$pred_TL]
  # res <- predict_CNV_model(sample_tab,trans_mat_list,cov_tab,NULL,library_type)
  # final_estimates <- res[[2]]
  # 
  # final_estimates[,cn_pred_norm := cn_pred]
  # #in WGS low coverage the LOH prediction based only on snps is off
  # final_estimates[cn_pred == 7,cn_pred_norm := 3]
  # final_estimates <- merge.data.table(cov_tab[,.(sample,region_id,chr,start,end,cov)],final_estimates,by = c("sample","region_id"))
  # rle_res <- final_estimates[,rle(cn_pred_norm),by = .(sample,chr)]
  # rle_res[,cn_id := seq_along(values),by = .(sample,chr)]
  # final_estimates[,cn_id := rep(rle_res$cn_id,rle_res$lengths)]
  # CNV_vars <- final_estimates[,.(start = min(start),end = max(end),cov = mean(cov),cn_pred_norm = cn_pred_norm[1]),by = .(sample,chr,cn_id)]
  
  
  
  
  
  # # res <- estimate_per_sample_tumor_load(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list)
  # # final_TL_estimates <- res[[1]]
  # # TL_estimates_table <- res[[2]]
  # 
  # load("res_sample_tab_list.Rdata")
  # 
  # 
  # 
  # 
  # TL_estimates_tab <- test_TL_from_N_longest_CNVs(final_estimates,20)
  # 
  # test <- TL_estimates_tab[rel_length > 1 & sample != "HC"]
  # 
  # 
  # final_cn_tab <-  final_estimates[,.(start = min(start),end = max(end),cov = mean(cov),cn_pred_id = cn_pred_norm[1]),by = .(sample,chr,cn_id)]
  # final_cn_tab[,cn_pred_vals := factor(cn_categories_vec[cn_pred_id],levels = cn_categories_vec)]
  # 
  # final_cn_tab[,cn_pred_count := structure(cn_categories_tab$cn_count,names = cn_categories_tab$cn_category)[cn_pred_vals]] 
  
  
}



# develop and test
input_var_call_dirs <- list.files("input_files/CNV/variant_calls/",full.names = T)
input_var_call_dirs <- setdiff(input_var_call_dirs,"input_files/CNV/variant_calls//all_samples")
cov_files <- paste0(input_var_call_dirs,"/jabCoNtool/normal.region_coverage.tsv")
snp_files <- paste0(input_var_call_dirs,"/jabCoNtool/normal.snpAF.tsv")
args <- c("results/CNV_result_tab.tsv",
          "input_files/CNV/variant_calls/all_samples/binned_genome_50000.bed",
          "no_use_snps",
          "tumor_only",
          "wgs",
          "input_files/CNV/variant_calls/all_samples/GC_profile_50000.cnp",
          "GRCh38_cytoBand.tsv",
          "True",
          "cov",cov_files)


#run as Rscript

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
print("start")
timestamp()
run_all(args)
print("end")
timestamp()
