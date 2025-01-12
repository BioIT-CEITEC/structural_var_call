library(data.table)
library(ggplot2)
library(circlize)
library(karyoploteR)
library(gtools)
library(cowplot)
library(grid)
library(gridExtra)
library(openxlsx)

#set_constants
normal_cn_value <- "2"
CNV_color_tab <- data.table(CNV = c("0","1","2","3","4","5","LOH"),
                            color = c("red","coral1","grey","skyblue2","royalblue3","navyblue","yellow3"))



plot_chromosome_circo <- function(){
  plot_tab <- final_estimates[,.(sample,chr,start,end,cov,cn_id,cn_pred)]


  print_single_sample_circo_plot <- function(plot_tab,select_sample,filename_posix = "_CNV_plot"){
    sample_plot_tab <- plot_tab[sample == select_sample]
    sample_plot_tab[,sample := NULL]
    sample_plot_tab[cov > mean(plot_tab$cov) * 2,cov := mean(plot_tab$cov) * 2]
    sample_plot_tab[,mean_cov := mean(cov),by = .(chr,cn_id)]
    sample_plot_tab[,cn_id := NULL]
    sample_plot_tab[,cn_color := colors[cn_pred]]
    sample_plot_tab[,cn_pred := NULL]
    sample_plot_tab[,chr := paste0("chr",chr)]


    output_pdf_filename <- paste0(select_sample,"_",filename_posix,".pdf")
    pdf(file = output_pdf_filename,width = 15,height = 15)

    circos.clear()
    circos.par(gap.after = c(rep(1,21),5))
    circos.initializeWithIdeogram(species = "hg38",chromosome.index = paste0("chr",1:22))


    circos.genomicTrack(sample_plot_tab,track.height = 0.3,
                        panel.fun = function(region, value, ...) {
                          circos.genomicPoints(region, value, pch = 16, cex = 0.3,numeric.column = "cov", col = "grey70" , ...)
                          circos.genomicPoints(region, value, pch = 16, cex = 0.2,numeric.column = "mean_cov", col = value$cn_color , ...)
                        })
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
      if(CELL_META$sector.numeric.index == 22) { # the last sector
        circos.text(CELL_META$cell.xlim[2] + convert_x(2.5, "mm"), CELL_META$ycenter,
                    select_sample, cex = 1, facing = "clockwise")
      }
    }, bg.border = NA)

    dev.off()
  }

  for(select_sample in unique(plot_tab$sample)){
    print_single_sample_circo_plot(plot_tab,select_sample)
  }
}


direct_cov_plot <- function(jCT_tab,out_filename_prefix,library_type,tumor_cell_fraction_table){
  
  pdf(file = paste0(out_filename_prefix,"_CNV_cov.pdf"),width = 10,height = 7)
  
  for(select_sample in unique(jCT_tab$sample)){
    # select_sample <- "PC_3_01"
    plot_tab <- jCT_tab[sample == select_sample & !is.na(cov)]

    if(!is.null(tumor_cell_fraction_table)){
      title <- get_plot_tittle_from_info(select_sample,tumor_cell_fraction_table)
    } else {
      title <- select_sample
    }
    
    plot_tab[,cn_status := ifelse(cn_pred == "2","normal","variant")]
    plot_tab[,cn_pred := factor(cn_pred,levels = CNV_color_tab$CNV)]
    plot_tab <- plot_tab[mixedorder(chr),]
    plot_tab <- plot_tab[,region_id := seq_len(nrow(plot_tab))]
    
    var_annot <- unique(plot_tab[cn_status == "variant"],by = c("chr","cn_pred","cn_id"))
    
    # Define a minimal constant for y
    min_constant <- 0
    # Compute the 1st percentile
    abs_limit <- ceiling(max(min_constant,sort(abs(plot_tab$cov - mean(plot_tab$cov)))[as.integer(length(plot_tab$cov) * 0.99)]))
    
    
    y_limits <- c(mean(plot_tab$cov)-abs_limit,mean(plot_tab$cov) + abs_limit)
    chr_tab <- plot_tab[,.(start = min(region_id),end= max(region_id)),by = chr]
    chr_tab[,mid := round(start + end) / 2]
    
    p <- ggplot(plot_tab,aes(region_id, cov)) +
      geom_point(aes(col = cn_pred),shape = 20,size = 0.2) +
      scale_color_manual(values=CNV_color_tab[match(sort(unique(plot_tab$cn_pred)),CNV)]$color) +
      scale_y_continuous(limits = y_limits) +
      geom_vline(xintercept = plot_tab[,min(region_id),by = chr]$V1[-1],linetype="dashed", color = "grey75",lwd = 0.1) +
      geom_text(data = chr_tab,aes(x = mid, y =  floor(y_limits[2]),label=chr),size = 2) +
      ggtitle(title) + ylab("standard deviation distance") +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
    
    # if(any(names(var_annot) == "region_name")){
    #   p <- p + geom_text(data = var_annot,aes(x = region_id, y =  sd_dist + (sign(sd_dist) * 0.5),label=region_name),size = 2)
    # }
    
    plot(p)
  }
  
  dev.off()
  
}

cov_plot <- function(jCT_tab,out_filename_prefix,library_type,tumor_cell_fraction_table){

  pdf(file = paste0(out_filename_prefix,"_CNV_rel_cov.pdf"),width = 10,height = 7)

  for(select_sample in unique(jCT_tab$sample)){
    # select_sample <- "PC_3_01"
    plot_tab <- jCT_tab[sample == select_sample & !is.na(cov)]
    if(library_type == "wgs"){
      plot_tab[cov > nbinom_mean * 3,cov := nbinom_mean * 3 ]
      plot_tab[,sd_dist := (cov - nbinom_mean) / sqrt(nbinom_var)]
    } else {
      plot_tab[,sd_dist := (cov - norm_dist_mean) / norm_dist_sd]
    }

    if(!is.null(tumor_cell_fraction_table)){
      title <- get_plot_tittle_from_info(select_sample,tumor_cell_fraction_table)
    } else {
      title <- select_sample
    }
    
    plot_tab[,cn_status := ifelse(cn_pred == "2","normal","variant")]
    plot_tab[,cn_pred := factor(cn_pred,levels = CNV_color_tab$CNV)]
    plot_tab <- plot_tab[mixedorder(chr),]
    plot_tab <- plot_tab[,region_id := seq_len(nrow(plot_tab))]

    var_annot <- unique(plot_tab[cn_status == "variant"],by = c("chr","cn_pred","cn_id"))

    # Define a minimal constant for y
    min_constant <- 15
    # Compute the 1st percentile
    abs_limit <- ceiling(max(min_constant,sort(abs(plot_tab$sd_dist))[as.integer(length(plot_tab$sd_dist) * 0.99)]))

    
    y_limits <- c(-abs_limit,abs_limit)
    chr_tab <- plot_tab[,.(start = min(region_id),end= max(region_id)),by = chr]
    chr_tab[,mid := round(start + end) / 2]

    p <- ggplot(plot_tab,aes(region_id, sd_dist)) +
      geom_point(aes(col = cn_pred),shape = 20,size = 0.2) +
      scale_color_manual(values=CNV_color_tab[match(sort(unique(plot_tab$cn_pred)),CNV)]$color) +
      scale_y_continuous(limits = y_limits, breaks = seq(ceiling(y_limits[1]),floor(y_limits[2]),1)) +
      geom_vline(xintercept = plot_tab[,min(region_id),by = chr]$V1[-1],linetype="dashed", color = "grey75",lwd = 0.1) +
      geom_text(data = chr_tab,aes(x = mid, y =  floor(y_limits[2]),label=chr),size = 2) +
      ggtitle(title) + ylab("standard deviation distance") +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())

    # if(any(names(var_annot) == "region_name")){
    #   p <- p + geom_text(data = var_annot,aes(x = region_id, y =  sd_dist + (sign(sd_dist) * 0.5),label=region_name),size = 2)
    # }

    plot(p)
  }

  dev.off()

}

prob_plot <- function(jCT_tab,out_filename_prefix,tumor_cell_fraction_table){
  pdf(file = paste0(out_filename_prefix,"_CNV_prob.pdf"),width = 10,height = 7)

  predicted_CNV_cols <- names(jCT_tab)[which(names(jCT_tab) %in% c(as.character(0:20),"LOH"))]

  for(select_sample in unique(jCT_tab$sample)){
    # select_sample <- "PC_3_01"
    plot_tab <- jCT_tab[sample == select_sample]

    # var_annot <- unique(plot_tab[cn_status == "variant"],by = c("chr","cn_pred","cn_id"))
    # if(!any(names(var_annot) == "region_name")){
    #   var_annot[,region_name := paste(chr,start,end,sep = "_")]
    # }

    if(!is.null(tumor_cell_fraction_table)){
      title <- get_plot_tittle_from_info(select_sample,tumor_cell_fraction_table)
    } else {
      title <- select_sample
    }
    
    # y_limits <- c(-max(abs(plot_tab$sd_dist)) - 1,max(abs(plot_tab$sd_dist)) + 1)
    chr_tab <- plot_tab[,.(start = min(region_id),end= max(region_id)),by = chr]
    chr_tab[,mid := round(start + end) / 2]

    plot_tab <- plot_tab[,c("region_id","cn_pred",predicted_CNV_cols),with =F]
    plot_tab <- melt.data.table(plot_tab,id.vars = c("region_id","cn_pred"),variable.name = "CNV",value.name = "prob")
    plot_tab[,CNV := as.factor(CNV)]
    plot_tab[,prob := exp(-prob)]
    plot_tab <- plot_tab[,.(prob = mean(prob)),by = c("region_id","cn_pred","CNV")]
    plot_tab[,prob := prob / sum(prob),by = c("region_id","cn_pred")]
    plot_tab[,predicted_CNV := F]
    plot_tab[CNV == cn_pred,predicted_CNV := T]



    p <- ggplot(plot_tab,aes(region_id, prob,fill = CNV)) +
      # geom_point(aes(col=CNV,shape=predicted_CNV),size = 0.75)+
      geom_bar(stat = "identity",colour = NA,width = 1)+
      scale_fill_manual(values=CNV_color_tab[match(levels(plot_tab$CNV),CNV)]$color) +
      # scale_shape_manual(values = c(20, 4)) +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,length.out = 11)) +
      # geom_vline(data = chr_tab,aes(xintercept = mid),linetype="dashed", color = "grey75",lwd = 0.1) +
      # # geom_text(data = chr_tab,aes(x = mid, y =  floor(y_limits[2]),label=chr),size = 2) +
      ylab("Probability") +
      xlab("Region ID") +
      ggtitle(title) +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())

    plot(p)
  }

  dev.off()

}


#TODO add reference selection
plot_chromosome_lines <- function(jCT_tab_CNVs,out_filename_prefix,reference = "GRCh38-p10"){

  karyotype <- NULL
  if(grepl("GRCh38",reference)){
    karyotype <- "hg38"
  }
  if(grepl("GRCh37",reference)){
    karyotype <- "hg19"
  }


  pdf(file = paste0(out_filename_prefix,"_chromosome_CNVs.pdf"),width = 10,height = 7)

  for(select_sample in unique(jCT_tab_CNVs$sample)){
    # select_sample <- "PC_3_01"
    plot_tab <- jCT_tab_CNVs[sample == select_sample]
    if(karyotype == "hg38" | karyotype == "hg19"){
      plot_tab[,chr := paste0("chr",chr)]
    }

    kp <- plotKaryotype(genome="hg38",chromosomes = paste0("chr",1:22))
    for(selected_cn_pred in unique(plot_tab$cn_pred)){
      kpRect(kp, data = toGRanges(plot_tab[cn_pred == selected_cn_pred,.(chr,start,end)]), y0=0, y1=0.8,col=CNV_color_tab[CNV == selected_cn_pred]$color,border=CNV_color_tab[CNV == selected_cn_pred]$color)
    }
    kpAddMainTitle(kp, main=select_sample)

    my_hist <- ggplot(CNV_color_tab, aes(CNV, fill = CNV)) + geom_bar() + scale_fill_manual(values=CNV_color_tab$color)

    # Using the cowplot package
    legend <- cowplot::get_legend(my_hist)

    legend$vp$x <- unit(.75, 'npc')
    legend$vp$y <- unit(.33, 'npc')
    grid.draw(legend)
  }


  dev.off()


}


get_jCT_CNV_tab <- function(jCT_tab,library_type){

  if(library_type == "wgs"){
    jCT_tab[,sd_dist := (cov - nbinom_mean) / sqrt(nbinom_var)]
  } else {
    jCT_tab[,sd_dist := (cov - norm_dist_mean) / norm_dist_sd]
  }

  if(any(names(jCT_tab) == "pop_HET_probability")){
    jCT_tab_CNVs <-  jCT_tab[,.(cn_pred = cn_pred[1],
                                start = min(start),
                                end = max(end),
                                len = max(end) - min(start),
                                cov = median(cov,na.rm = T),
                                cov_rel_dist = median(sd_dist,na.rm = T),
                                SNP_BAF = weighted.mean(alt_count/ref_count,pop_HET_probability),
                                SNP_weight = sum(pop_HET_probability),
                                region_names = paste(unique(region_name),collapse = ",")),by = .(sample,chr,cn_id)]
  } else {
    jCT_tab_CNVs <-  jCT_tab[,.(cn_pred = cn_pred[1],
                                start = min(start),
                                end = max(end),
                                len = max(end) - min(start),
                                cov = median(cov,na.rm = T),
                                cov_rel_dist = median(sd_dist,na.rm = T),
                                region_names = paste(unique(region_name),collapse = ",")),by = .(sample,chr,cn_id)]
  }

  if(length(unique(jCT_tab$sample)) == 1){
    jCT_tab_CNVs[,sample := NULL]
  }

  jCT_tab_CNVs[,cn_id := NULL]
  jCT_tab_CNVs <- jCT_tab_CNVs[cn_pred != normal_cn_value,]

  return(jCT_tab_CNVs)

}




get_plot_tittle_from_info <- function(select_sample,tumor_cell_fraction_table){
  title <- paste0("sample: ",select_sample,
                  "        filter: ", 
                  tumor_cell_fraction_table[sample == select_sample]$filter,
                  "        pred ctDNA: ", 
                  tumor_cell_fraction_table[sample == select_sample]$pred_TF_perc,
                  "%       perm_test_pvalue: ",
                  tumor_cell_fraction_table[sample == select_sample]$random_perm_test_p_value)
  return(title)
}



run_all <- function(args){
  final_CNV_vars_outfilename <- args[1]
  result_dir <- dirname(final_CNV_vars_outfilename)
  per_sample_results_dir <- paste0(result_dir,"/per_sample_results")
  dir.create(per_sample_results_dir,recursive = T,showWarnings = F)

  minimal_CNV_overlap <- as.numeric(args[2])
  region_bedfile <- args[3]

  #load panel of intervals
  panel_intervals <- fread(region_bedfile)
  if(length(panel_intervals) > 4){
    panel_intervals <- panel_intervals[,1:4,with = F]
  }
  if(length(panel_intervals) == 3){
    setnames(panel_intervals,c("chr","start","end"))
    panel_intervals[,region_name := paste(chr,start,sep = "_")]
  }
  setnames(panel_intervals,c("chr","start","end","region_name"))
  setkey(panel_intervals,chr,start,end)

  library_type <- args[4]

  # for now process and plot jabcontools
  if(any(args == "jabCoNtool")){
    jCT_tab <- fread(args[which(args == "jabCoNtool") + 1])
    jCT_tab[,chr := as.character(chr)]
    jCT_tab <- merge(jCT_tab,panel_intervals,by = c("chr","start","end"))
    # rle_res <- jCT_tab[,rle(cn_pred),by = .(sample,chr)]
    # rle_res[,cn_id := seq_along(values),by = .(sample,chr)]
    # jCT_tab[,cn_id := rep(rle_res$cn_id,rle_res$lengths)]
    jCT_tab_CNVs <- get_jCT_CNV_tab(jCT_tab,library_type)
    
    
    # TL_tab <- fread("structural_varcalls/all_samples/jabCoNtool/tc_fraction_prediction.tsv")
    # jCT_tab[]
    
    if(file.exists(paste0(dirname(args[which(args == "jabCoNtool") + 1]),"/tc_fraction_prediction.tsv"))){
      tumor_cell_fraction_table <- fread(paste0(dirname(args[which(args == "jabCoNtool") + 1]),"/tc_fraction_prediction.tsv"))
      # tumor_cell_fraction_table[,corrected_neg_log10_p_value := corrected_neg_log10_p_value - 3]
      # tumor_cell_fraction_table[corrected_neg_log10_p_value < 0,corrected_neg_log10_p_value := 0]
      write.xlsx(tumor_cell_fraction_table,paste0(result_dir,"/jabCoNtool_tc_fraction_prediction.xlsx"))
    } else {
      tumor_cell_fraction_table <- NULL
    }
    
    #plot res
    # prob_plot(jCT_tab,paste0(result_dir,"/jabCoNtool_all"),tumor_cell_fraction_table)
    # cov_plot(jCT_tab,paste0(result_dir,"/jabCoNtool_all"),library_type,tumor_cell_fraction_table)
    direct_cov_plot(jCT_tab,paste0(result_dir,"/jabCoNtool_all"),library_type,tumor_cell_fraction_table)

    for(select_sample in unique(jCT_tab$sample)){
      # prob_plot(jCT_tab[sample == select_sample],paste0(per_sample_results_dir,"/jabCoNtool_",select_sample),tumor_cell_fraction_table)
      # cov_plot(jCT_tab[sample == select_sample],paste0(per_sample_results_dir,"/jabCoNtool_",select_sample),library_type,tumor_cell_fraction_table)
      direct_cov_plot(jCT_tab[sample == select_sample],paste0(per_sample_results_dir,"/jabCoNtool_",select_sample),library_type,tumor_cell_fraction_table)
    }



    write.xlsx(jCT_tab_CNVs,paste0(result_dir,"/jabCoNtool_all_CNV_tab.xlsx"))
    plot_chromosome_lines(jCT_tab_CNVs,paste0(result_dir,"/jabCoNtool_all"))
    for(select_sample in unique(jCT_tab$sample)){
      plot_chromosome_lines(jCT_tab_CNVs[sample == select_sample],paste0(per_sample_results_dir,"/jabCoNtool_",select_sample))
      write.xlsx(jCT_tab_CNVs[sample == select_sample],paste0(per_sample_results_dir,"/jabCoNtool_",select_sample,"_CNV_tab.xlsx"))
    }

  }

  fwrite(jCT_tab_CNVs,file = final_CNV_vars_outfilename,sep = "\t")

}


# develop and test

# script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(paste0(script_dir,"/../.."))
# args <- readLines(con = "logs/process_and_format_CNV.log_Rargs")
# args <- strsplit(args,split = " ")[[1]]

#run as Rscript

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)















#
#
#
#
#
# run_all <- function(args){
#
#   args <- c("structural_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv",F,"tumor_only","wgs")
#
#   jabcontools_res <- args[1]
#   snps_used <- args[2]
#   calling_type <- args[3] #tumor_only, tumor_normal, germline
#   library_type <- args[4] #wgs, panel
#   # panel_intervals_filename <- args[2]
#   # panel_snps_filename <- args[3] #filename or "no_use_snps"
#   # calling_type <- args[4] #tumor_only, tumor_normal, germline
#   # library_type <- args[5] #wgs, panel
#   # GC_normalization_file <- args[6] #filename or "no_GC_norm"
#   # cytoband_file <- args[7] #filename or "no_cytoband"
#   # prior_est_tumor_ratio <- as.logical(args[8])
#   # max_CNV_frequency_in_cohort <- as.numeric(args[9])
#   # cov_tab_filenames <- args[(which(args == "cov") + 1):length(args)]
#
#   jCT_tab <- fread(jabcontools_res)
#   # jCT_tab <- copy(final_cn_pred_info_table)
#
#   missed_region_tab <- unique(jCT_tab,by = c("chr","start","end"))[,.(start = head(end, -1) + 1,end = tail(start,-1),length = tail(start,-1) - head(end, -1)),by = chr]
#   missed_region_tab <- missed_region_tab[length > 1000000]
#   missed_region_tab[,length := NULL]
#
#   jCT_tab_CNVs <-  jCT_tab[,.(cn_pred = cn_pred[1],
#                               start = min(start),
#                               end = max(end),
#                               len = max(end) - min(start),
#                               cov = median(cov)),by = .(sample,chr,cn_id)]
#
#   # jCT_tab_CNVs <- jCT_tab_CNVs[cn_pred != normal_cn_value]
#   #
#   # jCT_tab_CNVs[,cn_id := NULL]
#   # jCT_tab_CNVs[,dist := c(tail(start,-1),Inf) - end,by = .(sample,chr,cn_pred)]
#   # jCT_tab_CNVs[,cn_id := dist]
#   # jCT_tab_CNVs[cn_id < 5*10^6,cn_id := 0]
#   # rle_res <- jCT_tab_CNVs[,rle(cn_id),by = .(sample,chr)]
#   # rle_res[,cn_id := seq_along(values),by = .(sample,chr)]
#   # jCT_tab_CNVs[,cn_id := rep(rle_res$cn_id,rle_res$lengths)]
#   #
#   # jCT_tab_CNVs <- jCT_tab_CNVs[,.(cn_pred = cn_pred[1],
#   #                                                   WT_cn = WT_cn[1],
#   #                                                   start = min(start),
#   #                                                   end = max(end),
#   #                                                   len = max(end) - min(start),
#   #                                                   cov = median(cov),
#   #                                                   WT_cov = median(WT_cov),
#   #                                                   in_sample_count = mean(in_sample_count)),by = .(sample,chr,cn_id)]
#
#
#
#
#   #project_specific
#   jCT_tab[cn_pred == 0 & cov > nbinom_mean,cn_pred := "5"]
#   jCT_tab <-  merge.data.table(jCT_tab,jCT_tab[sample == "PC_3_WT",.(region_id,WT_cn = cn_pred,WT_cov = cov)],by = "region_id")
#
#   average_read_count <- mean(jCT_tab[cn_pred == "2"]$cov)
#
#   diff_CNV_tab <- jCT_tab[cn_pred != WT_cn]
#   diff_CNV_tab[,in_sample_count := .N,by = region_id]
#   diff_CNV_tab <-  diff_CNV_tab[,.(cn_pred = cn_pred[1],
#                                    WT_cn = WT_cn[1],
#                                    start = min(start),
#                                    end = max(end),
#                                    len = max(end) - min(start),
#                                    cov = median(cov),
#                                    WT_cov = median(WT_cov),
#                                    in_sample_count = mean(in_sample_count)),by = .(sample,chr,cn_id)]
#
#   diff_CNV_tab_filtered <- diff_CNV_tab[in_sample_count != 3]
#   diff_CNV_tab_filtered[,rel_diff := (cov - WT_cov) / average_read_count * 2]
#   diff_CNV_tab_filtered <- diff_CNV_tab_filtered[abs(rel_diff) > 0.9]
#   diff_CNV_tab_filtered <- diff_CNV_tab_filtered[len > 2*10^5]
#
#   diff_CNV_tab_filtered[,dist := c(tail(start,-1),Inf) - end,by = .(sample,chr)]
#   diff_CNV_tab_filtered[,cn_id_tmp := dist]
#   diff_CNV_tab_filtered[cn_id_tmp < 5*10^6,cn_id_tmp := 0]
#   rle_res <- diff_CNV_tab_filtered[,rle(cn_id_tmp),by = .(sample,chr)]
#   rle_res[,cn_id := seq_along(values),by = .(sample,chr)]
#   diff_CNV_tab_filtered[,cn_id := rep(rle_res$cn_id,rle_res$lengths)]
#
#   diff_CNV_tab_filtered <- diff_CNV_tab_filtered[,.(cn_pred = cn_pred[1],
#                                                     WT_cn = WT_cn[1],
#                                                     start = min(start),
#                                                     end = max(end),
#                                                     len = max(end) - min(start),
#                                                     cov = median(cov),
#                                                     WT_cov = median(WT_cov),
#                                                     in_sample_count = mean(in_sample_count)),by = .(sample,chr,cn_id)]
#
#   diff_CNV_tab_filtered[,cn_id := NULL]
#   diff_CNV_tab_filtered[,in_sample_count := NULL]
#   setorder(diff_CNV_tab_filtered,sample,chr,start)
#   setcolorder(diff_CNV_tab_filtered,c("sample","chr","start","end","len"))
#   write.xlsx(list(CNV_change = diff_CNV_tab_filtered),file = "CNVs_change_tab.xlsx")
#
#
#   pdf(file = "chromosome_CNVs_change.pdf",width = 10,height = 7)
#
#   karyotype = "hg38"
#
#   for(select_sample in unique(diff_CNV_tab_filtered$sample)){
#     # select_sample <- "PC_3_01"
#     plot_tab <- diff_CNV_tab_filtered[sample == select_sample & !is.na(cov)]
#     if(karyotype == "hg38" | karyotype == "hg19"){
#       plot_tab[,chr := paste0("chr",chr)]
#     }
#
#     kp <- plotKaryotype(genome="hg38",chromosomes = paste0("chr",1:22))
#     for(selected_cn_pred in unique(plot_tab[cn_pred != normal_cn_value]$cn_pred)){
#       kpRect(kp, data = toGRanges(plot_tab[cn_pred == selected_cn_pred,.(chr,start,end)]), y0=0, y1=0.4,col=CNV_color_tab[CNV == selected_cn_pred]$color)
#     }
#     for(selected_cn_pred in unique(plot_tab$WT_cn)){
#       kpRect(kp, data = toGRanges(plot_tab[WT_cn == selected_cn_pred,.(chr,start,end)]), y0=0.41, y1=0.8,col=CNV_color_tab[CNV == selected_cn_pred]$color)
#     }
#     kpAddMainTitle(kp, main=select_sample)
#
#     my_hist <- ggplot(CNV_color_tab, aes(CNV, fill = CNV)) + geom_bar() + scale_fill_manual(values=CNV_color_tab$color)
#
#     # Using the cowplot package
#     legend <- cowplot::get_legend(my_hist)
#
#     legend$vp$x <- unit(.75, 'npc')
#     legend$vp$y <- unit(.33, 'npc')
#     grid.draw(legend)
#   }
#
#   dev.off()
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# ############
# ###TODO
# ###########
#
# # CNV_tab <- data.table(caller = character(),chr = character(),start = integer(),end = integer(),cn = character(),weight = numeric(),log2 = numeric())
# #
# # #cnvkit result loading
# # if("cnvkit" %in% callers){
# #   cnv_kit_var_tab <- fread(CNV_files[which(callers == "cnvkit")])
# #   cnv_kit_var_tab[,chromosome := as.character(chromosome)]
# #   setnames(cnv_kit_var_tab,"chromosome","chr")
# #   cnv_kit_var_tab[,caller := "cnvkit"]
# #   cnv_kit_var_tab <- cnv_kit_var_tab[cn >= 0,]
# #   cnv_kit_var_tab <- cnv_kit_var_tab[cn != 2,]
# #   cnv_kit_var_tab <- cnv_kit_var_tab[cn < 10 | probes > 25,]
# #
# #   CNV_tab <- cnv_kit_var_tab[,.(caller,chr,start,end,cn,weight,log2)]
# #
# # }
# #
# # #jabCoNtool result loading
# # if("jabCoNtool" %in% callers){
# #   jabCoNtool_var_tab <- fread(CNV_files[which(callers == "jabCoNtool")])
# #   jabCoNtool_var_tab[,chr := as.character(chr)]
# #   jabCoNtool_var_tab[,cn_pred := as.character(cn_pred)]
# #   jabCoNtool_var_tab <- jabCoNtool_var_tab[sample == sample_name]
# #   jabCoNtool_var_tab[,select_nlog_prob := get(.BY[[3]]),by = .(chr,cn_id,cn_pred)]
# #   jabCoNtool_var_tab <- jabCoNtool_var_tab[,.(start = start[1],
# #                                               end = tail(end,1),
# #                                               cn = cn_pred[1],
# #                                               weight = mean(select_nlog_prob),
# #                                               log2 = mean(log2(cov / norm_dist_mean),na.rm = T),
# #                                               norm_weight = mean(get("2"))),by = .(chr,cn_id)]
# #   #cnv sanity check
# #   jabCoNtool_var_tab[norm_weight <= weight,cn := "2"]
# #   jabCoNtool_var_tab[,norm_weight := "2"]
# #   jabCoNtool_var_tab[,caller := "jabCoNtool"]
# #   jabCoNtool_var_tab <- jabCoNtool_var_tab[,.(caller,chr,start,end,cn,weight,log2)]
# #   jabCoNtool_var_tab <- jabCoNtool_var_tab[cn != 2,]
# #
# #   CNV_tab <- rbind(CNV_tab,jabCoNtool_var_tab)
# #
# # }
# #
# # setkey(CNV_tab,chr,start,end)
# # CNV_tab <- foverlaps(panel_intervals,CNV_tab,nomatch = NULL)
# # if(any(names(panel_intervals) == "region_name")){
# #   CNV_tab <- CNV_tab[,.(regions = paste(unique(region_name),collapse = ",")),by = .(caller,chr,start,end,cn,weight,log2)]
# # } else {
# #   CNV_tab <- unique(CNV_tab[,.(caller,chr,start,end,cn,weight,log2)])
# # }
# #
# # CNV_tab <- CNV_tab[,.(chr,start,end,cn,callers = caller,weight,log2,regions)]
# #
# # fwrite(CNV_tab,out_filename_prefix,sep = "\t")
#
#
#
#
#
#
#
