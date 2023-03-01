suppressMessages(library(data.table))
suppressMessages(library(VennDiagram))



run_all <- function(args){
  out_filename <- args[1]
  minimal_CNV_overlap <- as.numeric(args[2])  
  sample_name <- args[3]  
  region_bedfile <- args[4]  
  callers <- args[(which(args == "callers") + 1):(which(args == "CNV_files") - 1)] 
  CNV_files <- args[(which(args == "CNV_files") + 1):length(args)] 
  
  panel_intervals <- fread(region_bedfile)
  if(length(panel_intervals) > 4){
    panel_intervals <- panel_intervals[,1:4,with = F]
  }
  if(length(panel_intervals) == 3){
    setnames(panel_intervals,c("chr","start","end"))
  } else {
    setnames(panel_intervals,c("chr","start","end","region_name"))
  }

  setkey(panel_intervals,chr,start,end)
  
  
  CNV_tab <- data.table(caller = character(),chr = character(),start = integer(),end = integer(),cn = character(),weight = numeric(),log2 = numeric())
  
  #cnvkit result loading
  if("cnvkit" %in% callers){
    cnv_kit_var_tab <- fread(CNV_files[which(callers == "cnvkit")])
    cnv_kit_var_tab[,chromosome := as.character(chromosome)]
    setnames(cnv_kit_var_tab,"chromosome","chr")
    cnv_kit_var_tab[,caller := "cnvkit"]
    cnv_kit_var_tab <- cnv_kit_var_tab[cn >= 0,]
    cnv_kit_var_tab <- cnv_kit_var_tab[cn != 2,]
    cnv_kit_var_tab <- cnv_kit_var_tab[cn < 10 | probes > 25,]
    
    CNV_tab <- cnv_kit_var_tab[,.(caller,chr,start,end,cn,weight,log2)]
    
  }
  
  #jabCoNtool result loading
  if("jabCoNtool" %in% callers){
    jabCoNtool_var_tab <- fread(CNV_files[which(callers == "jabCoNtool")])
    jabCoNtool_var_tab[,chr := as.character(chr)]
    jabCoNtool_var_tab[,cn_pred := as.character(cn_pred)]
    jabCoNtool_var_tab <- jabCoNtool_var_tab[sample == sample_name]
    jabCoNtool_var_tab[,select_nlog_prob := get(.BY[[3]]),by = .(chr,cn_id,cn_pred)]
    jabCoNtool_var_tab <- jabCoNtool_var_tab[,.(start = start[1],
                                                end = tail(end,1),
                                                cn = cn_pred[1],
                                                weight = mean(select_nlog_prob),
                                                log2 = mean(log2(cov / norm_dist_mean),na.rm = T),
                                                norm_weight = mean(get("2"))),by = .(chr,cn_id)]
    #cnv sanity check
    jabCoNtool_var_tab[norm_weight <= weight,cn := "2"]
    jabCoNtool_var_tab[,norm_weight := "2"]
    jabCoNtool_var_tab[,caller := "jabCoNtool"]
    jabCoNtool_var_tab <- jabCoNtool_var_tab[,.(caller,chr,start,end,cn,weight,log2)]
    jabCoNtool_var_tab <- jabCoNtool_var_tab[cn != 2,]
    
    CNV_tab <- rbind(CNV_tab,jabCoNtool_var_tab)
    
  }
  
  setkey(CNV_tab,chr,start,end)
  CNV_tab <- foverlaps(panel_intervals,CNV_tab,nomatch = NULL)
  if(any(names(panel_intervals) == "region_name")){
    CNV_tab <- CNV_tab[,.(regions = paste(unique(region_name),collapse = ",")),by = .(caller,chr,start,end,cn,weight,log2)]
  } else {
    CNV_tab <- unique(CNV_tab[,.(caller,chr,start,end,cn,weight,log2)])
  }
  
  CNV_tab <- CNV_tab[,.(chr,start,end,cn,callers = caller,weight,log2,regions)]
  
  fwrite(CNV_tab,out_filename,sep = "\t")
  
}

# develop and test
# script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(paste0(script_dir,"/../.."))
# args <- c("structural_varcalls/BR-1365/CNVs.merged.tsv","0.6","BR-1365","/mnt/data/ceitec_cfg/shared/CFBioinformatics/references_backup/homo_sapiens/GRCh37-p13/intervals/BRONCO21/BRONCO21.bed","callers","cnvkit","jabCoNtool","CNV_files","structural_varcalls/BR-1365/cnvkit/CNV_calls.cns","structural_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv")

#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)
