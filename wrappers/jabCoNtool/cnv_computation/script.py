#############################################################
# wrapper for rule: cnv_computation
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cnv_computation \n##\n")
f.close()

shell.executable("/bin/bash")

command = "Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/cnv_computation.R" \
                       + " " + snakemake.output.all_res_prob_tab \
                       + " " + snakemake.input.region_bed\
                       + " " + snakemake.input.snp_bed\
                       + " norm_cov " + " ".join(snakemake.input.normal_sample_cov) \
                       + " tumor_cov " + " ".join(snakemake.input.tumor_sample_cov)\
                       + " tumor_snp_AF " + " ".join(snakemake.input.tumor_snp_AF)\
                       + " 2>> " + log_filename
f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

command = "touch " + snakemake.output.all_res_prob_tab
shell(command)
# shell(command)
