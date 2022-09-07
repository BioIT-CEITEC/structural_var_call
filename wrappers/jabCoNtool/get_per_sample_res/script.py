#############################################################
# wrapper for rule: cnv_computation
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: per_sample_results \n##\n")
f.close()

command = "Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/get_per_sample_res.R" \
                       + " " + snakemake.output.cnv_res \
                       + " " + snakemake.input.region_bed[0]\
                       + " " + snakemake.input.snp_bed[0]\
                       + " " + snakemake.params.all_kit_results[0]\
                       + " " + " ".join(snakemake.input.cov_tabs)\
                       + " snps " + " ".join(snakemake.input.snp_tabs) \
                       + " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
f.write("## args: \nc(\"" + "\",\"".join(command.split(" ")[3:]) + "\")\n")
f.close()
# shell(command)

command = "touch " + snakemake.output.all_res_prob_tab
shell(command)