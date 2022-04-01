#############################################################
# wrapper for rule: cnv_computation
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: cnv_computation \n##\n")
f.close()


command = "Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/cnv_computation.R" \
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
