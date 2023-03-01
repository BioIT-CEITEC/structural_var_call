#############################################################
# wrapper for rule: merge_CNVs
#############################################################
import os
import sys
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_CNVs \n##\n")
f.close()

shell.executable("/bin/bash")

CNV_caller_list = []
CNV_call_file_list = []

if hasattr(snakemake.input, "gatk_cnv_variants"):
    CNV_caller_list.append("gatk_cnv")
    CNV_call_file_list.append(snakemake.input.gatk_cnv_variants)

if hasattr(snakemake.input, "cnvkit_variants"):
    CNV_caller_list.append("cnvkit")
    CNV_call_file_list.append(snakemake.input.cnvkit_variants)

if hasattr(snakemake.input, "jabCoNtool_variants"):
    CNV_caller_list.append("jabCoNtool")
    CNV_call_file_list.append(snakemake.input.jabCoNtool_variants)

if hasattr(snakemake.input, "control_freec_variants"):
    CNV_caller_list.append("control_freec")
    CNV_call_file_list.append(snakemake.input.control_freec_variants)

if len(CNV_caller_list) > 0:
    command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_after_merge.R" \
              + " " + snakemake.output.merged \
              + " " + str(snakemake.params.overlap) \
              + " " + str(snakemake.wildcards.sample_name) \
              + " " + str(snakemake.input.region_bed) \
              + " callers " + " ".join(CNV_caller_list) \
              + " CNV_files " + " ".join(CNV_call_file_list) \
              + " >> " + log_filename + " 2>&1"

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
    f.close()

    shell(command)
else:
    print >> sys.stderr, "No callers set for a sample."
    sys.exit(1)


