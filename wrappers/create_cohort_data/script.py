######################################
# wrapper for rule: create_cohort_data
######################################
import os
import subprocess
from snakemake.shell import shell
import re

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: create_cohort_data \n##\n")
f.close()

shell.executable("/bin/bash")

if hasattr(snakemake.input, "previous_cohort_data"):
    command = "tar -xzf " + snakemake.input.previous_cohort_data + " -C cohort_data/ "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
else:
    command = "mkdir -p cohort_data/cohort_data/ "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

if hasattr(snakemake.input, "cnvkit_normal_coverage_inputs"):
    command = "mkdir -p cohort_data/cohort_data/cnvkit "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    f = open(log_filename, 'at')
    for filename in snakemake.input.cnvkit_normal_coverage_inputs:
        out_filename = filename.replace("structural_varcalls/","")
        out_filename = out_filename.replace("/cnvkit/", "_")
        command = "cp " + filename + " cohort_data/cohort_data/cnvkit/" + out_filename
        f.write("## COMMAND: " + command + "\n")
        shell(command)
    f.close()


if hasattr(snakemake.input, "jabCoNtool_all_res_prob_tab"):
    command = "mkdir -p cohort_data/cohort_data/jabCoNtool "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = 'awk -F"\\t" \'{{print $1 "\\t" $3 "\\t" $4 "\\t" $5 "\\t" $8 "\\t" $11 "\\t" $13 "\\t" $14}}\' ' + snakemake.input.jabCoNtool_all_res_prob_tab + " > cohort_data/cohort_data/jabCoNtool/region_info.tsv"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

command = "tar -czvf cohort_data/cohort_cnv_info.tar.gz cohort_data/cohort_data"
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "touch " + snakemake.output.update_finished_checkfile
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)




