######################################
# wrapper for rule: gatk_create_panel_of_normals
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: gatk_create_panel_of_normals \n##\n")
f.close()

shell.executable("/bin/bash")


callers = list()
counts_string = ""

command = "mkdir -p " + os.path.dirname(snakemake.output.hdf5PoN)
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)


#whoami ????

for input_read_count in snakemake.input.germinal_read_counts:
    if len(snakemake.input.germinal_read_counts) > 0:
        counts_string = counts_string + "-I " + input_read_count + " "

command = "gatk CreateReadCountPanelOfNormals " + \
            counts_string + \
            " -O " + snakemake.output.hdf5PoN + \
            " >> " + log_filename + " 2>&1 "


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)
