######################################
# wrapper for rule: control_freec_run
######################################
import os
from re import sub
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: control freec \n##\n")
f.close()

shell.executable("/bin/bash")

version = sub(" : a.*", "", str(subprocess.Popen("freec 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8'))
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()


if snakemake.params.calling_type == True: # tumor-normal paired
    command = "freec -conf " + snakemake.input.config + " -sample " + snakemake.input.tumor + " -control " + snakemake.input.normal + " >> " + log_filename + " 2>&1"
else:
    command = "freec -conf " + snakemake.input.config + " -sample " + snakemake.input.tumor + " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# remove workspace - lots of files messing copying speeds
shell("touch " + snakemake.output.CNVs_vcf)


