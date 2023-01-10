######################################
# wrapper for rule: cnvkit_fix_and_segment
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnv cnvkit fix and segment \n##\n")
f.close()

shell.executable("/bin/bash")

if snakemake.params.scope == "wgs" or snakemake.params.scope == "WGS":
    command = "cnvkit.py fix" + \
              " " + snakemake.input.targetcoverage + \
              " " + snakemake.input.antitargetcoverage + \
              " " + snakemake.input.cnv_reference + \
              " -o " + snakemake.output.fix + \
              " --no-edge " + \
              " >> " + log_filename + " 2>&1"
    else:
    command = "cnvkit.py fix" + \
              " " + snakemake.input.targetcoverage + \
              " " + snakemake.input.antitargetcoverage + \
              " " + snakemake.input.cnv_reference + \
              " -o " + snakemake.output.fix + \
              " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

if snakemake.params.scope == "wgs" or snakemake.params.scope == "WGS":
    command = "cnvkit.py segment" + \
              " " + snakemake.output.fix + \
              " -o " + snakemake.output.segments + \
              " -p " + snakemake.threads + \
              " -m cbs " + \
              " -t 1e-6 " + \
              " >> " + log_filename + " 2>&1"

    else:
    command = "cnvkit.py segment" + \
              " " + snakemake.output.fix + \
              " -o " + snakemake.output.segments + \
              " -p " + snakemake.threads + \
              " -m cbs " + \
              " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
