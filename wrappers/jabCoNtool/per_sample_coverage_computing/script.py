######################################
# wrapper for rule: per_sample_coverage_counting
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: per_sample_coverage_counting \n##\n")
f.close()

command = "bedtools coverage " + \
          " -a " + snakemake.input.region_bed + \
          " -b " + snakemake.input.bam + \
          " > " + snakemake.output.cov_tab + \
          " 2>> " + snakemake.log.run

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)
