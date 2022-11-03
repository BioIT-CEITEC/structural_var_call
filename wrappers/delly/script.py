######################################
import os
from snakemake.shell import shell
from snakemake_wrapper_utils.bcftools import get_bcftools_opts


bcftools_opts = get_bcftools_opts(snakemake, parse_ref=False, parse_memory=False)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


exclude = snakemake.input.get("exclude", "")
if exclude:
    exclude = f"-x {exclude}"


shell(
    "(OMP_NUM_THREADS={snakemake.threads} delly call"
    " -g {snakemake.input.ref}"
    " {exclude}"
    " {extra}"
    " {snakemake.input.alns} | "
    # Convert output to specified format
    "bcftools view"
    " {bcftools_opts}"
    ") {log}"
)


# command = "gatk MergeVcfs " + \
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz "+\
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.snvs.vcf.gz "+\
#         " -R "+ snakemake.input.ref +\
#         " -D " + snakemake.input.dict +\
#         " -O " + snakemake.output.vcf +\
#         " >> " + log_filename + " 2>&1"
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)


# command =  "cp " + snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# command =  "gunzip " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
