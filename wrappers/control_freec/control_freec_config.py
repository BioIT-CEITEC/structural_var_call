######################################
# wrapper for rule: control_freec_config
######################################
import sys, os
from re import sub
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: control freec make config file \n##\n")
f.close()

shell.executable("/bin/bash")

version = sub(" : a.*", "",
              str(subprocess.Popen("freec 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8'))
f = open(log_filename, 'at')
f.write("## VERSION: " + version + "\n")
f.close()

## FUNCTIONS
if snakemake.params.library_scope == "wgs" or snakemake.params.library_scope == "WGS":
    config_template = os.path.abspath(os.path.join(os.path.dirname(__file__), "config_WGS.txt"))
else:
    config_template = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "config_exome.txt"))

if snakemake.params.library_scope == "wgs" or snakemake.params.library_scope == "WGS":
    gemFile = str("/mnt/ssd/ssd_1/workspace/ROBIN2/CNV_DATA/data/GRCh38-p10_index.gem")  # only for WGS"
else:
    gemFile = str("")

# create dir
if not os.path.exists(snakemake.params.folder):
    os.makedirs(snakemake.params.folder)

# EXOME or WGS
if snakemake.params.library_scope == "wgs":
    ######################################################################################
    # Read in the config template
    f = open(config_template, 'at')
    with open(config_template, 'r') as file:
        filedata = file.read()
    # Replace specific parameters
    filedata = filedata.replace('X_window_X', str(snakemake.params.window_size))
    filedata = filedata.replace('X_ref_fai_X', str(snakemake.input.chrLenfile))
    filedata = filedata.replace('X_threads_X', str(snakemake.threads))
    filedata = filedata.replace('X_GCcontentProfile_X', str(snakemake.input.GC_profile_file))
    filedata = filedata.replace('X_SNPfile_X', str(snakemake.input.snp_bed))
    filedata = filedata.replace('X_outputDir_X', str(snakemake.params.folder))

    if len(gemFile) > 1:
        filedata = filedata.replace('X_gemMappabilityFile_X', gemFile)
    else:
        filedata = filedata.replace('gemMappabilityFile', "#gemMappabilityFile")

    # Write the specific output
    with open(snakemake.output.config, 'w') as file:
        file.write(filedata)
    ######################################################################################
else:
    ######################################################################################
    # Read in the config template
    f = open(config_template, 'at')
    with open(config_template, 'r') as file:
        filedata = file.read()
    # Replace specific parameters
    filedata = filedata.replace('X_window_X', str(snakemake.params.window_size))
    filedata = filedata.replace('X_ref_fai_X', str(snakemake.input.chrLenfile))
    filedata = filedata.replace('X_threads_X', str(snakemake.threads))
    filedata = filedata.replace('X_GCcontentProfile_X', str(snakemake.input.GC_profile_file))
    filedata = filedata.replace('X_SNPfile_X', str(snakemake.input.snp_bed))
    filedata = filedata.replace('X_outputDir_X', str(snakemake.params.folder))
    # jen exome
    filedata = filedata.replace('X_captureRegions_X', str(snakemake.input.region_bed))

    # Write the specific output
    with open(snakemake.output.config, 'w') as file:
        file.write(filedata)
    ######################################################################################
