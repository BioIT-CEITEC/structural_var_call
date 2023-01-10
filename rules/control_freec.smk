
def bam_inputs(wildcards):
    if config["tumor_normal_paired"] == True:
        return {'tumor': expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0],
                'normal':expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0]}
    else:
        return {'tumor': expand("mapped/{tumor_bam}.bam",tumor_bam=wildcards.sample_name)[0]}


rule control_freec_config:
    input:
        chrLenfile = expand("{ref_dir}/other/control_freec/{lib_ROI}.chrLenFile",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0], #cause only chromosomes in bed can be present in fa.fai file
        GC_profile_file = expand("variant_calls/all_samples/GC_profile_{window_size}.cnp",window_size=config["wgs_bin_size"])[0],
        capture_regions = expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        snp_bed = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output:
        config = "variant_calls/{sample_name}/control_freec/config.txt"
    log: "logs/{sample_name}/control_freec/control_freec_config.log"
    threads: 10 #will be written into config file! MUST be same as control_freec rule
    resources: mem=1
    params: folder = "variant_calls/{sample_name}/control_freec",
            library_scope = config["lib_ROI"], #?
            window_size= config["wgs_bin_size"],
            calling_type=config["tumor_normal_paired"],
    conda:  "../wrappers/control_freec/env.yaml"
    script: "../wrappers/control_freec/control_freec_config.py"


rule control_freec:
    input:
        unpack(bam_inputs),
        config= "variant_calls/{sample_name}/control_freec/config.txt",
    output:
        # info = "variant_calls/{sample_name}/control_freec/info.txt",
        CNVs = "variant_calls/{sample_name}/control_freec/{sample_name}.bam_CNVs",
        CNVs_vcf = "variant_calls/{sample_name}/control_freec/result_SV.vcf",
    params:
        calling_type=config["tumor_normal_paired"],
    log: "logs/{sample_name}/control_freec/control_freec.log"
    threads: 10
    resources: mem=6
    conda:  "../wrappers/control_freec/env.yaml"
    script: "../wrappers/control_freec/control_freec_run.py"



# def bam_inputs(wildcards):
#     if config["tumor_normal_paired"] == True:
#         return {'tumor': expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0],
#                 'normal':expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0]}
#     else:
#         return {'tumor': expand("mapped/{tumor_bam}.bam",tumor_bam=wildcards.sample_name)[0]}
#
#
# rule control_freec:
#     input:
#         unpack(bam_inputs),
#         ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#         ref_fai = expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
#         GC_profile_file = expand("variant_calls/all_samples/GC_profile_{window_size}.cnp",window_size=config["wgs_bin_size"])[0],
#         region_bed = expand("variant_calls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0],
#         snp_bed = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
#         config_template = "wrappers/control_freec/control_freec_config_template_WGS.txt"
#     output:
#         # info = "variant_calls/{sample_name}/control_freec/info.txt",
#         # CNVs = "variant_calls/{sample_name}/control_freec/CNV.bed",
#         config = "variant_calls/{sample_name}/control_freec/config.txt",
#         CNVs_vcf = "variant_calls/{sample_name}/control_freec/result_SV.vcf",
#     log: "logs/{sample_name}/control_freec/control_freec.log"
#     threads: 5
#     resources: mem=6
#     params: sample_name = "variant_calls/{sample_name}/manta",
#             library_scope = config["lib_ROI"],
#             calling_type = config["tumor_normal_paired"],
#             window_size= config["wgs_bin_size"]
#     conda:  "../wrappers/control_freec/env.yaml"
#     script: "../wrappers/control_freec/script.py"
#
