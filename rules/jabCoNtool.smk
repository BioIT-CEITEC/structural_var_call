
def get_bam_input(wildcards):
    if config["tumor_normal_paired"] == True:
        input_bam_name = sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample_name), "sample_name"]
        return "mapped/" + input_bam_name + ".bam"
    else:
        return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample_name)[0]

def get_region_bed_input(wildcards):
    if config["lib_ROI"] != "wgs":
        return expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    else:
        return expand("variant_calls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0],

rule jabCoNtool_per_sample_coverage:
    input:  bam= get_bam_input,
            region_bed = get_region_bed_input,
            ref_dict= expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output: cov_tab = "variant_calls/{sample_name}/jabCoNtool/{tumor_normal}.region_coverage.tsv",
    log:    "logs/{sample_name}/jabCoNtool/{tumor_normal}_get_coverage.log"
    threads: 1
    resources: mem=10
    conda:  "../wrappers/jabCoNtool/per_sample_coverage_computing/env.yaml"
    script: "../wrappers/jabCoNtool/per_sample_coverage_computing/script.py"

rule jabCoNtool_per_sample_snp_AF:
    input:  bam = get_bam_input,
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            snp_tsv = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.tsv",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: snp_tab = "variant_calls/{sample_name}/jabCoNtool/{tumor_normal}.snpAF.tsv",
    log:    "logs/{sample_name}/jabCoNtool/{tumor_normal}_get_snpAF.log"
    threads: 1
    resources: mem=10
    conda:  "../wrappers/jabCoNtool/per_sample_snp_AF_computing/env.yaml"
    script: "../wrappers/jabCoNtool/per_sample_snp_AF_computing/script.py"


def jabCoNtool_cnv_computation_inputs(wildcards):
    input_dict = {}
    if config["tumor_normal_paired"] == True:
        input_dict["normal_sample_cov"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.region_coverage.tsv",sample_name=sample_tab.loc[
            sample_tab.tumor_normal == "normal", "donor"].tolist()))
        input_dict["tumor_sample_cov"] = set(expand("variant_calls/{sample_name}/jabCoNtool/tumor.region_coverage.tsv",sample_name=
            sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()))
        input_dict["tumor_snp_AF"] = set(expand("variant_calls/{sample_name}/jabCoNtool/tumor.snpAF.tsv",sample_name=
            sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()))
        input_dict["normal_snp_AF"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.snpAF.tsv",sample_name=
            sample_tab.sample_name.tolist()))
    else:
        input_dict["normal_sample_cov"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.region_coverage.tsv",sample_name=sample_tab.sample_name.tolist()))
        input_dict["normal_snp_AF"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.snpAF.tsv",sample_name=sample_tab.sample_name.tolist()))
    if config["use_cohort_data"] == True:
        input_dict["cohort_data"] = "cohort_data/cohort_cnv_info.tar.gz"
    if config["lib_ROI"] != "wgs":
        input_dict["GC_profile_file"] = expand("variant_calls/all_samples/GC_profile_{window_size}.cnp",window_size=config["wgs_bin_size"])[0],
    return input_dict


rule jabCoNtool_cnv_computation:
    input: unpack(jabCoNtool_cnv_computation_inputs),
           region_bed = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
           snp_bed = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: all_res_prob_tab="variant_calls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    params: tumor_normal_paired = config["tumor_normal_paired"]
    log:    "logs/all_samples/jabCoNtool/cnv_computation.log",
    threads: workflow.cores
    conda:  "../wrappers/jabCoNtool/cnv_computation/env.yaml"
    script: "../wrappers/jabCoNtool/cnv_computation/script.py"


rule jabCoNtool_get_per_sample_res:
    input:  all_res_prob_tab="variant_calls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    output: vcf="variant_calls/{sample_name}/jabCoNtool/result_SV.vcf",
    log:    "logs/{sample_name}/jabCoNtool/get_per_sample_res.log",
    shell:
        "touch {output.vcf}"
    # conda:  "../wrappers/jabCoNtool/get_per_sample_res/env.yaml"
    # script: "../wrappers/jabCoNtool/get_per_sample_res/script.py"
