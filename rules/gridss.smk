
def get_bam_input(wildcards):
    if config["tumor_normal_paired"] == True:
        input_bam_name = sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample_name), "sample_name"].min()
        return "mapped/" + input_bam_name + ".bam"
    else:
        return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample_name)[0]

def gridss_cnv_computation_inputs(wildcards):
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
        input_dict["bam_inputs"] = set(expand("mapped/{sample_name}.bam",sample_name=sample_tab.sample_name.tolist()))
        input_dict["normal_snp_AF"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.snpAF.tsv",sample_name=sample_tab.sample_name.tolist()))

    return input_dict


rule gridss_computation:
    input:  bams = expand("mapped/{sample_name}.bam",sample_name=sample_tab.sort_values(by=['tumor_normal']).sample_name.tolist())
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            ref_fai=expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
            bwa_index = expand("{ref_dir}/index/BWA/{ref}.{suffix}", ref_dir=reference_directory, ref=config["reference"], suffix=["amb", "ann", "bwt", "pac", "sa"])
    output: vcf="variant_calls/all_samples/gridss/all_sample_variants.vcf"
    params: tumor_normal_paired = config["tumor_normal_paired"]
    log:    "logs/all_samples/gridss/gridss_computation.log",
    threads: 16
    conda:  "../wrappers/gridss/cnv_computation/env.yaml"
    script: "../wrappers/gridss/cnv_computation/script.py"


rule gridss_get_per_sample_res:
    input:  all_sample_vcf="variant_calls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    output: vcf="variant_calls/{sample_name}/gridss/result_SV.vcf",
    log:    "logs/{sample_name}/gridss/get_per_sample_res.log",
    shell:
        "touch {output.vcf}"
