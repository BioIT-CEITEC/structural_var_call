def final_report_inputs(wildcards):
    input = {}
    if config["organism"] == "homo_sapiens":
        tag = "svdb_query"
    else:
        tag = "merged"

    if config["tumor_normal_paired"] == True:
        if len(used_cnv_callers) > 0:
            input['svdb'] = expand("final_SV_calls/{sample_name}.{tag}.vcf",tag = tag,sample_name=sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist())
    else:
        if len(used_cnv_callers) > 0:
            input['svdb'] = expand("final_SV_calls/{sample_name}.{tag}.vcf",tag = tag,sample_name=sample_tab.sample_name)
    return input

rule final_report:
    input:
        unpack(final_report_inputs)
    output: "cnv_sv/final_report.html"
    shell:
        "touch {output}"

# vcf_cnvkit = lambda wildcards: expand("cnv_sv/cnvkit_vcf/{sample_name}.vcf",sample_name=sample_tab.sample_name),
# segment_regions = lambda wildcards: expand("cnv_sv/cnvkit_batch/{sample_name}.cnr",sample_name=sample_tab.sample_name),
# vcf_gatk = lambda wildcards: expand("cnv_sv/gatk_cnv_vcf/{sample_name}.vcf",sample_name=sample_tab.sample_name),
# manta_som_sv_vcf = lambda wildcards: expand("cnv_sv/manta/{sample_name}/results/variants/somaticSV.vcf.gz",donor=sample_tab.sample_name),
# pindel_vcf = lambda wildcards: expand("cnv_sv/pindel/{sample_name}.vcf",donor=sample_tab.sample_name),
# svdb = lambda wildcards: expand("cnv_sv/svdb_query/{sample_name}.svdb_query.vcf",sample_name=sample_tab.sample_name)

