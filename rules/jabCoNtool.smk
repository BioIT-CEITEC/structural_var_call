
def get_bam_input(wildcards):
    if config["tumor_normal_paired"] == True:
        return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0],
    else:
        return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample_name)[0]

rule per_sample_coverage_computing:
    input:  bam= get_bam_input,
            region_bed = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: cov_tab = "variant_calls/{sample_name}/jabCoNtool/{tumor_normal}.region_coverage.tsv",
    log:    "logs/{sample_name}/jabCoNtool/{tumor_normal}_get_coverage.log"
    threads: 2
    conda:  "../wrappers/jabCoNtool/per_sample_coverage_computing/env.yaml"
    script: "../wrappers/jabCoNtool/per_sample_coverage_computing/script.py"

rule per_sample_snp_AF_computing:
    input:  bam= get_bam_input,
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            snp_bed = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: snp_tab = "variant_calls/{sample_name}/jabCoNtool/{tumor_normal}.snpAF.tsv",
    log:    "logs/{sample_name}/jabCoNtool/{tumor_normal}_get_snpAF.log"
    threads: 2
    conda:  "../wrappers/jabCoNtool/per_sample_snp_AF_computing/env.yaml"
    script: "../wrappers/jabCoNtool/per_sample_snp_AF_computing/script.py"


def cnv_computation_inputs(wildcards):
    input_dict = {}
    if config["tumor_normal_paired"] == True:
        input_dict["normal_sample_cov"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.region_coverage.tsv",sample_name=sample_tab.loc[
            sample_tab.tumor_normal == "normal", "donor"].tolist()))
        input_dict["tumor_sample_cov"] = set(expand("variant_calls/{sample_name}/jabCoNtool/tumor.region_coverage.tsv",sample_name=
            sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()))
        input_dict["tumor_snp_AF"] = set(expand("variant_calls/{sample_name}/jabCoNtool/tumor.snpAF.tsv",sample_name=
        sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()))
    else:
        input_dict["normal_sample_cov"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.region_coverage.tsv",sample_name=sample_tab.sample_name.tolist()))
        input_dict["normal_snp_AF"] = set(expand("variant_calls/{sample_name}/jabCoNtool/normal.snpAF.tsv",sample_name=sample_tab.sample_name.tolist()))
    if config["use_cohort_data"] == True:
        input_dict["cohort_data"] = "cohort_data/cohort_cnv_info.tar.gz"


rule cnv_computation:
    input: cnv_computation_inputs,
           region_bed = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
           snp_bed = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: all_res_prob_tab="variant_calls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    log:    "logs/all_samples/jabCoNtool/cnv_computation.log",
    threads: 20
    conda:  "../wraps/jabCoNtool/cnv_computation/env.yaml"
    script: "../wraps/jabCoNtool/cnv_computation/script.py"


rule get_per_sample_res:
    input:  all_res_prob_tab="variant_calls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    output: vcf="variant_calls/{sample_name}/jabCoNtool/CNV.vcf",
    log:    "logs/{sample_name}/jabCoNtool/get_per_sample_res.log",
    conda:  "../wrappers/jabCoNtool/get_per_sample_res/env.yaml"
    script: "../wrappers/jabCoNtool/get_per_sample_res/script.py"

# rule cnv_visualization:
#     input:  cnv_res = ADIR+"/results/CNV_result_tab.tsv"
#     output: viz_pdf = ADIR+"/results/all_samples_cnv.pdf"
#     log:    run = ADIR+"/results/cnv_computation.log"
#     threads: 1
#     conda:  "../wraps/CNV_analysis/cnv_vizualization/env.yaml"
#     script: "../wraps/CNV_analysis/cnv_vizualization/script.py"


# def normal_coverage_inputs(wildcards):
#     if config["tumor_normal_paired"] == True:
#         return {'normal_coverage_inputs': set(expand("variant_calls/{sample_name}/cnvkit/normal.{tag}targetcoverage.cnn",sample_name=sample_tab.loc[
#             sample_tab.tumor_normal == "normal", "donor"].tolist(),tag = ["","anti"]))}
#     else:
#         if len(sample_tab.index) > 4:
#             return {'normal_coverage_inputs': set(expand("variant_calls/{sample_name}/cnvkit/tumor.{tag}targetcoverage.cnn",sample_name=sample_tab.sample_name.tolist(),tag = ["","anti"]))}
#         else:
#             return {'target': "variant_calls/all_samples/cnvkit/target.bed",
#                 'antitarget': "variant_calls/all_samples/cnvkit/antitarget.bed"}
#
#
# rule cnv_computation:
#     input:
#         bam= get_bam_input,
#         reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#         target="variant_calls/all_samples/cnvkit/target.bed",
#         antitarget="variant_calls/all_samples/cnvkit/antitarget.bed"
#     output:
#         targetcoverage="variant_calls/{sample_name}/cnvkit/{tumor_normal}.targetcoverage.cnn",
#         antitargetcoverage="variant_calls/{sample_name}/cnvkit/{tumor_normal}.antitargetcoverage.cnn",
#     log:
#         "logs/{sample_name}/cnvkit/{tumor_normal}_get_coverage.log"
#     threads: workflow.cores
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     shell:
#         """
#         cnvkit.py coverage {input.bam} {input.target} -o {output.targetcoverage} &> {log}
#         cnvkit.py coverage {input.bam} {input.antitarget} -o {output.antitargetcoverage} &>> {log}
#         """
#
#
# rule cnvkit_prepare_reference:
#     input:
#         unpack(normal_coverage_inputs),
#         reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#     output:
#         reference_cnn="variant_calls/all_samples/cnvkit/normal_reference.cnn",
#     log:
#         "logs/all_samples/cnvkit_prepare_reference.log"
#     threads: workflow.cores
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     script:
#         "../wrappers/cnvkit/cnvkit_reference.py"
#     # shell:
#     #     """
#     #     cnvkit.py reference {params.folder}/*.{{,anti}}targetcoverage.cnn --fasta {input.reference} -o {output.reference_cnn} &> {log}
#     #     """
# #jak to udelat jinak cnv_sv/cnvkit_prepare_reference/*.{{,anti}}targetcoverage.cnn pres snakemake, alespon cestu - {params.folder}/*.{{,anti}}targetcoverage.cnn
#
# rule cnvkit_fix_and_segment:
#     input:
#         targetcoverage="variant_calls/{sample_name}/cnvkit/tumor.targetcoverage.cnn",
#         antitargetcoverage="variant_calls/{sample_name}/cnvkit/tumor.antitargetcoverage.cnn",
#         cnv_reference="variant_calls/all_samples/cnvkit/normal_reference.cnn",
#     output:
#         fix="variant_calls/{sample_name}/cnvkit/fixed_cov.cnr",
#         segments="variant_calls/{sample_name}/cnvkit/segmented_cov.cns",
#     params:
#         outdir=lambda wildcards, output: os.path.dirname(output[0]),
#         method="hybrid",
#         extra="",
#     log:
#         "logs/{sample_name}/cnvkit/cnvkit_fix_and_segment.log"
#     threads: 8
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     script:
#         "../wrappers/cnvkit/cnvkit_fix_and_segment.py"
#
#
# rule vardict:
#     input:  bam = "mapped/{bam_name}.bam",
#             ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#             refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
#             regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
#     output: vcf="variant_calls/{sample_name}/cnvkit/vardict_SNV_{bam_name}.vcf",
#     log: "logs/{sample_name}/cnvkit/{bam_name}_vardict.log"
#     threads: 8
#     resources:
#         mem_mb=8000
#     params:
#         AF_threshold=config["min_variant_frequency"]
#     conda: "../wrappers/vardict/env.yaml"
#     script: "../wrappers/vardict/script.py"
#
#
# def vardict_SNV_vcf_input(wildcards):
#     if config["tumor_normal_paired"] == True:
#         return expand("variant_calls/{sample_name}/cnvkit/vardict_SNV_{bam_name}.vcf",sample_name = wildcards.sample_name,input_bam=sample_tab.loc[(sample_tab["donor"] == wildcards.sample_name) & (sample_tab["tumor_normal"] == "normal"), "sample_name"])[0],
#     else:
#         return expand("variant_calls/{sample_name}/cnvkit/vardict_SNV_{sample_name}.vcf",sample_name = wildcards.sample_name)[0],
#
#
# rule cnvkit_call:
#     input:
#         vcf = vardict_SNV_vcf_input,
#         segment="variant_calls/{sample_name}/cnvkit/segmented_cov.cns",
#     output:
#         calls="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
#     params:
#         TC=0.5, #lambda wildcards: sample_tab.loc[wildcards.sample_name, 'donor'], #tumor content?
#         extra=""
#     log:
#         "logs/{sample_name}/cnvkit/cnvkit_call.log"
#     threads: 8
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     shell:
#         "(cnvkit.py call -y -m clonal {input.segment} -v {input.vcf} -o {output.calls} --purity {params.TC} {params.extra}) &> {log}"
#
#
#
# rule cnvkit_diagram:
#     input:
#         cns="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
#         cnr="variant_calls/{sample_name}/cnvkit/fixed_cov.cnr",
#     output:
#         pdf="variant_calls/{sample_name}/cnvkit/cnvkit_diagram.pdf",
#     params:
#         extra="",
#     log:
#         "logs/{sample_name}/cnvkit/cnvkit_diagram.log"
#     threads: 8
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     shell:
#         "(cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} {params.extra}) &> {log}"
#
# rule cnvkit_scatter:
#     input:
#         vcf = vardict_SNV_vcf_input,
#         cns="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
#         cnr="variant_calls/{sample_name}/cnvkit/fixed_cov.cnr",
#     output:
#         plot="variant_calls/{sample_name}/cnvkit/cnvkit_scatter.png",
#     params:
#         extra="",
#     log:
#         "logs/{sample_name}/cnvkit/cnvkit_scatter.log"
#     threads: 8
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     shell:
#         "(cnvkit.py scatter {input.cnr} -s {input.cns} -v {input.vcf} -o {output.plot} {params.extra}) &> {log}"
#
#
# rule cnvkit_convert_to_vcf:
#     input:
#         segment="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
#     output:
#         vcf="variant_calls/{sample_name}/cnvkit/CNV.vcf",
#     params:
#         sample_name="{sample_name}",
#         hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
#         het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
#         dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
#     log:
#         "logs/{sample_name}/cnvkit/convert_to_vcf.log",
#     threads: 8
#     conda:
#         "../wrappers/cnvkit/env_python.yaml"
#     script:
#         "../wrappers/cnvkit/cnvkit_vcf.py"
#
#
