
def get_bam_input(wildcards):
    if config["tumor_normal_paired"] == True:
        return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0],
    else:
        return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample_name)[0]


rule cnvkit_prepare_region_beds:
    input:
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output:
        reference_bed="variant_calls/cnvkit_prepare_reference/reference_bed.bed",
        target="variant_calls/all_samples/cnvkit/target.bed",
        antitarget="variant_calls/all_samples/cnvkit/antitarget.bed"
    params:
        normal_bams="mapped/*.bam",
        target=expand("{lib_ROI}.target.bed",lib_ROI=config["lib_ROI"]),
        antitarget=expand("{lib_ROI}.antitarget.bed",lib_ROI=config["lib_ROI"])
    log:
        "logs/cnvkit_prepare_reference/prepare_reference_1.log",
    threads: workflow.cores
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py access {input.reference} -o {output.reference_bed} &> {log}
        cnvkit.py autobin {params.normal_bams} -t {input.regions} -g {output.reference_bed} &>> {log}    
        echo "$PWD" 
        mv $PWD/{params.target} {output.target}
        mv $PWD/{params.antitarget} {output.antitarget}
        """

rule cnvkit_get_coverage:
    input:
        bam= get_bam_input,
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        target="variant_calls/all_samples/cnvkit/target.bed",
        antitarget="variant_calls/all_samples/cnvkit/antitarget.bed"
    output:
        targetcoverage="variant_calls/{sample_name}/cnvkit/{tumor_normal}.targetcoverage.cnn",
        antitargetcoverage="variant_calls/{sample_name}/cnvkit/{tumor_normal}.antitargetcoverage.cnn",
    log:
        "logs/{sample_name}/cnvkit/{tumor_normal}_get_coverage.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py coverage -p {threads} {input.bam} {input.target} -o {output.targetcoverage} &> {log}
        cnvkit.py coverage -p {threads} {input.bam} {input.antitarget} -o {output.antitargetcoverage} &>> {log}
        """

def normal_coverage_inputs(wildcards):
    if config["tumor_normal_paired"] == True: #normal reference
        return {'normal_coverage_inputs': set(expand("variant_calls/{sample_name}/cnvkit/normal.{tag}targetcoverage.cnn",sample_name=sample_tab.loc[
            sample_tab.tumor_normal == "normal", "donor"].tolist(),tag = ["","anti"]))}
    else:
        if len(sample_tab.index) > 4: # pool of all samples in the run if >4
            return {'normal_coverage_inputs': set(expand("variant_calls/{sample_name}/cnvkit/tumor.{tag}targetcoverage.cnn",sample_name=sample_tab.sample_name.tolist(),tag = ["","anti"]))}
        else: # FlatReference of neutral copy number (
            return {'target': "variant_calls/all_samples/cnvkit/target.bed",
                'antitarget': "variant_calls/all_samples/cnvkit/antitarget.bed"}


rule cnvkit_prepare_reference:
    input:
        unpack(normal_coverage_inputs),
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output:
        reference_cnn="variant_calls/all_samples/cnvkit/normal_reference.cnn",
    params:
        scope=config["lib_ROI"]
    log:
        "logs/all_samples/cnvkit_prepare_reference.log"
    threads: workflow.cores
    conda:
        "../wrappers/cnvkit/env.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_reference.py"


rule cnvkit_fix_and_segment:
    input:
        targetcoverage="variant_calls/{sample_name}/cnvkit/tumor.targetcoverage.cnn",
        antitargetcoverage="variant_calls/{sample_name}/cnvkit/tumor.antitargetcoverage.cnn",
        cnv_reference="variant_calls/all_samples/cnvkit/normal_reference.cnn",
    output:
        fix="variant_calls/{sample_name}/cnvkit/fixed_cov.cnr",
        segments="variant_calls/{sample_name}/cnvkit/segmented_cov.cns",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        method="hybrid",
        extra="",
        scope=config["lib_ROI"]
    log:
        "logs/{sample_name}/cnvkit/cnvkit_fix_and_segment.log"
    threads: 10
    resources: mem=8
    conda:
        "../wrappers/cnvkit/env.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_fix_and_segment.py"


rule vardict:
    input:  bam = "mapped/{bam_name}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
            regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: vcf="variant_calls/{sample_name}/cnvkit/vardict_SNV_{bam_name}.vcf",
    log: "logs/{sample_name}/cnvkit/{bam_name}_vardict.log"
    threads: 5
    resources: mem=8
    params:
        AF_threshold=0.05
    conda: "../wrappers/vardict/env.yaml"
    script: "../wrappers/vardict/script.py"


def vardict_SNV_vcf_input(wildcards):
    if config["tumor_normal_paired"] == True:
        return expand("variant_calls/{sample_name}/cnvkit/vardict_SNV_{input_bam}.vcf",sample_name = wildcards.sample_name,input_bam=sample_tab.loc[(sample_tab["donor"] == wildcards.sample_name) & (sample_tab["tumor_normal"] == "normal"), "sample_name"])[0],
    else:
        return expand("variant_calls/{sample_name}/cnvkit/vardict_SNV_{sample_name}.vcf",sample_name = wildcards.sample_name)[0],


rule cnvkit_call:
    input:
        vcf = vardict_SNV_vcf_input,
        segment="variant_calls/{sample_name}/cnvkit/segmented_cov.cns",
    output:
        calls="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
    params:
        TC=0.5, #lambda wildcards: sample_tab.loc[wildcards.sample_name, 'donor'], #tumor content?
        scope=config["lib_ROI"]
    log:
        "logs/{sample_name}/cnvkit/cnvkit_call.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py call -y -m clonal {input.segment} -v {input.vcf} -o {output.calls} --purity {params.TC} {params.extra}) &> {log}"



rule cnvkit_diagram:
    input:
        cns="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample_name}/cnvkit/fixed_cov.cnr",
    output:
        pdf="variant_calls/{sample_name}/cnvkit/cnvkit_diagram.pdf",
    params:
        extra="",
    log:
        "logs/{sample_name}/cnvkit/cnvkit_diagram.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} {params.extra}) &> {log}"

rule cnvkit_scatter:
    input:
        vcf = vardict_SNV_vcf_input,
        cns="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample_name}/cnvkit/fixed_cov.cnr",
    output:
        plot="variant_calls/{sample_name}/cnvkit/cnvkit_scatter.png",
    params:
        extra="",
    log:
        "logs/{sample_name}/cnvkit/cnvkit_scatter.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py scatter {input.cnr} -s {input.cns} -v {input.vcf} -o {output.plot} {params.extra}) &> {log}"


rule cnvkit_convert_to_vcf:
    input:
        segment="variant_calls/{sample_name}/cnvkit/CNV_calls.cns",
    output:
        vcf="variant_calls/{sample_name}/cnvkit/result_SV.vcf",
    params:
        sample_name="{sample_name}",
        hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
    log:
        "logs/{sample_name}/cnvkit/convert_to_vcf.log",
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env_python.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_vcf.py"



# rule cnvkit_batch:
#     input:
#         bam=single_bam_input,
#         bai=single_bam_bai_input,
#         cnv_reference="cnv_sv/cnvkit_prepare_reference/reference_cnn.cnn",
#     output:
#         regions="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cnr",
#         segments="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cns",
#         segments_called="cnv_sv/cnvkit_batch/{donor}/{sample_name}.call.cns",
#         bins="cnv_sv/cnvkit_batch/{donor}/{sample_name}.bintest.cns",
#         target_coverage="cnv_sv/cnvkit_batch/{donor}/{sample_name}.targetcoverage.cnn",
#         antitarget_coverage="cnv_sv/cnvkit_batch/{donor}/{sample_name}.antitargetcoverage.cnn",
#     params:
#         outdir=lambda wildcards, output: os.path.dirname(output[0]),
#         method="hybrid",
#         extra="",
#     log:
#         "logs/cnvkit_batch/{donor}/{sample_name}.log",
#     threads: 8
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     shell:
#         "(cnvkit.py batch {input.bam} "
#         "-r {input.cnv_reference} "
#         "-d {params.outdir} "
#         "-m {params.method} "
#         "{params.extra}) &> {log}"

