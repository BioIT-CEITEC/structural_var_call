
def single_bam_input(wildcards):
    return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])

def single_bam_bai_input(wildcards):
    return expand("mapped/{input_bam}.bam.bai",input_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])

def normal_bam_inputs(wildcards):
    return set(expand("mapped/{input_bam}.bam",zip\
            ,input_bam=sample_tab.loc[sample_tab.tumor_normal == "normal" , "sample_name"].tolist()))


rule cnvkit_prepare_reference_first:
    input:
        normal_bams=normal_bam_inputs,
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output:
        reference_bed="cnv_sv/cnvkit_prepare_reference/reference_bed.bed",
        target=expand("cnv_sv/cnvkit_prepare_reference/{lib_ROI}.target.bed",lib_ROI=config["lib_ROI"]),
        antitarget=expand("cnv_sv/cnvkit_prepare_reference/{lib_ROI}.antitarget.bed",lib_ROI=config["lib_ROI"])
    params:
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
        cnvkit.py autobin {input.normal_bams} -t {input.regions} -g {output.reference_bed} &>> {log}    
        echo "$PWD" 
        mv $PWD/{params.target} {output.target}
        mv $PWD/{params.antitarget} {output.antitarget}
        """

rule cnvkit_prepare_reference_second:
    input:
        normal_bam=lambda wildcards: expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample_name) & (sample_tab["tumor_normal"] == "normal"), "sample_name"]),
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        target=expand("cnv_sv/cnvkit_prepare_reference/{lib_ROI}.target.bed",lib_ROI=config["lib_ROI"]),
        antitarget=expand("cnv_sv/cnvkit_prepare_reference/{lib_ROI}.antitarget.bed",lib_ROI=config["lib_ROI"])
    output:
        targetcoverage="cnv_sv/cnvkit_prepare_reference/{sample_name}.targetcoverage.cnn",
        antitargetcoverage="cnv_sv/cnvkit_prepare_reference/{sample_name}.antitargetcoverage.cnn",
    log:
        "logs/cnvkit_prepare_reference/prepare_reference_{sample_name}.log",
    threads: workflow.cores
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py coverage {input.normal_bam} {input.target} -o {output.targetcoverage} &> {log}
        cnvkit.py coverage {input.normal_bam} {input.antitarget} -o {output.antitargetcoverage} &>> {log}
        """

def anti_targetcoverage(wildcards):
    return    {'targetcoverage': expand("cnv_sv/cnvkit_prepare_reference/{sample_name}.targetcoverage.cnn",sample_name=sample_tab.loc[sample_tab.tumor_normal == "normal" , "sample_name"]),\
        'antitargetcoverage':expand("cnv_sv/cnvkit_prepare_reference/{sample_name}.antitargetcoverage.cnn",sample_name=sample_tab.loc[sample_tab.tumor_normal == "normal" , "sample_name"])}

rule cnvkit_prepare_reference_third:
    input:
        unpack(anti_targetcoverage),
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output:
        reference_cnn="cnv_sv/cnvkit_prepare_reference/reference_cnn.cnn",
    params:
        folder="cnv_sv/cnvkit_prepare_reference",
    log:
        "logs/cnvkit_prepare_reference/prepare_reference_3.log",
    threads: workflow.cores
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py reference {params.folder}/*.{{,anti}}targetcoverage.cnn --fasta {input.reference} -o {output.reference_cnn} &> {log}
        """
#jak to udelat jinak cnv_sv/cnvkit_prepare_reference/*.{{,anti}}targetcoverage.cnn pres snakemake, alespon cestu - {params.folder}/*.{{,anti}}targetcoverage.cnn

rule cnvkit_batch:
    input:
        bam=single_bam_input,
        bai=single_bam_bai_input,
        cnv_reference="cnv_sv/cnvkit_prepare_reference/reference_cnn.cnn",
    output:
        regions="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cnr",
        segments="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cns",
        segments_called="cnv_sv/cnvkit_batch/{donor}/{sample_name}.call.cns",
        bins="cnv_sv/cnvkit_batch/{donor}/{sample_name}.bintest.cns",
        target_coverage="cnv_sv/cnvkit_batch/{donor}/{sample_name}.targetcoverage.cnn",
        antitarget_coverage="cnv_sv/cnvkit_batch/{donor}/{sample_name}.antitargetcoverage.cnn",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        method="hybrid",
        extra="",
    log:
        "logs/cnvkit_batch/{donor}/{sample_name}.log",
    threads: 8
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py batch {input.bam} "
        "-r {input.cnv_reference} "
        "-d {params.outdir} "
        "-m {params.method} "
        "{params.extra}) &> {log}"


rule vardict:
    input:  bam = lambda wildcards: expand("mapped/{input_bam}.bam", input_bam=sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample_name) & (sample_tab["tumor_normal"] == "normal"), "sample_name"])[0],
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
            regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: vcf="cnv_sv/vardict_germline_vcf/{sample_name}.germline.vcf",
    log: "logs/vardict/{sample_name}_vardict.log"
    threads: 8
    resources:
        mem_mb=8000
    params:
        AF_threshold=config["min_variant_frequency"]
    conda: "../wrappers/vardict/env.yaml"
    script: "../wrappers/vardict/script.py"


rule cnvkit_call:
    input:
        segment=lambda wildcards:expand("cnv_sv/cnvkit_batch/{donor}/{sample_name}.cns",sample_name = sample_tab.loc[sample_tab["sample_name"] == wildcards.sample_name, "sample_name"], donor=sample_tab.loc[sample_tab["sample_name"] == wildcards.sample_name, "donor"])[0],
        vcf=lambda wildcards: expand("cnv_sv/vardict_germline_vcf/{sample_name}.germline.vcf",sample_name = sample_tab.loc[sample_tab["sample_name"] == wildcards.sample_name, "germinal"])[0],
    output:
        segment="cnv_sv/cnvkit_call/{sample_name}.loh.cns",
    params:
        TC=0.5, #lambda wildcards: sample_tab.loc[wildcards.sample_name, 'donor'], #tumor content?
        extra=""
    log:
        "logs/cnvkit_call/{sample_name}.loh.cns.log",
    threads: 8
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py call {input.segment} -v {input.vcf} -o {output.segment} --purity {params.TC} {params.extra}) &> {log}"



rule cnvkit_diagram:
    input:
        cns="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cns",
        cnr="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cnr",
    output:
        pdf="cnv_sv/cnvkit_diagram/{sample_name}.pdf",
    params:
        extra="",
    log:
        "logs/cnvkit_diagram/{sample_name}.pdf.log",
    threads: 8
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} {params.extra}) &> {log}"

rule cnvkit_scatter:
    input:
        segments="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cns",
        segment_regions="cnv_sv/cnvkit_batch/{donor}/{sample_name}.cnr",
        vcf=lambda wildcards: expand("cnv_sv/vardict_germline_vcf/{sample_name}.germline.vcf",sample_name=sample_tab.loc[sample_tab["sample_name"] == wildcards.sample_name, "germinal"])[0],
    output:
        plot="cnv_sv/cnvkit_scatter/{sample_name}.png",
    params:
        extra="",
    log:
        "logs/cnvkit_scatter/{sample_name}.log",
    threads: 8
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py scatter {input.segment_regions} -s {input.segments} -v {input.vcf} -o {output.plot} {params.extra}) &> {log}"


rule cnvkit_vcf:
    input:
        segment="cnv_sv/cnvkit_call/{sample_name}.loh.cns",
    output:
        vcf="cnv_sv/cnvkit_vcf/{sample_name}.vcf",
    params:
        sample_name="{sample_name}",
        hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
    log:
        "logs/cnvkit_vcf/{sample_name}.vcf.log",
    threads: 8
    conda:
        "../wrappers/cnvkit/env_python.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_vcf.py"
