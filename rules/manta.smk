def normal_bam_input(wildcards):
    return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.donor), "sample_name"])

def normal_bam_bai_input(wildcards):
    return expand("mapped/{input_bam}.bam.bai",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.donor), "sample_name"])

def tumor_bam_input(wildcards):
    return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.donor), "sample_name"])

def tumor_bam_bai_input(wildcards):
    return expand("mapped/{input_bam}.bam.bai",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.donor), "sample_name"])

rule config_manta_tn:
    input:
        bam_t = tumor_bam_input,
        bai_t = tumor_bam_bai_input,
        bam_n = normal_bam_input,
        bai_n = normal_bam_bai_input,
        ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output:
        scrpt="cnv_sv/manta/{donor}/runWorkflow.py",
    log:
        "logs/manta/{donor}/run_workflow.py.log",
    threads: 8
    conda:
        "../wrappers/manta/env.yaml"
    shell:
        "configManta.py "
        "--tumorBam={input.bam_t} "
        "--normalBam={input.bam_n} "
        "--referenceFasta={input.ref} "
        "--runDir=cnv_sv/manta/{wildcards.donor} &> {log}"

rule manta_run_workflow_tn:
    input:
        bam_t=tumor_bam_input,
        bai_t=tumor_bam_bai_input,
        bam_n=normal_bam_input,
        bai_n=normal_bam_bai_input,
        ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        scrpt="cnv_sv/manta/{donor}/runWorkflow.py",
    output:
        cand_si_vcf="cnv_sv/manta/{donor}/results/variants/candidateSmallIndels.vcf.gz",
        cand_si_tbi="cnv_sv/manta/{donor}/results/variants/candidateSmallIndels.vcf.gz.tbi",
        cand_sv_vcf="cnv_sv/manta/{donor}/results/variants/candidateSV.vcf.gz",
        cand_sv_tbi="cnv_sv/manta/{donor}/results/variants/candidateSV.vcf.gz.tbi",
        dipl_sv_vcf="cnv_sv/manta/{donor}/results/variants/diploidSV.vcf.gz",
        dipl_sv_tbi="cnv_sv/manta/{donor}/results/variants/diploidSV.vcf.gz.tbi",
        som_sv_vcf="cnv_sv/manta/{donor}/results/variants/somaticSV.vcf.gz",
        som_sv_tbi="cnv_sv/manta/{donor}/results/variants/somaticSV.vcf.gz.tbi",
        wrk_dir=directory("cnv_sv/manta/{donor}/workspace"),
    log:
        "logs/manta/{donor}/manta_tn.log",
    threads: 8
    conda:
        "../wrappers/manta/env.yaml"
    shell:
        "{input.scrpt} "
        "-j {threads} "
        "-g unlimited &> {log}"
