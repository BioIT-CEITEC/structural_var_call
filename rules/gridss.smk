

rule gridss_computation:
    input:  bams = expand("mapped/{sample_name}.bam",sample_name=sample_tab.sort_values(by=['tumor_normal']).sample_name.tolist()),
            ref = config["organism_fasta"], #defined in bioroots utilities
            ref_fai = config["organism_fasta"] + ".fai", #defined in bioroots utilities
            bwa_index = config["organism_bwa"] #defined in bioroots utilities
    output: vcf="structural_varcalls/all_samples/gridss/all_sample_variants.vcf"
    params: tumor_normal_paired = config["tumor_normal_paired"]
    log:    "logs/all_samples/gridss/gridss_computation.log",
    threads: 16
    conda:  "../wrappers/gridss/env.yaml"
    script: "../wrappers/gridss/script.py"

rule gridss_get_per_sample_res:
    input:  all_sample_vcf="structural_varcalls/all_samples/gridss/all_sample_variants.vcf"
    output: vcf="structural_varcalls/{sample_name}/gridss/result_SV.vcf",
    log:    "logs/{sample_name}/gridss/get_per_sample_res.log",
    shell:
        "touch {output.vcf}"