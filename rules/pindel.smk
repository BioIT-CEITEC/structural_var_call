
rule picard_insert_size:
    input:
        bam=tumor_bam_input,
    output:
        metrics="cnv_sv/picard_collect_multiple_metrics/{donor}.insert_size_metrics.txt",
        histogram="cnv_sv/picard_collect_multiple_metrics/{donor}.insert_size_histogram.pdf",
    log: "logs/picard/{donor}/CollectInsertSizeMetrics.log"
    threads: 5
    conda:  "../wrappers/picard_insert_size/env.yaml"
    script: "../wrappers/picard_insert_size/script.py"

rule generate_pindel_config:
    input:
        bam=tumor_bam_input,
        metrics="cnv_sv/picard_collect_multiple_metrics/{donor}.insert_size_metrics.txt",
    output:
        config="cnv_sv/pindel/{donor}.cfg",
    log:
        "logs/pindel/{donor}_generate_pindel_config.log",
    threads: 8
    conda:
        "../wrappers/pindel/env.yaml"
    script:
        "../wrappers/pindel/generate_pindel_config.py"

rule pindel:
    input:
        config="cnv_sv/pindel/{donor}.cfg",
        bam=tumor_bam_input,
        bai=tumor_bam_bai_input,
        ref=config["organism_fasta"], #defined in bioroots utilities
    output:
        pindel=expand("cnv_sv/pindel/{{donor}}_{ext}",ext=["BP","CloseEndMapped","D","INT_final","INV","LI","RP","SI","TD",],),
    params:
        prefix=lambda wildcards: "cnv_sv/pindel/%s" % (wildcards.donor),
    log:
        "logs/pindel/{donor}_pindel.log",
    threads: 8
    conda:
        "../wrappers/pindel/env.yaml"
    script:
        "../wrappers/pindel/pindel_call.py"

rule pindel2vcf:
    input:
        pindel=expand("cnv_sv/pindel/{{donor}}_{ext}",ext=["BP","CloseEndMapped","D","INT_final","INV","LI","RP","SI","TD",],),
        ref=config["organism_fasta"], #defined in bioroots utilities
    output:
        vcf="cnv_sv/pindel/{donor}.vcf",
    params:
        refname="hg19", #hg19-GRCh37 WTF dodelat
        refdate="200902", ## WTF dodelat
    log:
        "logs/pindel/{donor}_pindel2vcf.log",
    threads: 1
    conda:
        "../wrappers/pindel/env.yaml"
    script:
        "../wrappers/pindel/pindel2vcf.py"