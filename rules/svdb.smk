rule svdb_merge:
    input:
        vcfs=expand("cnv_sv/{cnv_caller}_vcf/{{sample_name}}.vcf", cnv_caller=["cnvkit","gatk_cnv"], ),
    output:
        vcf="cnv_sv/svdb_merge/{sample_name}.merged.vcf",
    params:
        overlap=0.6, #config.get("svdb_merge", {}).get("overlap", 0.6),
    log:
        "logs/svdb_merge/{sample_name}.vcf.log",
    threads: 8
    conda:
        "../wrappers/svdb/env.yaml"
    shell:
        "(svdb --merge --vcf {input.vcfs} > {output.vcf}) 2> {log}"
#

# svdb_vcf from DBVAR = https://www.ncbi.nlm.nih.gov/dbvar/content/ftp_manifest/
rule svdb_query:
    input:
        vcf="cnv_sv/svdb_merge/{sample_name}.merged.vcf",
        svdb_vcf=expand("{ref_dir}/other/dbvar/{ref_name}.variant_call.all.vcf",ref_dir=reference_directory,ref_name=config["reference"])[0], # "reference/normal_26_svdb_0.8.vcf"
    output:
        vcf="cnv_sv/svdb_query/{sample_name}.svdb_query.vcf",
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0][:-6],
    log:
        "logs/svdb_query/{sample_name}.log",
    threads: 8
    conda:
        "../wrappers/svdb/env.yaml"
    shell:
        "(svdb --query --query_vcf {input.vcf} --db {input.svdb_vcf} --prefix {params.prefix}) &> {log}"
