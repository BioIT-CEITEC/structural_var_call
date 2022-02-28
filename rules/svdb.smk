rule svdb_merge:
    input:
        vcfs=expand("variant_calls/{{sample_name}}/{cnv_caller}/CNV.vcf", cnv_caller=["cnvkit"]),
    output:
        vcf="final_CNV_calls/merged/{sample_name}.merged.vcf",
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
        vcf="final_CNV_calls/merged/{sample_name}.merged.vcf",
        svdb_vcf=expand("{ref_dir}/other/svdb/gnomad_v2.1_sv.sites.vcf",ref_dir=reference_directory)[0], # "reference/normal_26_svdb_0.8.vcf"
    output:
        vcf="final_CNV_calls/{sample_name}.svdb_query.vcf",
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0][:-6],
    log:
        "logs/svdb_query/{sample_name}.log",
    threads: 8
    conda:
        "../wrappers/svdb/env.yaml"
    shell:
        "(svdb --query --query_vcf {input.vcf} --db {input.svdb_vcf} --prefix {params.prefix}) &> {log}"
