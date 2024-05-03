rule svdb_merge:
    input:
        vcfs=expand("structural_varcalls/{{sample_name}}/{cnv_caller}/result_SV.vcf", cnv_caller=used_SV_callers),
    output:
        vcf="final_SV_calls/{sample_name}.merged.vcf",
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
        vcf="final_SV_calls/{sample_name}.merged.vcf",
        svdb_vcf = config["organism_svdb"],
    output:
        vcf="final_SV_calls/{sample_name}.svdb_query.vcf",
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0][:-6],
    log:
        "logs/svdb_query/{sample_name}.log",
    threads: 8
    conda:
        "../wrappers/svdb/env.yaml"
    shell:
        "(svdb --query --query_vcf {input.vcf} --in_occ AC --in_frq AC --db {input.svdb_vcf} --prefix {params.prefix}) &> {log}"
