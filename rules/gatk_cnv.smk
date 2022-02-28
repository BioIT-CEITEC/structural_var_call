#
# rule gatk_cnv_collect_allelic_counts:
#     input:
#         bam=single_bam_input,
#         bai=single_bam_bai_input,
#         interval=expand("{ref_dir}/other/snp/{ref_name}.{lib_ROI}.af-only-gnomad.intervals",ref_dir=reference_directory,ref_name=config["reference"],lib_ROI=config["lib_ROI"])[0], # "reference/gnomad_SNP_0.001_target.annotated.interval_list" - /mnt/ssd/ssd_3/references/homsap/GRCh37-p13/other/snp/GRCh37-p13.snp.bed.new
#         ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#     output:
#         "cnv_sv/gatk_cnv_collect_allelic_counts/{sample_name}.clean.allelicCounts.tsv",
#     params:
#         extra="",
#     log:
#         "logs/gatk_cnv_collect_allelic_counts/{sample_name}.clean.allelicCounts.tsv.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env.yaml"
#     shell:
#         "(gatk --java-options '-Xmx8g' CollectAllelicCounts "
#         "-L {input.interval} "
#         "-I {input.bam} "
#         "-R {input.ref} "
#         "-O {output}"
#         "{params.extra}) &> {log}"
#
# rule gatk_cnv_collect_read_counts:
#     input:
#         bam=single_bam_input,
#         bai=single_bam_bai_input,
#         interval=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
#     output:
#         "cnv_sv/gatk_cnv_collect_read_counts/{sample_name}.counts.hdf5",
#     params:
#         mergingRule="OVERLAPPING_ONLY",
#         extra="",
#     log:
#         "logs/gatk_cnv_collect_read_counts/{sample_name}.counts.hdf5.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env.yaml"
#     shell:
#         "(gatk --java-options '-Xmx8g' CollectReadCounts "
#         "-I {input.bam} "
#         "-L {input.interval} "
#         "--interval-merging-rule {params.mergingRule} "
#         "{params.extra} "
#         "-O {output}) &> {log}"
#
# rule gatk_cnv_call_copy_ratio_segments:
#     input:
#         "cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.cr.seg",
#     output:
#         segments="cnv_sv/gatk_cnv_call_copy_ratio_segments/{sample_name}.clean.calledCNVs.seg",
#         igv_segments="cnv_sv/gatk_cnv_call_copy_ratio_segments/{sample_name}.clean.calledCNVs.igv.seg",
#     params:
#         extra="",
#     log:
#         "logs/gatk_cnv_call_copy_ratio_segments/{sample_name}.clean.calledCNVs.seg.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env.yaml"
#     shell:
#         "(gatk --java-options '-Xmx8g' CallCopyRatioSegments "
#         "--input {input} "
#         "--output {output.segments} "
#         "{params.extra}) &> {log}"
#
#
# def normal_read_counts_input(wildcards):
#     return set(expand("cnv_sv/gatk_cnv_collect_read_counts/{sample_name}.counts.hdf5",zip\
#             ,sample_name=sample_tab.loc[sample_tab.tumor_normal == "normal" , "sample_name"].tolist()))
#
# rule gatk_create_panel_of_normals:
#     input:
#         germinal_read_counts = normal_read_counts_input,
#     output:
#         hdf5PoN="cnv_sv/gatk_create_panel_of_normals/panel_of_normals.hdf5"
#     log:
#         "logs/gatk_create_panel_of_normals/gatk_create_panel_of_normals.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env.yaml"
#     script:
#         "../wrappers/gatk/create_panel_of_normals.py"
#
#
# rule gatk_cnv_denoise_read_counts:
#     input:
#         hdf5PoN="cnv_sv/gatk_create_panel_of_normals/panel_of_normals.hdf5",
#         hdf5Tumor="cnv_sv/gatk_cnv_collect_read_counts/{sample_name}.counts.hdf5",
#     output:
#         denoisedCopyRatio="cnv_sv/gatk_cnv_denoise_read_counts/{sample_name}.clean.denoisedCR.tsv",
#         stdCopyRatio="cnv_sv/gatk_cnv_denoise_read_counts/{sample_name}.clean.standardizedCR.tsv",
#     params:
#         extra="",
#     log:
#         "logs/gatk_cnv_denoise_read_counts/{sample_name}.clean.denoisedCR.tsv.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env.yaml"
#     shell:
#         "(gatk --java-options '-Xmx8g' DenoiseReadCounts -I {input.hdf5Tumor} "
#         "--count-panel-of-normals {input.hdf5PoN} "
#         "--standardized-copy-ratios {output.stdCopyRatio} "
#         "--denoised-copy-ratios {output.denoisedCopyRatio} "
#         "{params.extra}) &> {log}"
#
#
# rule gatk_cnv_model_segments:
#     input:
#         denoisedCopyRatio="cnv_sv/gatk_cnv_denoise_read_counts/{sample_name}.clean.denoisedCR.tsv",
#         allelicCounts="cnv_sv/gatk_cnv_collect_allelic_counts/{sample_name}.clean.allelicCounts.tsv",
#     output:
#         "cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelFinal.seg",
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.cr.seg"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.af.igv.seg"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.cr.igv.seg"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.hets.tsv"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelBegin.cr.param"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelBegin.af.param"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelBegin.seg"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelFinal.af.param"),
#         temp("cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelFinal.cr.param"),
#     params:
#         outdir=lambda wildcards, output: os.path.dirname(output[0]),
#         outprefix="{sample_name}.clean",
#         extra="",
#     log:
#         "logs/gatk_cnv_model_segments/{sample_name}.clean.modelFinal.seg.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env.yaml"
#     shell:
#         "(gatk --java-options '-Xmx8g' ModelSegments "
#         "--denoised-copy-ratios {input.denoisedCopyRatio} "
#         "--allelic-counts {input.allelicCounts} "
#         "--output {params.outdir} "
#         "--output-prefix {params.outprefix}"
#         "{params.extra}) &> {log}"
#
# rule gatk_cnv_vcf:
#     input:
#         segment="cnv_sv/gatk_cnv_model_segments/{sample_name}.clean.modelFinal.seg",
#     output:
#         vcf="cnv_sv/gatk_cnv_vcf/{sample_name}.vcf",
#     params:
#         sample_id="{sample_name}",
#         hom_del_limit=0.5, #dat do workflow.json
#         het_del_limit=1.5, #dat do workflow.json
#         dup_limit=2.5, #dat do workflow.json
#         TC=0.5,#lambda wildcards: get_sample(samples, wildcards)["TC"],
#     log:
#         "logs/gatk_cnv_vcf/{sample_name}.vcf.log",
#     threads: 8
#     conda:
#         "../wrappers/gatk/env_python.yaml"
#     script:
#         "../wrappers/gatk/gatk_cnv_vcf.py"