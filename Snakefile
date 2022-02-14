import os
import pandas as pd

# GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_REF_PATH = "/home/rj/4TB/CEITEC/"

##### Config processing #####
#conversion from new SEQUIA json
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
# for CNV pipeline to get corresponding germinal sample column
for index, row in sample_tab.iterrows():
    donorvalue = sample_tab.loc[index,"donor"]
    sample_tab.loc[index,"germinal"] = sample_tab.loc[(sample_tab["donor"] == donorvalue) & (sample_tab["tumor_normal"]=="normal"),"sample_name"].to_string(index=False)    


##### Reference processing
#
config["material"] = "DNA"
if config["lib_ROI"] != "wgs" or config["lib_ROI"] != "RNA":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
else:
    if config["lib_ROI"] != "RNA":
        config["material"] = "RNA" # ??? nemelo by DNA
        config["lib_ROI"] = "wgs"

#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]


#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

####################################
# DEFAULT VALUES
if not "format" in config:
    config["format"] = "default"
if not "min_variant_frequency" in config:
    config["min_variant_frequency"] = 0
if not "not_use_merged" in config:
    config["not_use_merged"] = False



####################################
# SEPARATE RULES
include: "rules/cnvkit.smk"
include: "rules/gatk_cnv.smk"
include: "rules/manta.smk"
include: "rules/pindel.smk"
include: "rules/svdb.smk"
# include: "rules/vep.smk"

####################################
# RULE ALL
rule all:
    input:  
        vcf_cnvkit=lambda wildcards: expand("cnv_sv/cnvkit_vcf/{sample_name}.vcf", sample_name = sample_tab.sample_name),
        segment_regions=lambda wildcards: expand("cnv_sv/cnvkit_batch/{donor}/{sample_name}.cnr",sample_name = sample_tab.sample_name,donor=sample_tab.donor),
        vcf_gatk=lambda wildcards: expand("cnv_sv/gatk_cnv_vcf/{sample_name}.vcf", sample_name = sample_tab.sample_name),
        manta_som_sv_vcf=lambda wildcards: expand("cnv_sv/manta/{donor}/results/variants/somaticSV.vcf.gz",donor = sample_tab.donor),
        pindel_vcf=lambda wildcards: expand("cnv_sv/pindel/{donor}.vcf",donor = sample_tab.donor),
        svdb=lambda wildcards: expand("cnv_sv/svdb_query/{sample_name}.svdb_query.vcf",sample_name = sample_tab.sample_name),