import os
import pandas as pd


#testhh
configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]
# GLOBAL_REF_PATH = "/home/rj/4TB/CEITEC/"

##### Config processing #####
#conversion from json
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

# for CNV pipeline to get corresponding germinal sample column
# for index, row in sample_tab.iterrows():
#     donorvalue = sample_tab.loc[index,"donor"]
#     sample_tab.loc[index,"germinal"] = sample_tab.loc[(sample_tab["donor"] == donorvalue) & (sample_tab["tumor_normal"]=="normal"),"sample_name"].to_string(index=False)

#### Reference info processing

#### Setting reference from lib_ROI
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]

#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# ####################################
# # VARIALBES FROM CONFIG
used_cnv_callers = []
if config["use_gatk_cnv"]:
    used_cnv_callers.append("gatk_cnv")
if config["use_cnvkit"]:
    used_cnv_callers.append("cnvkit")
if config["use_jabCoNtool"]:
    used_cnv_callers.append("jabCoNtool")
if config["use_control_freec"]:
    used_cnv_callers.append("control_freec")

wildcard_constraints:
     tumor_normal = "tumor|normal",


####################################
# SEPARATE RULES
include: "rules/cnvkit.smk"
include: "rules/gatk_cnv.smk"
include: "rules/jabCoNtool.smk"
include: "rules/control_freec.smk"
include: "rules/manta.smk"
include: "rules/pindel.smk"
include: "rules/svdb.smk"
include: "rules/SV_postprocessing.smk"
# include: "rules/vep.smk"



####################################
# RULE ALL
rule all:
    input:
        final_report="cnv_sv/final_report.html",
        # all_res_prob_tab="variant_calls/all_samples/jabCoNtool/final_CNV_probs.tsv"