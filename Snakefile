from pathlib import Path

# Load configuration
configfile: "config/pipeline.yaml"

# Define Path objects for clean syntax
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path(config["directories"]["scripts"])

# Include sub-modules
include: "workflow/rules/preprocess_omics.smk"
include: "workflow/rules/preprocess_drugscreen.smk"
include: "workflow/rules/curate_omics.smk"
include: "workflow/rules/curate_drugscreen.smk"
include: "workflow/rules/curate_xevaset.smk"

rule all:
    input:
        "data/results/UHN_Breast_XevaSet_v2025.rds",
        "data/results/UHN_Breast_XevaSet_tsv"
    localrule: True

# This will be your final assembly rule later
# rule build_XevaSet:
#     input:
#         rna = rules.format_rnaseq.output.matrix,
#         cnv = rules.format_cnv.output.matrix,
#         mut = rules.format_mutation.output.matrix
#     output:
#         xeva = results / "Xeva_UHN_Breast.rds"
#     script:
#         "workflow/scripts/drug_screen/run_build.R"