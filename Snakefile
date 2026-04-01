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