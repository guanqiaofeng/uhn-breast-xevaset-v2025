from pathlib import Path

configfile: "config/pipeline.yaml"

# Configuration for Paths
RAW_DIR = "data/rawdata/drug_screen"
ADD_DIR = "data/rawdata/additional"
OMICS_DIR = "data/procdata/omics"
PROC_DIR = "data/procdata/drug_screen"
LOG_DIR = "logs/drug_screen"

rule all_drug_screen:
    input:
        expand(PROC_DIR + "/build/{file}.csv", file=["model", "drug", "experiment", "expDesign", "batch_metrics", "model_metrics"]),
        PROC_DIR + "/update/update_model_drug_batch_data.tsv"

rule scan_drug_screen:
    """Step 1: Scan raw Excel files and generate a checksum manifest."""
    input:
        raw_files = list(Path(RAW_DIR).glob("**/*.xlsx"))
    output:
        manifest = PROC_DIR + "/scan/raw_file_manifest.tsv"
    log:
        LOG_DIR + "/scan.log"
    shell:
        "python workflow/scripts/preprocessing/drug_screen/run_scan.py > {log} 2>&1"

rule process_drug_screen:
    """Step 2: Clean raw data and standardize column names."""
    input:
        manifest = PROC_DIR + "/scan/raw_file_manifest.tsv"
    output:
        processed = PROC_DIR + "/process/processed_data.tsv"
    log:
        LOG_DIR + "/process.log"
    shell:
        "python workflow/scripts/preprocessing/drug_screen/run_process.py > {log} 2>&1"

rule update_drug_screen:
    """Step 3: Apply Regex ID extraction and Drug Mapping."""
    input:
        processed = PROC_DIR + "/process/processed_data.tsv",
        drug_map = ADD_DIR + "/all_drug.xlsx",
        model_map = ADD_DIR + "/all_modelid_sampleid_modeltype_mapping.xlsx"
    output:
        updated = PROC_DIR + "/update/update_model_drug_batch_data.tsv",
        filtered = PROC_DIR + "/update/drug_filtered_out.tsv"
    params:
        model_type = config["drug_screen"]["update"]["model_type"]  # Options: both, resistant, parental
    log:
        LOG_DIR + "/update.log"
    shell:
        "python workflow/scripts/preprocessing/drug_screen/run_update.py -m {params.model_type} > {log} 2>&1"

rule build_drug_screen:
    """Step 4: Generate final XevaSet input CSVs and Omics Mapping."""
    input:
        updated = PROC_DIR + "/update/update_model_drug_batch_data.tsv",
        omics = OMICS_DIR
    output:
        model = PROC_DIR + "/build/model.csv",
        drug = PROC_DIR + "/build/drug.csv",
        exp = PROC_DIR + "/build/experiment.csv",
        design = PROC_DIR + "/build/expDesign.csv",
        mod_map = PROC_DIR + "/build/modToBiobaseMap.csv"
    log:
        LOG_DIR + "/build.log"
    shell:
        "python workflow/scripts/preprocessing/drug_screen/run_build.py > {log} 2>&1"

rule calculate_drug_metrics:
    """Step 5: Calculate TGI (Batch-level) and mRECIST/Response (Model-level)."""
    input:
        experiment = PROC_DIR + "/build/experiment.csv",
        design = os.path.join(PROC_DIR, "build/expDesign.csv")
    output:
        batch = PROC_DIR + "/build/batch_metrics.csv",
        model = PROC_DIR + "/build/model_metrics.csv"
    params:
        tgi_window = config["drug_screen"]["curve_calculation"]["tgi_window"],
        recist_window = config["drug_screen"]["curve_calculation"]["recist_window"]
    log:
        LOG_DIR + "/calculate_metrics.log"
    shell:
        """
        python workflow/scripts/preprocessing/drug_screen/run_calculate_metrics.py \
            --exp {input.experiment} \
            --design {input.design} \
            --out_batch {output.batch} \
            --out_model {output.model} \
            --tgi_win {params.tgi_window} \
            --recist_win {params.recist_window} > {log} 2>&1
        """