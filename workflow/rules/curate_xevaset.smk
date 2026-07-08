rule build_XevaSet:
    input:
        # 1. The Multi-Omics Backbone
        multiAssayExperiment = rules.build_MultiAssayExperiment.output.mae,
        
        # 2. Metadata & Mapping
        modelInfo = rules.build_drug_screen.output.model,
        treatmentMetadata = rules.annotate_drugs.output,
        modToBiobaseMap = rules.build_drug_screen.output.mod_map,
        
        # 3. Growth Curves & Design
        experiment = rules.build_drug_screen.output.exp,
        expDesign = rules.build_drug_screen.output.design,
        
        # 4. YOUR UPDATED METRICS (The Python Outputs)
        # We use your validated batch and model files here
        # batchMetrics = rules.calculate_drug_metrics.output.batch,
        # modelMetrics = rules.calculate_drug_metrics.output.model
    output:
        xeva = results / "UHN_Breast_XevaSet_v2025.rds",
        tsvDir = directory(results / "UHN_Breast_XevaSet_tsv")
    log:
        logs / "build_XevaSet.log"
    resources:
        mem_mb = 16000  # High memory for multi-modal data integration
    script:
        "../scripts/curation/xevaset/build_XevaSet.R"