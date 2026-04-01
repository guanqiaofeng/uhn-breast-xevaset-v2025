rule make_RNASeq_SE:
    input:
        raw = procdata / "omics/rna_tpm_hugo.tsv",
        feature_meta = "data/procdata/omics/feature_metadata_rna.tsv" # Point to your existing file
    output:
        se = procdata / "omics/RNASeq_SE.rds",
        matrix = procdata / "omics/RNASeq_expression.tsv"
    log:
        logs / "omics/make_RNASeq_SE.log"
    script:
        "../scripts/curation/omics/make_RNASeq_SE.R"

rule make_CNV_SE:
    input:
        raw = "data/procdata/omics/cnv_hugo.tsv",
        feature_meta = "data/procdata/omics/feature_metadata_cnv.tsv"
    output:
        se = "data/procdata/omics/CNV_SE.rds",
        matrix = "data/procdata/omics/CNV_matrix.tsv"
    log:
        logs / "omics/make_CNV_SE.log"
    script:
        "../scripts/curation/omics/make_CNV_SE.R"

rule make_mutation_SE:
    input:
        raw = "data/procdata/omics/mutation_hugo.tsv",
        feature_meta = "data/procdata/omics/feature_metadata_mutation.tsv"
    output:
        se = "data/procdata/omics/mutation_SE.rds",
        matrix = "data/procdata/omics/mutation_matrix.tsv"
    log:
        logs / "omics/make_mutation_SE.log"
    script:
        "../scripts/curation/omics/make_mutation_SE.R"

rule build_MultiAssayExperiment:
    input:
        # The individual SummarizedExperiment objects we just curated
        assays = [
            "data/procdata/omics/RNASeq_SE.rds",
            "data/procdata/omics/CNV_SE.rds",
            "data/procdata/omics/mutation_SE.rds"
        ]
    output:
        mae = "data/procdata/omics/UHN_Breast_MAE_v2025.rds"
    log:
        "logs/omics/build_MAE.log"
    script:
        "../scripts/curation/omics/build_MultiAssayExperiment.R"