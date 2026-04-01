# Target rule for this module
rule all_omics:
    input:
        procdata / "omics" / "rna_tpm_hugo.tsv",
        procdata / "omics" / "cnv_hugo.tsv",
        procdata / "omics" / "mutation_hugo.tsv"

rule format_rnaseq:
    input:
        rna = rawdata / "omics" / "20251007_gene_tpm_normalized_matrix.csv",
        model_meta = rawdata / "omics" / "omics_raw_metadata_v2025.xlsx",
        info = rawdata / "additional" / "geneInfo.tab"
    output:
        matrix = procdata / "omics" / "rna_tpm_hugo.tsv",
        features = procdata / "omics" / "feature_metadata_rna.tsv"
    log:
        logs / "omics" / "format_rnaseq.log"
    shell:
        "python {scripts}/preprocessing/omics/format_RNAseq.py"

rule format_cnv:
    input:
        cnv = rawdata / "omics" / "20251121_cnv_matrix.tsv",
        info = rawdata / "additional" / "geneInfo.tab"
    output:
        matrix = procdata / "omics" / "cnv_hugo.tsv",
        features = procdata / "omics" / "feature_metadata_cnv.tsv"
    log:
        logs / "omics" / "format_cnv.log"
    shell:
        "python {scripts}/preprocessing/omics/format_cnv.py"

rule format_mutation:
    input:
        maf = rawdata / "omics" / "cohort_20260217_gnomad0002_vaf04_alt5_updaterescue_allgermline.maf.oncokb.txt",
        info = rawdata / "additional" / "geneInfo.tab"
    output:
        matrix = procdata / "omics" / "mutation_hugo.tsv",
        features = procdata / "omics" / "feature_metadata_mutation.tsv"
    log:
        logs / "omics" / "format_mutation.log"
    shell:
        "python {scripts}/preprocessing/omics/format_mutation.py"