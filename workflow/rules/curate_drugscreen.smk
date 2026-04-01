rule annotate_drugs:
    input:
        "data/procdata/drug_screen/build/drug.csv"
    output:
        "data/procdata/drug_screen/build/drug_annotated.csv"
    script:
        "../scripts/curation/drug_screen/annotate_drugs.R"
