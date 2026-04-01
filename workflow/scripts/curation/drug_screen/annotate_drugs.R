# Standard R script for the curation stage
library(AnnotationGx)
library(data.table)
library(dplyr)

# Load the preprocessed (but unannotated) drug file
#drug <- fread(snakemake@input[[1]]) |> as.data.frame()
input_file <- "data/procdata/drug_screen/build/drug.csv"
drug <- fread(input_file) |> as.data.frame()

# Filter for single agents (AnnotationGx works best here)
drugs_to_map <- drug %>%
  filter(is.ADC == "no", is.PDC == "no", is.single == "yes") %>%
  pull(drugname.standardized)

# Map and Fetch
message("📡 Querying PubChem for ", length(drugs_to_map), " compounds...")
compound_2_cids <- mapCompound2CID(drugs_to_map, first = TRUE) %>%
  filter(!is.na(cids)) %>%
  rename(CID = cids)

compound_props <- mapCID2Properties(
  ids = compound_2_cids$CID,
  properties = c("MolecularFormula", "MolecularWeight", "IUPACName", "InChI", "InChIKey", "CanonicalSMILES")
)

# Merge back
compound_info <- left_join(compound_2_cids, compound_props, by = "CID")
drug_annotated <- left_join(drug, compound_info, by = c("drugname.standardized" = "name"))

# Write out the final curated file for the Python BUILD step
#fwrite(drug_annotated, snakemake@output[[1]])
output_file <- "data/procdata/drug_screen/build/drug_annotated.csv"
fwrite(drug_annotated, output_file)

