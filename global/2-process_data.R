
#+ Knits opts

library(knitr)
knitr::opts_chunk$set(cache=F, cache.lazy=FALSE, message=FALSE, warning=FALSE, eval=TRUE)

#+ Load packages

if(!require(MSnID))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSnID")

library(devtools)
if(!require(PlexedPiper))
  devtools::install_github("vladpetyuk/PlexedPiper", build_vignettes = TRUE)

library(dplyr)

source("../util/synapseUtil.R")
#+ Data processing parameters

data_package_num <- 3719
organism_name <- "Rattus norvegicus"
study_name <- "MoTrPAC_PASS1B_Lung_Global"

#+ Read study design tables

frac_id='syn25005573'
samp_id='syn25005574'
ref_id='syn25005575'
synapse_project_id <- 'syn25005572'

fractions <- read.delim(syn$get(frac_id)$path,
                        stringsAsFactors = FALSE, 
                        colClasses = "character")

message("   + Read samples.txt")
samples <- read.delim(syn$get(samp_id)$path,
                      stringsAsFactors = FALSE, 
                      colClasses = "character")

message("   + Read reference.txt")
references <- read.delim(syn$get(ref_id)$path,
                         stringsAsFactors = FALSE, 
                         colClasses = "character")

#+ Process MS/MS data

message("- Prepare MS/MS IDs")
message("   + Read the MS-GF+ output")
msnid <- read_msms_data_from_DMS(data_package_num)

message("   + Correct for isotope selection error")
msnid <- correct_peak_selection(msnid)

message("   + MS/MS ID filter and peptide level")
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)

message("   + MS/MS compute num peptides per 1000 amino acids")
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)


message("   + MS/MS ID filter at protein level")
msnid <- filter_msgf_data_protein_level(msnid, 0.01)

message("   + Inference of parsimonious protein set")
msnid1 <- infer_parsimonious_accessions(msnid, unique_only=T)
msnid2 <- infer_parsimonious_accessions(msnid)

message("   + Remove decoy accessions")
msnid <- apply_filter(msnid, "!isDecoy")

#+ Process MASIC data

message("- Prepare reporter ion intensities")
message("   + Read MASIC ouput")
masic_data <- read_masic_data_from_DMS(data_package_num,
                                       interference_score = TRUE)

message("   + Filtering MASIC data")
masic_data <- filter_masic_data(masic_data, 0.5, 0)

save.image(".RData")

#+ Make crosstabs

message("- Create Reporter Ion Intensity Results")
rii_peptide <- make_rii_peptide_gl(msnid, masic_data, 
                                   fractions, samples, references, 
                                   org_name = organism_name)

message("- Create Ratio Results")
ratio_results <- make_results_ratio_gl(msnid, masic_data, 
                                       fractions, samples, references, 
                                       org_name = organism_name)

message("- Save results")

write.table(rii_peptide, paste(study_name, "results_RII-peptide.txt", sep="_"),
            sep="\t", row.names = FALSE, quote = FALSE)

write.table(ratio_results, paste(study_name, "results_ratio.txt", sep="_"),
            sep="\t", row.names = FALSE, quote = FALSE)

message("- Save to synapse")
synStore(paste(study_name, "results_RII-peptide.txt", sep="_"),synapse_project_id)
synStore(paste(study_name, "results_ratio.txt", sep="_"),synapse_project_id)
message("- Done!")

