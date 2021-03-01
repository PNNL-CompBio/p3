
library(knitr)
knitr::opts_chunk$set(cache=FALSE)

#+ Setup, message=F, warning=F
library(MSnID)
library(PlexedPiper)
library(dplyr)

#+ Data processing parameters
data_package_num <- NULL
organism_name <- "Homo sapiens"
path_to_data <- "./data"

#+ Check connection
if (is.null(data_package_num)) {
  stop("Missing data package number.")
} else if(!is_PNNL_DMS_connection_successful())  {
  stop("There is no connection to PNNL DMS. This code won't be evaluated.")
}
  
#+ Create data directory
dir.create(path_to_data)

#+ Read MS-GF+ data from the DMS
msnid <- read_msms_data_from_DMS(data_package_num)

#+ MS/MS peptide level FDR filter
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
show(msnid)

#+ Remap accessions from Uniprot to Gene
msnid <- remap_accessions_uniprot_to_gene(msnid, 
                                          organism_name=organism_name)
show(msnid)

#+ MS/MS protein level FDR filter
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)
path_to_FASTA <- gsub("\\\\","/",path_to_FASTA)
path_to_FASTA_gene <- remap_accessions_uniprot_to_gene_fasta(
  path_to_FASTA)

msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA_gene)

msnid <- filter_msgf_data_protein_level(msnid, 0.01)
show(msnid)

#+ Inference of parsimonious accessions
msnid <- infer_parsimonious_accessions(msnid, unique_only=TRUE)
show(msnid)

#+ Remove decoys
msnid <- apply_filter(msnid, "!isDecoy")

#+ Save filtered MS-GF+ data to file
saveRDS(msnid, file.path(path_to_data, "msnid_filtered.Rds"))

#+ Read MASIC data

masic_data <- read_masic_data_from_DMS(data_package_num,
                                       interference_score = TRUE)
nrow(masic_data)

#+ Filter MASIC data
masic_data <- filter_masic_data(masic_data,
                                interference_score_threshold = 0.5,
                                s2n_threshold = 0)
nrow(masic_data)

#+ Save filtered MASIC data to file
saveRDS(masic_data, file.path(path_to_data, "masic_data_filtered.Rds"))