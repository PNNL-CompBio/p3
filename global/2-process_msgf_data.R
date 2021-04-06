
data_package_num <- NULL
data_folder <- "./data"
organism_name <- "Homo sapiens"

#+ Libraries, message=F, warning=F
library(PlexedPiper)

#+ Check connection
if (is.null(data_package_num)) {
  stop("Missing data package number.")
} else if(!is_PNNL_DMS_connection_successful())  {
  stop("There is no connection to PNNL DMS. This code won't be evaluated.")
}

#+ Read MS-GF+ data from the DMS
msnid <- read_msms_data_from_DMS(data_package_num)

save(msnid, file="data/msgfData_original.RData")

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

#+ Remove decoys
msnid <- apply_filter(msnid, "!isDecoy")

#+ Inference of parsimonious accessions
msnid <- infer_parsimonious_accessions(msnid, unique_only=FALSE)
show(msnid)

#+ Save filtered MS-GF+ data to file
save(msnid, file="data/msgfData_filtered.RData")
