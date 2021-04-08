
data_package_num <- NULL
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


#+ Correct peak selection
msnid <- correct_peak_selection(msnid)


#+ Get AScore results
ascore <- get_AScore_results(data_package_num)
msnid <- best_PTM_location_by_ascore(msnid, ascore)
msnid <- apply_filter(msnid, "grepl(\"\\\\*\", peptide)")


#+ Remap accessions from Uniprot to Gene
msnid <- remap_accessions_uniprot_to_gene(msnid, 
                                          organism_name=organism_name)
show(msnid)


#+ MS/MS peptide level FDR filter
msnid <- filter_msgf_data(msnid, level="peptide", fdr.max=0.01)
show(msnid)
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)


#+ Inference of parsimonious accessions
msnid <- infer_parsimonious_accessions(msnid, unique_only=FALSE)
show(msnid)


#+ Map mod sites
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)
path_to_FASTA <- gsub("\\\\","/",path_to_FASTA)
path_to_FASTA_gene <- remap_accessions_uniprot_to_gene_fasta(path_to_FASTA)

library(Biostrings)
fst <- readAAStringSet(path_to_FASTA_gene)
msnid <- map_mod_sites(msnid, fst, 
                       accession_col = "accession", 
                       peptide_mod_col = "Peptide", 
                       mod_char = "*",
                       site_delimiter = "lower")
msnid <- apply_filter(msnid, "!is.na(PepLocFirst")
head(psms(msnid))


#+ Save filtered MS-GF+ data to file
save(msnid, file="data/msgfData_filtered.RData")
