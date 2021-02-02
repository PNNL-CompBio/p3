

library(knitr)
knitr::opts_chunk$set(cache=F, cache.lazy=FALSE, message=FALSE, warning=FALSE, eval=TRUE)

#+

library(MSnID)
library(PlexedPiper)
library(dplyr)

#+ Data processing parameters

data_package_num <- 3718
organism_name <- "Rattus norvegicus"
data_package_name <- "MoTrPAC_PASS1B_Lung_Phospho"

path_to_crosstab_gl <- "./../global/MoTrPAC_PASS1B_Lung_Global_results_ratio.txt"

#+ 

fractions <- read.delim("fractions.txt",
                        stringsAsFactors = FALSE, 
                        colClasses = "character")

message("   + Read samples.txt")
samples <- read.delim("samples.txt",
                      stringsAsFactors = FALSE, 
                      colClasses = "character")

message("   + Read reference.txt")
references <- read.delim("references.txt",
                         stringsAsFactors = FALSE, 
                         colClasses = "character")

#+

message("- Prepare MS/MS IDs")
message("   + Read the MS-GF+ output")
msnid <- read_msms_data_from_DMS(data_package_num)

message("   + Read Ascore output")
ascore <- get_AScore_results(data_package_num)
msnid <- best_PTM_location_by_ascore(msnid, ascore)
msnid <- apply_filter(msnid, "grepl(\"\\\\*\", peptide)")

message("   + FDR filter")
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)

crosstab_gl <- read.table(path_to_crosstab_gl, header=T)

# phospho data has isoform, global does not
# solution: give priority to all isoforms
res <- psms(msnid) %>%
  dplyr::select(accession) %>%
  distinct() %>%
  mutate(protein_id = gsub("\\.\\d", "", accession))
res <- crosstab_gl %>%
  dplyr::select(protein_id) %>%
  left_join(res)
res <- res$accession


message("   + Inference of parsimonius set")
msnid <- infer_parsimonious_accessions(msnid, unique_only=F, prior=res)

message("   + Mapping sites to protein sequence")
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)
fst <- Biostrings::readAAStringSet(path_to_FASTA)
names(fst) <- sub("^([A-Z]P_\\d+\\.\\d+)\\s.*", "\\1", names(fst))



msnid <- map_mod_sites(msnid, 
                       fst, 
                       accession_col = "accession", 
                       peptide_mod_col = "Peptide", 
                       mod_char = "*",
                       site_delimiter = "lower")

message("   + Remove decoy sequences")
msnid <- apply_filter(msnid, "!isDecoy")

#+ 

message("- Prepare reporter ion intensities")
message("   + Read MASIC ouput")
masic_data <- read_masic_data_from_DMS(data_package_num,
                                       interference_score =TRUE)

message("   + Filtering MASIC data")
masic_data <- filter_masic_data(masic_data, 0.5, 0)

save.image(".RData")
#+

message("- Create Reporter Ion Intensity Results")
rii_peptide <- make_rii_peptide_ph(msnid, masic_data, 
                                   fractions, samples, references, 
                                   org_name = organism_name)

message("- Create Ratio Results")

library(tidyverse)
make_results_ratio_ph_debug <- function (msnid, masic_data, fractions, samples, references, 
          org_name = "Rattus norvegicus", sep = "_") 
{
  aggregation_level <- c("SiteID")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, 
                              fractions, samples, references)
  crosstab <- crosstab %>% as.data.frame() %>% rownames_to_column("ptm_id") %>% 
    mutate(protein_id = sub("(.P_.*)\\.\\d+-.*", "\\1", ptm_id))
  x <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  y <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID")
  conv <- inner_join(x, y) %>% dplyr::rename(protein_id = REFSEQ, 
                                             gene_symbol = SYMBOL, entrez_id = ENTREZID)
  crosstab <- inner_join(crosstab, conv)
  crosstab <- crosstab %>% dplyr::select(protein_id, ptm_id, 
                                         gene_symbol, entrez_id, everything())
  ascore <- dplyr::select(psms(msnid), protein_id = accession, 
                          ptm_id = SiteID, confident_score = maxAScore) %>% mutate(protein_id = sub("(.P_.*)\\.\\d+", 
                                                                                                    "\\1", protein_id)) %>% group_by(protein_id, ptm_id) %>% 
    summarize(confident_score = max(confident_score))
  crosstab <- left_join(crosstab, ascore)
  crosstab <- crosstab %>% mutate(confident_site = dplyr::case_when(confident_score >= 
                                                                      17 ~ TRUE, confident_score < 17 ~ FALSE),
                                  ptm_id = sub("-", sep, ptm_id)) %>%
    dplyr::select(ptm_id, protein_id, gene_symbol, entrez_id, confident_score, confident_site, everything())
  crosstab[, c(7:ncol(crosstab))] <- signif(crosstab[, c(7:ncol(crosstab))], 
                                            3)
  return(crosstab)
}

ratio_results <- make_results_ratio_ph_debug(msnid, masic_data, 
                                       fractions, samples, references, 
                                       org_name = organism_name)

message("- Save results")

write.table(rii_peptide, paste(data_package_name, "results_RII-peptide.txt", sep="_"),
            sep="\t", row.names = FALSE, quote = FALSE)

write.table(ratio_results, paste(data_package_name, "results_ratio.txt", sep="_"),
            sep="\t", row.names = FALSE, quote = FALSE)

message("- Done!")