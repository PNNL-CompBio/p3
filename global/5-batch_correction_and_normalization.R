
# # To install vp.misc:
# library(remotes)
# remotes::install_github("vladpetyuk/vp.misc", build_vignettes=F)

#+ Libraries, message=F, warning=F
library(vp.misc)
library(sva)
library(dplyr)
library(tibble)

#+ Load data
load("data/msnset_global_gene.RData")
load("data/msnset_global_peptide.RData")

#+ Normalize by sample medians
normalizeGlobalDataBySampleMedians <- function(msnset_global) {
  global_coeffs <- apply(exprs(msnset_global), 2, median, na.rm=T)
  exprs(msnset_global) <- sweep(exprs(msnset_global), 2, global_coeffs, "-")
  return(msnset_global)
}

msnset_global_gene <- normalizeGlobalDataBySampleMedians(msnset_global_gene)
msnset_global_peptide <- normalizeGlobalDataBySampleMedians(msnset_global_peptide)

#+ Missing value filter by proportion
filterByProportionMissingValues <- function(msnset, least_proportion_threshold = 0.5) {
  sufficiently_present_features <- which(apply(!is.na(exprs(msnset)), 1, mean) >= least_proportion_threshold)
  msnset <- msnset[sufficiently_present_features,]
}

msnset_global_gene <- filterByProportionMissingValues(msnset_global_gene)
msnset_global_peptide <- filterByProportionMissingValues(msnset_global_peptide)

#+ Missing value filter by plex
filterByMissingPerBatch <- function(msnset, batch_name, least_count_threshold=2) {
  batch_to_sample <- pData(msnset) %>%
    select(removed_cov_name) %>%
    rownames_to_column("sample_name")
  
  sufficiently_present_features <- exprs(msnset) %>%
    as.data.frame() %>% 
    rownames_to_column("feature_name") %>%
    gather(sample_name, abundance, -feature_name) %>%
    inner_join(batch_to_sample, by = "sample_name") %>%
    group_by_at(c(removed_cov_name, "feature_name")) %>% 
    summarize(cnt = sum(!is.na(abundance))) %>%
    group_by_at("feature_name") %>% 
    summarize(min_cnt = min(cnt)) %>% 
    filter(min_cnt >= least_count_threshold) %>%
    pull(feature_name)
  
  msnset <- msnset[sufficiently_present_features, ]
}

msnset_global_gene <- filterByMissingPerBatch(msnset_global_gene, "PlexID")
msnset_global_peptide <- filterByMissingPerBatch(msnset_global_peptide, "PlexID")


#+ Batch correction using ComBat
correctBatchEffectComBat <- function (msnset, batch_name, BPPARAM = bpparam(), ...) 
{
  modcombat <- model.matrix(~1, data = select(pData(msnset), batch_name))
  combat_edata <- ComBat(dat = exprs(msnset),
                         batch = pData(msnset)[, batch_name],
                         mod = modcombat,
                         par.prior = FALSE, prior.plots = TRUE, 
                         BPPARAM = BPPARAM, ...)
  exprs(msnset) <- combat_edata
  return(msnset)
}


msnset_global_gene <- correctBatchEffectComBat(msnset_global_gene, "PlexID")
msnset_global_peptide <- correctBatchEffectComBat(msnset_global_peptide, "PlexID")

#+ Medpolish
library(vp.misc)
msnset_global_gene <- normalizeByGlob(msnset_global_gene, method="medpolish")
msnset_global_peptide <- normalizeByGlob(msnset_global_peptide, method="medpolish")

#+ Save results
save(msnset_global_gene, file="data/msnset_global_gene_corrected.RData"))
save(msnset_global_peptide, file="data/msnset_global_peptide_corrected.RData"))

write.table(exprs(msnset_global_gene), file="crosstab_global_gene_corrected.txt",
            quote=F, sep="\t", row.names=F)
write.table(exprs(msnset_global_peptide), file="crosstab_global_peptide_corrected.txt",
            quote=F, sep="\t", row.names=F)

