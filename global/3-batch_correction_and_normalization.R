
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
study_name <- NULL

#+

library(PlexedPiper)
library(MSnbase)
library(dplyr)
library(vp.misc)
library(ggplot2)
library(readxl)
library(tidyverse)


msnset_gl <- readRDS(file.path(path_to_data, "msnset_gl.Rds"))

msnset_gl_peptide <- readRDS(file.path(path_to_data, "msnset_gl_peptide.Rds"))


global.coeffs <- apply(exprs(msnset_gl), 2, median, na.rm=T)

exprs(msnset_gl) <- sweep(msnset_gl, MARGIN=2,
                          STATS=apply(exprs(msnset_gl), 2, median, na.rm=T)
                          FUN="-")

exprs(msnset_gl_peptide) <- sweep(msnset_gl_peptide, MARGIN=2,
                          STATS=apply(exprs(msnset_gl_peptide), 2, median, na.rm=T)
                          FUN="-")

# Batch correct function
# Todo: push this version to vp.misc
correct_batch_effect_empiricalBayesLM <- function (x, removed_cov_name, retained_cov_name = NULL, least_count_threshold=2, 
                                                   ...) 
{
  batch_to_sample <- pData(x) %>% select(removed_cov_name) %>% rownames_to_column("sample_name")
  sufficiently_present_features <- exprs(x) %>% as.data.frame() %>% 
    rownames_to_column("feature_name") %>% gather(sample_name, 
                                                  abundance, -feature_name) %>% inner_join(batch_to_sample, 
                                                                                           by = "sample_name") %>% group_by_at(c(removed_cov_name, "feature_name")) %>% 
    summarize(cnt = sum(!is.na(abundance))) %>% group_by_at("feature_name") %>% 
    summarize(min_cnt = min(cnt)) %>% filter(min_cnt >= 
                                               least_count_threshold) %>% pull(feature_name)
  x <- x[sufficiently_present_features, ]
  
  e <- exprs(x)
  e <- as.data.frame(t(e))
  if (!all(c(removed_cov_name, retained_cov_name) %in% names(pData(x)))) {
    stop("The covariates are not recognized")
  }
  if (is.null(retained_cov_name)) {
    soln <- WGCNA::empiricalBayesLM(data = e, removedCovariates = pData(x)[, 
                                                                           removed_cov_name])
  }
  else {
    soln <- WGCNA::empiricalBayesLM(data = e, removedCovariates = pData(x)[, 
                                                                           removed_cov_name], retainedCovariates = pData(x)[, 
                                                                                                                            retained_cov_name])
  }
  e <- soln[["adjustedData"]]
  exprs(x) <- t(as.matrix(e))
  return(x)
}

# Write unnormalized crosstabs to file


write.table(exprs(msnset_gl), 
            file.path(path_to_data, paste(c(study_name,"crosstab_gl_unnormalized.txt"),collapse="_")),
            sep="\t", quote=F)

write.table(exprs(msnset_gl_peptide),
            file.path(path_to_data, paste(c(study_name,"crosstab_gl_peptide_unnormalized.txt"),collapse="_")),
            sep="\t", quote=F)



# 50% missing value filter

msnset_gl <- msnset_gl[which(apply(!is.na(exprs(msnset_gl)), 1, mean) >= 0.5),]

mnset_gl_peptide <- msnset_gl_peptide[which(apply(!is.na(exprs(msnset_gl_peptide)), 1, mean) >= 0.5),]

# main batch correction step

msnset_gl <- correct_batch_effect_empiricalBayesLM(msnset_gl,
                                           removed_cov_name = "PlexID"))

msnset_gl_peptide <- correct_batch_effect_empiricalBayesLM(msnset_gl_peptide,
                                                   removed_cov_name = "PlexID"))

# Median polish

msnset_gl <- normalizeByGlob(msnset_gl, method="medpolish")
msnset_gl_peptide <- normalizeByGlob(msnset_gl_peptide, method="medpolish")

# Write normalized crosstabs to file

write.table(exprs(msnset_gl),
            file.path(path_to_data, paste(c(study_name,"crosstab_gl_unnormalized.txt"),collapse="_")),
            sep="\t", quote=F)

write.table(exprs(msnset_gl_peptide),
            file.path(path_to_data, paste(c(study_name,"crosstab_gl_peptide_unnormalized.txt"),collapse="_")),
            sep="\t", quote=F)

