
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

#+ Load processed MS-GF+ data
msnid <- readRDS(file.path(path_to_data, "msnid_filtered.Rds"))

#+ Load processed MASIC data
msnid <- readRDS(file.path(path_to_data, "msnid_filtered.Rds"))

#+ Load study design tables
study_design <- get_study_design_by_dataset_package(data_package_num)

fractions <- study_design$fractions
samples <- study_design$samples
references <- study_design$references

#+ Create accession-level crosstab
crosstab <- create_crosstab(msnid, masic_data,
                            aggregation_level = "accession",
                            fractions, samples, references)

m <- create_msnset(crosstab, samples)

m <- m[,which(pData(m)$MeasurementName != "AML_Reference_F")]

#+ Save to file
saveRDS(m, file = file.path(path_to_data, "msnset_gl.Rds"))

#+ Create peptide-level crosstab
crosstab <- create_crosstab(msnid, masic_data,
                            aggregation_level = c("accession","peptide"),
                            fractions, samples, references)

m <- create_msnset(crosstab, samples)

m <- m[,which(pData(m)$MeasurementName != "AML_Reference_F")]

#+ Save to file
saveRDS(m, file = file.path(path_to_data, "msnset_gl_peptide.Rds"))

