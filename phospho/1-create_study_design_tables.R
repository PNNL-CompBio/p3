

library(knitr)
knitr::opts_chunk$set(cache=TRUE, cache.lazy=FALSE, message=FALSE, warning=FALSE, eval=TRUE)

#+

if(!require(MSnID))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSnID")

library(devtools)
if(!require(PlexedPiper))
  devtools::install_github("vladpetyuk/PlexedPiper", build_vignettes = TRUE)

library(dplyr)
library(stringr)

source("../util/synapseUtil.R")

##login to synapse
syn <- synapseLogin()


#+ 

data_package_num <- 3718
metadata_file_id <- ''
synapse_project_id <- 'syn25005572'

fractions <- get_job_records_by_dataset_package(data_package_num) %>%
  dplyr::select(Dataset) %>%
  arrange(Dataset) %>%
  unique() %>%
  mutate(PlexID = str_extract(Dataset, "^\\d+"))

# 
metadata_file <-syn$get(metadata_file_id)$path

samples <- readxl::read_xlsx(metadata_file) %>% 
  dplyr::select(vialLabel, plex, channel)  %>%
  mutate(PlexID = paste0("0", plex),
         ReporterAlias = vialLabel,
         ReporterName = ifelse(channel=="126C", "126", channel),
         QuantBlock = 1,
         MeasurementName = vialLabel) %>%
  dplyr::select(-vialLabel, -plex, -channel)

references <- samples %>% 
  dplyr::select(PlexID) %>%
  distinct() %>%
  mutate(ReporterAlias = "ref",
         ReporterName = "131C",
         QuantBlock = 1,
         MeasurementName = NA)

samples <- rbind(samples, references)

references <- references %>% 
  dplyr::select(PlexID, QuantBlock, ReporterAlias) %>%
  dplyr::rename(Reference = ReporterAlias)

#+

write.table(fractions,  "phospho_fractions.txt",  quote=F, sep="\t", row.names=F)
write.table(samples,    "phospho_samples.txt",    quote=F, sep="\t", row.names=F)
write.table(references, "phoshpo_references.txt", quote=F, sep="\t", row.names=F)

synapseStore('phospho_fractions.txt',synapse_project_id)
synapseStore('phospho_samples.txt',synapse_project_id)
synapseStore('phospho_references.txt',synapse_project_id)