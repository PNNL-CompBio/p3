

library(knitr)
knitr::opts_chunk$set(cache=TRUE, cache.lazy=FALSE, message=FALSE, warning=FALSE, eval=TRUE)

#+

library(MSnID)
library(PlexedPiper)
library(dplyr)
library(stringr)

#+ 

data_package_num <- 3718

fractions <- get_job_records_by_dataset_package(data_package_num) %>%
  dplyr::select(Dataset) %>%
  arrange(Dataset) %>%
  unique() %>%
  mutate(PlexID = str_extract(Dataset, "^\\d+"))

# 
samples <- readxl::read_xlsx("./../PASS1B_T66_lung_metadata.xlsx") %>% 
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

write.table(fractions,  "fractions.txt",  quote=F, sep="\t", row.names=F)
write.table(samples,    "samples.txt",    quote=F, sep="\t", row.names=F)
write.table(references, "references.txt", quote=F, sep="\t", row.names=F)
