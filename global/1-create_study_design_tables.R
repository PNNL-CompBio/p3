

data_package_num <- NULL

#+ Libraries, warning=F, message=F
library(PlexedPiper)
library(dplyr)
library(stringr)

# Below is example code. It can vary between studies.

# 
# fractions <- get_job_records_by_dataset_package(data_package_num) %>%
#   dplyr::select(Dataset) %>%
#   arrange(Dataset) %>%
#   unique() %>%
#   mutate(PlexID = str_extract(Dataset, "^\\d+"))
# 
# samples <- readxl::read_xlsx("./../PASS1B_T66_lung_metadata.xlsx") %>% 
#   dplyr::select(vialLabel, plex, channel)  %>%
#   mutate(PlexID = paste0("0", plex),
#          ReporterAlias = vialLabel,
#          ReporterName = ifelse(channel=="126C", "126", channel),
#          QuantBlock = 1,
#          MeasurementName = vialLabel) %>%
#   dplyr::select(-vialLabel, -plex, -channel)
# 
# references <- samples %>% select(PlexID) %>%
#   distinct() %>%
#   mutate(ReporterAlias = "ref",
#          ReporterName = "131C",
#          QuantBlock = 1,
#          MeasurementName = NA)
# 
# samples <- rbind(samples, references)
# 
# # creating references is optional
#
# references <- references %>% select(PlexID, QuantBlock, ReporterAlias) %>%
#   dplyr::rename(Reference = ReporterAlias)


# Save to file

write.table(fractions,  
            file="study_design/fractions.txt", 
            quote=F, sep="\t", row.names=F)
write.table(samples,  
            file="study_design/samples.txt", 
            quote=F, sep="\t", row.names=F)

if (exists("references")) {
  write.table(references,  
              file="study_design/references.txt", 
              quote=F, sep="\t", row.names=F)
}


