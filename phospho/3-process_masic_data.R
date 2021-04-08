
data_package_num <- NULL

#+ Libraries, message=F, warning=F
library(PlexedPiper)

#+ Check connection
if (is.null(data_package_num)) {
  stop("Missing data package number.")
} else if(!is_PNNL_DMS_connection_successful())  {
  stop("There is no connection to PNNL DMS. This code won't be evaluated.")
}

#+ Read MASIC data

masic_data <- read_masic_data_from_DMS(data_package_num,
                                       interference_score = TRUE)
nrow(masic_data)

save(masic_data, file="data/masicData_original.RData")

#+ Filter MASIC data
masic_data <- filter_masic_data(masic_data,
                                interference_score_threshold = 0.5,
                                s2n_threshold = 0)
nrow(masic_data)

#+ Save filtered MASIC data to file
save(masic_data, file="data/masicData_filtered.RData")
