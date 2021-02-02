
# this script checks burn in period for global normalization coefficients

library(knitr)
knitr::opts_chunk$set(cache=F, cache.lazy=FALSE, message=FALSE, warning=FALSE, eval=TRUE)

#+

load(".RData")

library(PlexedPiper)
library(stringr)
library(dplyr)

msnid$Fraction <- as.numeric(sub("f", "", str_extract(msnid$Dataset, "f\\d\\d")))
results <- data.frame()
for (i in 1:13) {
  
  msnid_small <- apply_filter(msnid, "Fraction <= i")
  masic_data_small <- masic_data %>% filter(Dataset %in% msnid_small$Dataset)
  
  crosstab <- create_crosstab(msnid_small, masic_data_small, aggregation_level=c("accession"),
                       fractions, samples, references)
  
  medians <- apply(crosstab, 2, median, na.rm=T)
  results <- rbind(results, medians)
}
library(tidyverse)


results_norm <- sweep(results, 2, as.numeric(results[13,]), "-")

x <- results_norm %>%
  rownames_to_column("Fractions") %>%
  mutate(Fractions = as.numeric(Fractions)) %>%
  pivot_longer(!Fractions)

ggplot(x,aes(Fractions, value, group=1)) +
  ggtitle("pass1b lung global norm. coeff. burn-in test") +
  xlab("n") + ylab("sample median of fractions 1-n") +
  stat_summary(geom = "line", fun.y = mean) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3)



