
#+ Libraries, message=F, warning=F
library(PlexedPiper)

#+ Read study design
study_design <- read_study_design("study_design/")

fractions <- study_design$fractions
samples <- study_design$samples
references <- study_design$references

#+ Load MS-GF+ data
load("data/msgfData_filtered.RData"))

#+ Load MASIC data
load("data/masicData_filtered.RData"))

#+ Create accession-level crosstab
crosstab <- create_crosstab(msnid, masic_data,
                            aggregation_level=c("accession"),
                            fractions, samples, references)

write.table(crosstab, file="crosstab_global_gene_original.txt"),
            quote=F, sep="\t", row.names=F)

#+ Create accession-level MSnSet
msnset_global_gene <- create_msnset(as.matrix(crosstab), samples)
save(msnset_global_gene, file="msnsets/msnset_global_gene.RData"))

#+ Create peptide-level crosstab
crosstab <- create_crosstab(msnid, masic_data,
                            aggregation_level=c("accession", "peptide"),
                            fractions, samples, references)

write.table(crosstab, file="crosstab_global_peptide_original.txt",
            quote=F, sep="\t", row.names=F)

msnset_global_peptide <- create_msnset(as.matrix(crosstab), samples)
save(msnset_global_peptide, file="msnsets/msnset_global_peptide.RData"))

