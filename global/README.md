Global proteomics pipeline using PlexedPiper
================
Michael Nestor
4/7/2021

# Global Proteomics analysis

The global proteomics pipeline uses the R package
[PlexedPiper](https://github.com/vladpetyuk/PlexedPiper). It also
requires a connection to the DMS to access data packages.

    ## Warning: package 'knitr' was built under R version 4.0.4

    library(PlexedPiper)
    data_package_num <- 3606

    if (!is_PNNL_DMS_connection_successful()) {
      stop("No connection to DMS.")
    }

## 1. Read study design information

Study design information in PlexedPiper is encoded in three tables:
fractions, samples, and references. These tables can be made using
metadata and should be stored on the DMS before processing.

    study_design <- get_study_design_by_dataset_package(data_package_num)

    fractions <- study_design$fractions
    samples <- study_design$samples
    references <- study_design$references

## 2 Processing MS-GF+ data

MS-GF+ data is processed in several steps. First, read MS-GF+ output
from the DMS. (This step can take a while).

    msnid <- read_msms_data_from_DMS(data_package_num)

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    show(msnid)

    ## MSnID object
    ## Working directory: "."
    ## #Spectrum Files:  4 
    ## #PSMs: 96612 at 30 % FDR
    ## #peptides: 53098 at 54 % FDR
    ## #accessions: 72917 at 84 % FDR

### 2.1 Remap accessions

This function remaps UniProt protein accessions to gene symbol.

    msnid <- remap_accessions_refseq_to_gene(msnid,
                                              organism_name="Rattus norvegicus")
    show(msnid)

    ## MSnID object
    ## Working directory: "."
    ## #Spectrum Files:  4 
    ## #PSMs: 96612 at 30 % FDR
    ## #peptides: 53098 at 54 % FDR
    ## #accessions: 23565 at 80 % FDR

### 2.2 FDR filter

We use the target-decoy search strategy method described in [(Elias
2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922680/). Filtering
is done first at peptide level, then at protein level, both with max FDR
of 1%.

    msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
    show(msnid)

    ## MSnID object
    ## Working directory: "."
    ## #Spectrum Files:  4 
    ## #PSMs: 29489 at 0.78 % FDR
    ## #peptides: 10776 at 0.99 % FDR
    ## #accessions: 4069 at 2.8 % FDR

    path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)
    path_to_FASTA <- gsub("\\\\", "/", path_to_FASTA)
    path_to_FASTA_gene <- remap_accessions_refseq_to_gene_fasta(path_to_FASTA,
                                                                organism_name="Rattus norvegicus")

    msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA_gene)
    msnid <- filter_msgf_data_protein_level(msnid, 0.01)
    show(msnid)

    ## MSnID object
    ## Working directory: "."
    ## #Spectrum Files:  4 
    ## #PSMs: 27086 at 0.18 % FDR
    ## #peptides: 9105 at 0.29 % FDR
    ## #accessions: 2628 at 1 % FDR

    msnid <- apply_filter(msnid, "!isDecoy")
    show(msnid)

    ## MSnID object
    ## Working directory: "."
    ## #Spectrum Files:  4 
    ## #PSMs: 27036 at 0 % FDR
    ## #peptides: 9079 at 0 % FDR
    ## #accessions: 2602 at 0 % FDR

### 2.3 Parsimonious inference

To reduce number of protein identifications, we use a parsimonious
inference algorithm described in [(Zhang et
al. 2007)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2810678/).

    msnid <- infer_parsimonious_accessions(msnid)
    show(msnid)

    ## MSnID object
    ## Working directory: "."
    ## #Spectrum Files:  4 
    ## #PSMs: 26251 at 0 % FDR
    ## #peptides: 8916 at 0 % FDR
    ## #accessions: 2135 at 0 % FDR

## 3 Process MASIC data

Output from the MASIC software is read from DMS, then filtered by
inteference score.

    masic_data <- read_masic_data_from_DMS(data_package_num,
                                           interference_score = TRUE)

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    nrow(masic_data)

    ## [1] 133310

    masic_data <- filter_masic_data(masic_data,
                                    interference_score_threshold = 0.5,
                                    s2n_threshold = 0)
    nrow(masic_data)

    ## [1] 127440

## 4 Create crosstab

The quantitative crosstab combines MS/MS identifications with reporter
ion intensities. Abundances are taken relative to the reference channel
and then log-transformed.

    aggregation_level <- c("accession")
    crosstab <- create_crosstab(msnid, masic_data,
                                aggregation_level,
                                fractions, samples, references)

    head(crosstab)

    ##                 R_01        R_02       R_03       R_04       R_05       R_06
    ## Aar2     -0.42583395 -0.22750410 -0.6327384 -0.4048138 -0.7398460 -0.4532292
    ## Aars      0.42828710  0.17127004 -0.1019011  0.2566366 -0.0480728 -0.5617621
    ## Aarsd1   -0.12023845  0.01511791 -0.3802445 -0.1853858 -0.4032233 -0.9480688
    ## Aasdhppt -0.20639226 -0.14694454 -0.4930121 -0.2691141 -0.4189392 -0.3202074
    ## Abca8a    0.28351420 -0.39489529 -0.4669246  0.4997377 -0.1558518 -0.9563149
    ## Abcf2    -0.08554908 -0.07049877 -0.3473433 -0.1480327 -0.4063104 -0.6928486
    ##                R_07       R_08        R_09       S_01       S_02       S_03
    ## Aar2     -0.4393845 -0.1819555  0.01647107 -0.7029618 -0.9735167 -0.6317493
    ## Aars     -0.7707195 -0.2250906 -0.06907425 -0.7724941 -2.3295992 -0.9400888
    ## Aarsd1   -0.4499616 -0.5698316  0.11402759 -0.1131844 -0.7855627 -0.5890304
    ## Aasdhppt -0.6276055 -0.3094378 -0.01016151 -0.2567013 -0.9917061 -0.5689407
    ## Abca8a   -0.4161911 -1.0309648 -0.25237278 -0.7653302 -0.9532367 -0.8305449
    ## Abcf2    -0.5848583 -0.3399608 -0.08172162 -0.3018273 -0.9694549 -0.4949420
    ##                S_04        S_05       S_06         S_07         S_08       S_09
    ## Aar2     -1.1097444 -0.03755291 -0.1883340 -0.009765108 -0.688394321 -0.1993067
    ## Aars     -1.9347177 -0.52464810 -0.2261498  0.196205252 -0.780449091 -0.2027615
    ## Aarsd1   -0.9976901  0.06778754 -0.2433397 -0.322310988 -0.621658332 -0.5719310
    ## Aasdhppt -1.1919043 -0.09455174 -0.3415784 -0.030111499 -0.598307560 -0.3868711
    ## Abca8a   -1.2029846  0.66504287 -0.6686938 -0.694667060 -0.001703917 -0.8668528
    ## Abcf2    -0.8742557  0.10895016 -0.1140746 -0.632424701 -0.613297398 -0.4806694
