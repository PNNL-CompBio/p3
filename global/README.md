# Global Proteomics analysis

The global proteomics pipeline uses the R package [PlexedPiper](https://github.com/vladpetyuk/PlexedPiper). The pipeline consists of three major steps.

## 1. Processing MS-GF+ data

First, [MS-GF+](https://omics.pnl.gov/software/ms-gf) provides a list of peptide identifications found in MS/MS spectra. The false discovery rate (FDR) is estimated using target-decoy search strategy. We filter the identifications by q-value and mass measurement error to maximize the number of peptide IDs with the constraint that FDR does not exceed 1%. A similar filter is then applied to the number of peptides per 1000 amino acids at the protein aggregation level with the same 1% FDR threshold. 

## 2. Processing MASIC data

Next, reporter ion intensities are taken from [MASIC](https://omics.pnl.gov/software/masic) output. They are filtered by intereference score with minimum value of 0.5. 

## 3. Relative quantification at the protein level

In this step, MS/MS identifications from step 1 and reporter ion intensities from step 2 are combined to produce a quantitative crosstab containing raw intensities at the protein level in each sample. Raw intensities are divided by intensities in the reference channel to get relative abundances. Then log2 transform is applied.

## 4. Batch correction and normalization

In multiplex design, relative abundances in different plexes can vary greatly. Therefore we need to apply batch correction. We use the following steps:

4a. Median center the quantative crosstab from step 3.

4b. Use an empirical Bayesian linear model to correct for differences due to plex.

4c. Median center the crosstab.

## 5. Uploading to Synapse.

This step uploads the final crosstab from step 4 to Synapse using the R interface.
