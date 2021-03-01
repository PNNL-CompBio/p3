# Global Proteomics analysis

The global proteomics pipeline consists of three major steps:

1. Processing and filtering of MS/MS identifications and reporter ion intensities.

2. Relative quantification at the protein level.

3. Batch correction and normalization.

4. Linking with Synapse.

The primary inputs are MS/MS identifications from MS-GF+, reporter ion intensities from MASIC, and study design tables provided by the user. The output is a quantitative crosstab giving the log-fold change in abundance relative to a reference sample. 

Steps 1 and 2 rely on the package PlexedPiper. The entire pipeline requires a connection to PNNL's Data Management System (DMS).
