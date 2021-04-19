Global proteomics pipeline using PlexedPiper
================
4/19/2021

The global proteomics pipeline can be run in three scripts:

- `1-create_study_design_tables.Rmd`

- `2-process_global_data.Rmd`

- `3-batch_correction_and_normalization.Rmd`

The R package PlexedPiper is required to run the first two scripts, along with a connection to the DMS. All output files, along with intermediate data, are stored in the `data/` folder.
