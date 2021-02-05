# The PNNL Proteomics Pipeline: P3
This project provides the basic analysis and tools for carrying out end-to-end proteomic analysis for data generated at PNNL. This repository will serve as a resource upon which we can build additional tools in a reusable fashion.

## Overview

#### To contribute to this pipeline
Check out the [CONTRIBUTING.md](./CONTRIBUTING.md) document. We welcome all additional tool and methods. Then you can create a pull request and merge your code in.

#### To use this pipeline
Fork this repository into your own repository. You can rename it, make it private, and do what you need for your project. Then you can share the repository with others to serve as a code base for project development. The fork button is at the top: [fork](./img/fork.jpg).

## Pipeline overview
This pipeline links together a series of data and code repositories to ensure that any data is saved for future use. Furthermore we hope to be able to reproduce the analysis for future projects. A brief overview of each tool is described here.
[arch](./img/arch.jpg)

### Proteomics and phosphoproteomics processing
This processing relies on a series of customized `R` scripts that are in the following directories:
- [global](./global): global proteomics analysis documentation is [here](./global/README.md)
- [phospo](./phospho): phosphoproteomics analysis documentation is [here](./global/README.md)
Once this analysis is complete, the data should be stored on Synapse

### Synapse processing
Synapse is used for:
- data storage: data is backed up for free
- data processing: UUIDs and permission modifications allow for easy analysis
- data release: data can be moved to the PDC or shared

### Functional analysis


## Getting started
