# ExplainingDelta

This repository provides R code and data sets for the experiments underlying

> Evert, S., Proisl, T., Jannidis, F., Reger, I., Pielström, S., Schöch, C., and Vitt, T. (2017). [Understanding and explaining Delta measures for authorship attribution](https://doi.org/10.1093/llc/fqx023). _Digital Scholarship in the Humanities_, **22**(suppl_2): ii4--ii16.

as well as many further investigations that did not fit into the paper.

## Prerequisites

- a recent version of [**R**](https://www.r-project.org/) 

- [RStudio](https://posit.co/download/rstudio-desktop/) is strongly recommended for exceuting the RMarkdown notebooks

- these add-on packages

  - ape
  - cluster (should be pre-installed with R)
  - colorspace
  - corpora
  - ggplot2
  - knitr
  - MASS (should be pre-installed with R)
  - pvclust
  - matrixStats
  - tsne
  - wordspace

- you can install all necessary packages with the R command 

  ```R
  install.packages(c("ape", "cluster", "colorspace", "corpora", "ggplot2", "knitr", "MASS", "pcvlust", "matrixStats", "tsne", "wordspace"))
  ```



## Overview of notebooks

- `delta_tools.Rmd`: utility functions for experiments (automatically converted to `delta_tools.R` and imported by other notebooks)
- `dsh_2016.Rmd`: replication and plots for Evert et al. (2017)
- `delta_theory.Rmd`: mathematical notation and some background on Delta measures
- `delta_scaling.Rmd`: experiments on feature selection and feature scaling
- `delta_normalization.Rmd`: experiments on vector normalization
- `delta_clustering.Rmd`: experiments with different clustering algorithms
- `delta_understanding.Rmd`: understanding Burrows Delta and its variants
- `delta_significance.Rmd`: measuring the statistical significance of Delta-based authorshop attribution

