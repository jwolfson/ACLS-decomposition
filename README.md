# ACLS-decomposition

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2574138.svg)](https://doi.org/10.5281/zenodo.2574138)

## Description

Code supporting "Decomposing the Effects of Physical Activity and Cardiorespiratory Fitness on Mortality".

## Instructions for use

1. Clone the repository, or download the files `fourway-decomp-functions.R`, `tables-and-figures.Rmd`, and the `data` subfolder onto your computer.
2. Open the R Markdown file `tables-and-figures.Rmd`.
3. Change the `dat.folder` and `code.folder` arguments at the top to correspond to the appropriate directories on your computer.
4. Knit the file.

Note that the first compilation of the document will take several minutes as the decomposition analysis is being run (bootstrapping the standard errors is the time-consuming aspect). After that, the results are cached so that future compilations will be much faster.
