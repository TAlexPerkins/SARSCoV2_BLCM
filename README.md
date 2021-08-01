SARSCoV2_BLCM
====================

This repository contains code used in the following preprint.

Perkins TA, Stephens M, Alvarez Barrios W, Cavany SM, Rulli L, Pfrender ME. (2021) **Performance of three molecular tests for SARS-CoV-2 on a university campus estimated jointly with Bayesian latent class modeling**. *medRxiv* doi:. [url](url)


### Data

The entirety of the data used for this analysis is located in `data.RData`. This is the same information contained in Table 2 in the preprint.


### Code

Code for this analysis was written and executed in the R programming language. It made use of the following R packages: BayesianTools, lubridate.

`saliva_properties.R` contains all code used in the analysis of empirical data, including the generation of associated figures.

`simulation_analysis.R` contains all code used in the analysis of simulated data.


### Outputs

The figures in the preprint are all contained here and can be generated by running the corresponding R scripts. `mcmc.RData` contains outputs from the analysis of empirical data, and `mcmc_simulated.RData` contains outputs from the analysis of simulated data.