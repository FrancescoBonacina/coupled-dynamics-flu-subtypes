# coupled-dynamics-flu-subtypes
This repository provides the code used for the analysis of the paper _Characterization and forecast of global influenza (sub)type dynamics_

**Content of the repository**
- _forecasts_average-VAR-hVAR.R_ : R script to forecast the subtype compositions one year in advance. The forecast is made via a Hierarchical Vector AutoRegressive model, where the model's parameters are optimized via Monte Carlo Gibb Sampling;
- _hVAR-GibbSampling.R_ : R script with functions to perform the Monte Carlo Gibb sampling;
- _hVAR-AnalyzeResults.R_ : R script containing functions to analyze the samples obtained from the Monte Carlo procedure;
- _data_subtype_abundances_2010-2019_th50cases_ilr-map_clustering.csv_ : file of data. Data points consist of country-year proportions of infections by (sub)type.

Our code is based on the code of _Lu, F., Zheng, Y., Cleveland, H., Burton, C. & Madigan, D. Bayesian Hierarchical
Vector Autoregressive Models for Patient-Level Predictive Modeling. PLOS ONE 13, e0208082. https://doi.org/10.1371/journal.pone.0208082 (2018)._
