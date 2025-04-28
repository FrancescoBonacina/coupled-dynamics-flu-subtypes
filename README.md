# coupled-dynamics-flu-subtypes
This repository provides the code used for the analysis of the paper _Characterization and forecast of global influenza (sub)type dynamics_ doi: https://doi.org/10.1101/2024.08.01.24311336


**Content of the repository**

<ins> Forecast of (sub)type trajectores </ins>
- _input_forecast_subtype_compositions_2010-19.csv_ : file of data. Data points consist of country-year proportions of infections by (sub)type;
- _forecasts_M3average-M4VAR-M5hVAR.R_ : R script to forecast the (sub)type compositions one year in advance. The forecast is made via a Hierarchical Vector AutoRegressive model, where the model's parameters are optimized via Monte Carlo Gibb Sampling;
- _hVAR-GibbSampling.R_ : R script with functions to perform the Monte Carlo Gibb sampling;
- _hVAR-AnalyzeResults.R_ : R script containing functions to analyze the samples obtained from the Monte Carlo procedure.

<ins> Drivers of country-pair (sub)type distances </ins>
- _input_spatial_analysis.csv_ : file of data. List of country pairs, with country distances in terms of (sub)type compositions, temperatures, relative humidities, and epidemic synchrony. Air-traffic distances are omitted since IATA data is commercially available and restrictions apply to the redistribution of this data;
- _spatial_analysis_drivers_of_country-pair_subtype_distances.ipynb_ : Jupyter Notebook with the Python script for the linear regression analysis of the drivers of country-pair (sub)type distances.

Our code for the forecast is based on the code of _Lu, F., Zheng, Y., Cleveland, H., Burton, C. & Madigan, D. Bayesian Hierarchical
Vector Autoregressive Models for Patient-Level Predictive Modeling. PLOS ONE 13, e0208082. https://doi.org/10.1371/journal.pone.0208082 (2018)._
