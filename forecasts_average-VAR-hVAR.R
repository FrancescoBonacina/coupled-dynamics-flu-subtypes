################################################################################
# forecast of subtype trajectories from 2010 to 2019
#
# forecast algorithms:
# - average
# - VAR
# - hVAR
#
# code for the VAR and hVAR predictions is based on: 
# Lu et al., Bayesian hierarchical vector autoregressive models for 
# patient-level predictive modeling, PLOS ONE (2018), 
# DOI 10.1371/journal.pone.0208082
#
# The main difference from the code of Lu et al., is that we consider a VAR
# process with the intercept term. That is, we consider:
#    y_t = nu + A_1*y_t-1 + ... + A_l * y_t-l + eps_t
# rather than:
#    y_t = A_1*y_t-1 + ... + A_l * y_t-l + eps_t
#
################################################################################

rm(list=ls())

### libraries
#library(mvnfast)
#library(png)
library(scoringRules)    # proper scoring rules for evaluating probabilistic forecasting
#library(robCompositions) # to apply log-ratio transformations (from simplex to R² and back)

### set directory for output files
th_annual_cases <- 50   # 50, 500
coda_map <- 'ilr'       # 'ilr', 'alr'
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
main_directory <- gsub('/5_forecast_of_subtype_trajectories_COPIA', '', wd)
outdir <- '4_out'
dir.create(file.path(wd, outdir), showWarnings = FALSE)
outdir_path <- paste0(wd, '/', outdir)

### other sources
source("hVAR-GibbSampling.R")
source("hVAR-AnalyzeResults.R")

### useful functions
composition_R2S <- function(v, map){
  # Function to compute coordinates of compositions from the 2D Euclidean space to the 3-parts Simplex.
  # v is a matrix of size Nx2, where N is the number of observations and the 2 columns corresponds to the
  # p,q coordinates in the Euclidean space.

  # coordinates in the Euclidean space
  p = v[,1]
  q = v[,2]
  # parts in the simplex
  if (map=='ilr'){
    x <- exp(sqrt(3/2)*p + q/sqrt(2))
    y <- exp(sqrt(2)*q)
    z <- rep(1, nrow(v))
  } else if (map=='alr'){
    x <- exp(p)
    y <- exp(q)
    z <- rep(1, nrow(v))
  } else {
    print ("ERROR: map not known. Chose a map between 'ilr' and 'alr'")
  }
  
  s <- cbind(x,y,z)
  return(s/rowSums(s))
}

my_dss <- function(y, mu, sig){
  # function to compute the Dawid-Sebastiani Score (DSS) for evaluation of 
  # multivariate probabilistic forecast.
  # Input:
  # y: vector. The observation
  # mu: vector. The point prediction, it corresponds to the center of the forecast distribution.
  # sig: matrix. The covariance matrix of the forecast distribution.
  # Output:
  # dss: scalar. The Dawid-Sebastiani Score, defined as:
  #      DSS = +log(|sig|) + (y-mu)'*sig⁻1*(y-mu)
  # formula (4) in Ziel and Berk arXiv:1910.07325
  
  # /!\ the covariance matrix must be invertible (not singular)
  
  return( log(det(sig)) + t(y-mu)%*%chol2inv(chol(sig))%*%(y-mu))
}

## ============== input data ================== ###
data_file <- paste0(main_directory, sprintf('/3_countries_with_similar_subtype_alternation/3_out/3a_df_subtype_abundances_2010-2019_th%dcases_%s-map_clustering.csv',th_annual_cases, coda_map))
#data_file <- paste0(wd, sprintf('/3a_df_subtype_abundances_2010-2019_th%dcases_%s-map_clustering.csv',th_annual_cases, coda_map))
data <- read.csv(data_file)
ncountries <- length(unique(data$country))
str(data)

## ==============   parameters of the analysis   ===============================
##  parameters related to data
start_year <- 2010
years_to_predict <- c(2017, 2018, 2019)
lag <- c(1,2)                                # lag of the autoregressive model
B <- 2                                       # nb of coordinates of the series
## parameters for MC sampling (hVAR model)
n.iter <- 20000          
n.thin <- 20
n.warmup <- 15000
n.sim <- (n.iter-n.warmup)/n.thin
nsamples <- 4*n.sim
save_diagnostic_plots <- TRUE
if (save_diagnostic_plots){dir.create(file.path(outdir_path, sprintf('MC_diagnostic_plots_th%dcases_%s-map', th_annual_cases, coda_map)), showWarnings = FALSE)}

## matrices to store results
#  AVERAGE
#cols_av <- c('country', 'group', 'year', 'p_obs', 'q_obs', 'p_hat', 'q_hat', 
#                'S_pp', 'S_pq', 'S_qp', 'S_qq', 'ES', 'DSS', 'VS05', 'VS1', 
#                'B%_hat', 'AH1%_hat', 'AH3%_hat', 'X%_hat')
cols_av <- c('country', 'group', 'year', 'p_obs', 'q_obs', 's_obs',
             'p_hat', 'q_hat', 'S_pp', 'S_pq', 'S_qp', 'S_qq',
             'ES', 'DSS', 'VS05', 'VS1',
             's_hat', 'prob_s_is_B', 'prob_s_is_H1', 'prob_s_is_H3', 'prob_s_is_X',
             'B_dom_hat', 'H1_dom_hat', 'H3_dom_hat', 'B_prob_to_be_dom', 'H1_prob_to_be_dom', 'H3_prob_to_be_dom',
             'B_neg_hat', 'H1_neg_hat', 'H3_neg_hat', 'B_prob_to_be_neg', 'H1_prob_to_be_neg', 'H3_prob_to_be_neg')

mat_av <- matrix(NA, nrow=ncountries*length(years_to_predict), ncol=length(cols_av))
#  VAR
#cols_VAR <- c('country', 'group', 'year', 'lag', 'p_obs', 'q_obs', 'p_hat', 'q_hat', 
#              'S_pp', 'S_pq', 'S_qp', 'S_qq',  'ES', 'DSS', 'VS05', 'VS1', 
#              'B%_hat', 'AH1%_hat', 'AH3%_hat', 'X%_hat',
#              'nu_p', 'nu_q', 'A1_pp', 'A1_qp', 'A1_pp', 'A1_qq', 'A2_pp', 'A2_qp', 'A2_pq', 'A2_qq')
cols_VAR <- c('country', 'group', 'year', 'lag', 'p_obs', 'q_obs', 's_obs',
              'p_hat', 'q_hat', 'S_pp', 'S_pq', 'S_qp', 'S_qq',
              'ES', 'DSS', 'VS05', 'VS1',
              's_hat', 'prob_s_is_B', 'prob_s_is_H1', 'prob_s_is_H3', 'prob_s_is_X',
              'B_dom_hat', 'H1_dom_hat', 'H3_dom_hat', 'B_prob_to_be_dom', 'H1_prob_to_be_dom', 'H3_prob_to_be_dom',
              'B_neg_hat', 'H1_neg_hat', 'H3_neg_hat', 'B_prob_to_be_neg', 'H1_prob_to_be_neg', 'H3_prob_to_be_neg',
              'nu_p', 'nu_q', 'A1_pp', 'A1_qp', 'A1_pp', 'A1_qq', 'A2_pp', 'A2_qp', 'A2_pq', 'A2_qq')
mat_VAR <- matrix(NA, nrow=ncountries*length(lag)*length(years_to_predict), ncol=length(cols_VAR))
#  hVAR
w_cols <- c()
v_cols <- c()
for (i in 1:10){
  w_cols <- c(w_cols, sprintf('w_%d',i))
  v_cols <- c(v_cols, sprintf('v_%d',i))
}
#cols_hVAR <- c('country', 'group', 'year', 'lag', 'p_obs', 'q_obs', 'p_hat', 'q_hat', 
#               'S_pp', 'S_pq', 'S_qp', 'S_qq', 'ES', 'DSS', 'VS05', 'VS1',
#               'B%_hat', 'AH1%_hat', 'AH3%_hat', 'X%_hat', 
#               w_cols, v_cols, 'L_pp', 'L_pq', 'L_qp', 'L_pp', 'R_max')
cols_hVAR <- c('country', 'group', 'year', 'lag', 'p_obs', 'q_obs', 's_obs',
               'p_hat', 'q_hat', 'S_pp', 'S_pq', 'S_qp', 'S_qq',
               'ES', 'DSS', 'VS05', 'VS1',
               's_hat', 'prob_s_is_B', 'prob_s_is_H1', 'prob_s_is_H3', 'prob_s_is_X',
               'B_dom_hat', 'H1_dom_hat', 'H3_dom_hat', 'B_prob_to_be_dom', 'H1_prob_to_be_dom', 'H3_prob_to_be_dom',
               'B_neg_hat', 'H1_neg_hat', 'H3_neg_hat', 'B_prob_to_be_neg', 'H1_prob_to_be_neg', 'H3_prob_to_be_neg',
                w_cols, v_cols, 'L_pp', 'L_pq', 'L_qp', 'L_pp', 'R_max')
mat_hVAR <- matrix(NA, nrow=ncountries*length(lag)*length(years_to_predict), ncol=length(cols_hVAR))

crows_av <- 1
crows_ar <- 1
start_time <- Sys.time()
for (gr in 1:2){
  #gr <- 1
  for (end_year in years_to_predict){
    #end_year <- 2019
    ### selection of trajectories: select years and countries in the group)
    my_group <- sprintf('group_2010%d', end_year-2001)
    data1 <- data[data[,c(my_group)]==gr & data$year>=start_year & data$year<=end_year,]
    data1$year <- data1$year - start_year + 1   # set years from 1 to T+1
    N <- length(unique(data1$country))          # nb of series
    group_name <- sprintf('group %d', gr)
    
    ### set training and test
    # list of bivariate series, one for each country
    Y <- split(data1, data1$country)
    Y.train <- lapply(Y, function(x){x[-nrow(x),c('p','q')]})     # all but the last point
    Y.test <- lapply(Y, function(x){x[nrow(x),c('p','q')]})       # the last point for each series
    TT <- end_year - start_year                                   # nb of points in the training set --> /!\ here TT=nb. obs in the training set. That is different from the definition used in Lutkepohl, where T=(nb. obs in the training set - lag)
    s_obs <- lapply(Y, function(x){ ifelse(x[nrow(x),'B.']>=0.5, 'B', ifelse(x[nrow(x),'AH1.']>=0.5, 'AH1', ifelse(x[nrow(x),'AH3.']>=0.5, 'AH3', 'X'))) })   # the observed dominance state, among {B, H1, H3, X}
    
    ### print info about the current loop
    cat('Fitting', group_name, ': ', 'year to predict =', end_year, '\n')
    
    ### AVERAGE estimator:   ----------------------------------------------------------------------------------------------------------------
    # predict the composition: point forecast and prediction uncertainties
    y_av <- lapply(Y.train, function(x){apply(x[,c('p','q')], 2, mean)})
    sig_av <- lapply(Y.train, function(x){cov(x)/TT})
    # scores to evaluate the forecast of the composition
    samples_av <- lapply(1:N, function(n){ t(mvrnorm(1000, mu=y_av[[n]], Sigma=sig_av[[n]])) })         # samples from the forecast distribution
    ES_av <- sapply(1:N, function(n){ es_sample(y=unlist(Y.test[[n]]), dat=samples_av[[n]]) })          # energy score
    #DSS_av <- sapply(1:N, function(n){ dss_sample(y=unlist(Y.test[[n]]), dat=samples_av[[n]]) })   # Dawid-Sebastiani score
    #DSS_av <- sapply(1:N, function(n){ dss(y=unlist(Y.test[[n]]), mu=y_av[[n]], Sigma=sig_av[[n]]) })   # Dawid-Sebastiani score
    DSS_av <- sapply(1:N, function(n){ tryCatch(my_dss(y=unlist(Y.test[[n]]), mu=y_av[[n]], sig=sig_av[[n]]), error=function(e) NA) })   # Dawid-Sebastiani score
    VS05_av <- sapply(1:N, function(n){ vs_sample(y=unlist(Y.test[[n]]), dat=samples_av[[n]], p=0.5) }) # variogram score of order p=0.5
    VS1_av <- sapply(1:N, function(n){ vs_sample(y=unlist(Y.test[[n]]), dat=samples_av[[n]], p=1) })    # variogram score of order p=1
    # prediction of the dominance state: point forecast and prediction uncertainties
    #y_av_simplex <- lapply(1:N, function(n){ pivotCoordInv(t(y_av[[n]])) })
    y_av_simplex <- lapply(1:N, function(n){ composition_R2S(v = t(y_av[[n]]), map = coda_map) })
    s_hat_av <- sapply(y_av_simplex, function(y){ ifelse(y[1]>=0.5, 'B', ifelse(y[2]>=0.5, 'AH1', ifelse(y[3]>=0.5, 'AH3', 'X'))) })
    #samples_av_simplex <- lapply(1:N, function(n){ pivotCoordInv(t(samples_av[[n]])) })
    samples_av_simplex <- lapply(1:N, function(n){ composition_R2S(v = t(samples_av[[n]]), map = coda_map) })
    prob_s_is_B_av <- sapply(samples_av_simplex, function(s){ sum(s[,1]>=0.5) / 1000 })                 # percentage of times the predicted dominance state is B
    prob_s_is_H1_av <- sapply(samples_av_simplex, function(s){ sum(s[,2]>=0.5) / 1000 })                # ... same for AH1, AH3 and X
    prob_s_is_H3_av <- sapply(samples_av_simplex, function(s){ sum(s[,3]>=0.5) / 1000 }) 
    prob_s_is_X_av <- sapply(samples_av_simplex, function(s){ sum(apply(s<0.5,1,prod)) / 1000 })
    #test <- lapply(1:N, function(n){prob_s_is_B_av[[n]] + prob_s_is_H1_av[[n]] + prob_s_is_H3_av[[n]] + prob_s_is_X_av[[n]] })
    # prediction of the dominance (yes or no) of B, H1, and H3
    B_dom_hat_av <- sapply(y_av_simplex, function(y){ y[1]>=0.5 })                                      # do the prediction falls in the dominance region of B?
    H1_dom_hat_av <- sapply(y_av_simplex, function(y){ y[2]>=0.5 })
    H3_dom_hat_av <- sapply(y_av_simplex, function(y){ y[3]>=0.5 })
    B_prob_to_be_dom_av <- prob_s_is_B_av                                                               # prob. that B will be dominant = prob. of the prediction to fall in the dominance region of B
    H1_prob_to_be_dom_av <- prob_s_is_H1_av
    H3_prob_to_be_dom_av <- prob_s_is_H3_av
    # prediction of the negligibility (yes or no) of B, H1, and H3
    B_neg_hat_av <- sapply(y_av_simplex, function(y){ y[1]<0.1 })                                       # do the prediction falls in the region where B is negligible (B<10% of cases)?
    H1_neg_hat_av <- sapply(y_av_simplex, function(y){ y[2]<0.1 })
    H3_neg_hat_av <- sapply(y_av_simplex, function(y){ y[3]<0.1 })
    B_prob_to_be_neg_av <- sapply(samples_av_simplex, function(s){ sum(s[,1]<0.1) / 1000 })             # percentage of times the prediction falls in the region of negligibility of B
    H1_prob_to_be_neg_av <- sapply(samples_av_simplex, function(s){ sum(s[,2]<0.1) / 1000 })
    H3_prob_to_be_neg_av <- sapply(samples_av_simplex, function(s){ sum(s[,3]<0.1) / 1000 })
    
    # fill matrix of results
    mat_av[crows_av:(crows_av+N-1),1] <- labels(Y)                     # country     
    mat_av[crows_av:(crows_av+N-1),2] <- rep(gr, N)                    # group
    mat_av[crows_av:(crows_av+N-1),3] <- rep(end_year, N)              # predicted year
    mat_av[crows_av:(crows_av+N-1),c(4,5)] <- t(matrix(unlist(Y.test), 2))         # observation: p_obs, q_obs
    mat_av[crows_av:(crows_av+N-1),6] <- t(matrix(unlist(s_obs), 1))              # the observed dominance state, among (B, AH1, AH3, X)
    mat_av[crows_av:(crows_av+N-1),c(7,8)] <- t(matrix(unlist(y_av), 2))           # prediction: p_hat, q_hat     
    mat_av[crows_av:(crows_av+N-1),c(9,10,11,12)] <- t(matrix(unlist(sig_av), 4))  # covariance matrix: S_pp, ...
    mat_av[crows_av:(crows_av+N-1),c(13,14,15,16)] <- cbind(ES_av, DSS_av, VS05_av, VS1_av)                                   # performance scores to evaluate the prediction of the composition
    mat_av[crows_av:(crows_av+N-1),17] <- s_hat_av                                                                            # the predicted dominance state, among (B, AH1, AH3, X)
    mat_av[crows_av:(crows_av+N-1),c(18,19,20,21)] <- cbind(prob_s_is_B_av, prob_s_is_H1_av, prob_s_is_H3_av, prob_s_is_X_av) # predicted probability for each dominance state (B, AH1, AH3, X)
    mat_av[crows_av:(crows_av+N-1),c(22,23,24)] <- cbind(B_dom_hat_av, H1_dom_hat_av, H3_dom_hat_av)                          # binary variables for the predicted dominance (yes/no) of each subtype
    mat_av[crows_av:(crows_av+N-1),c(25,26,27)] <- cbind(B_prob_to_be_dom_av, H1_prob_to_be_dom_av, H3_prob_to_be_dom_av)     # probabilities of dominance of each subtype
    mat_av[crows_av:(crows_av+N-1),c(28,29,30)] <- cbind(B_neg_hat_av, H1_neg_hat_av, H3_neg_hat_av)                          # binary variables for the predicted negligibility (yes/no) of each subtype
    mat_av[crows_av:(crows_av+N-1),c(31,32,33)] <- cbind(B_prob_to_be_neg_av, H1_prob_to_be_neg_av, H3_prob_to_be_neg_av)     # probabilities of negligibility of each subtype
    crows_av <- crows_av + N
    
    # try several lags for autoregressive models:
    for (l in lag){
      #l <- 1
      ### VAR:   y_t = nu + A1*y_t-1 + ... + Al * y_t-l    ----------------------------------------------------------------------------------------
      ### prepare matrices for calculations
      Z <- lapply(1:N, function(n){H <- NULL; for(i in 1:(TT-l)){H <- cbind(H, c(1, t(Y.train[[n]][(i+l-1):i,])))}; H})
      # list of N matrices, each one of shape (B*l+1) x TT-l.
      # for each matrix: the first column contains the B*l coord. necessary to predict the B coord. of point l+1
      #                  the last column contains the B*l coord. necessary to predict the B coord. of points TT
      Z.test <- lapply(1:N, function(n){c(1, t(Y.train[[n]][TT:(TT-l+1),]))})
      # list of N lists, each one containing the (B*l+1) coordinates necessary to predict the point at time TT+1 (the point in the test set)
      Y0 <- lapply(1:N, function(n){t(Y.train[[n]][(l+1):TT,])})
      # list of N matrices, each one containing B rows and TT-l columns. Columns contain coord. of points from l+1 to TT
      ZZ <- lapply(Z, tcrossprod)
      # list of the N matrices ZZ^t, of shape ((B*l+1) x (TT-l)) x ((T-l) x (B*l+1)) = (B*l+1) x (B*l+1)
      
      ### LS estimator of the VAR(l) model (formula 3.2.10, p.64, Lutkepohl (1991), Introduction to Multiple Time Series Analysis)
      B_hat <- lapply(1:N, function(n){c(Y0[[n]] %*% t(Z[[n]]) %*% chol2inv(chol(ZZ[[n]])))}) # chol2inv(...) = (ZZ^t)^-1
                                                                   
      ### prediction on the test set and prediction of uncertainties
      #   y_t = nu + A1*y_t-1 + ... + Al * y_t-l 
      yVAR.test <- lapply(1:N, function(n){matrix(B_hat[[n]],B) %*% Z.test[[n]]})
      #   covariance matrix of the predicted point corrected for small sample size
      #   look at section 3.5.4., p. 91, Lutkpohl (1991): combination of formula 3.5.17 (p. 90) and formula 3.2.23 (p. 72)
      coef_corr_sig <- (TT-l + B*l + 1)/((TT-l)*(TT-l-B*l-1))
      sig_yVAR <- lapply(1:N, function(n){Bn <- matrix(B_hat[[n]], B);
                                          coef_corr_sig * ((Y0[[n]] - Bn %*% Z[[n]]) %*% t(Y0[[n]])) })
      # prediction of the dominance state: point forecast
      #y_VAR_simplex <- lapply(1:N, function(n){ pivotCoordInv(t(yVAR.test[[n]])) })
      y_VAR_simplex <- lapply(1:N, function(n){ composition_R2S(v = t(yVAR.test[[n]]), map = coda_map) })
      s_hat_VAR <- sapply(y_VAR_simplex, function(y){ ifelse(y[1]>=0.5, 'B', ifelse(y[2]>=0.5, 'AH1', ifelse(y[3]>=0.5, 'AH3', 'X'))) })
      # prediction of the dominance (yes or no) of B, H1, and H3: point forecast
      B_dom_hat_VAR <- sapply(y_VAR_simplex, function(y){ y[1]>=0.5 })                  # do the prediction falls in the dominance region of B?
      H1_dom_hat_VAR <- sapply(y_VAR_simplex, function(y){ y[2]>=0.5 })
      H3_dom_hat_VAR <- sapply(y_VAR_simplex, function(y){ y[3]>=0.5 })
      # prediction of the negligibility (yes or no) of B, H1, and H3: point forecast
      B_neg_hat_VAR <- sapply(y_VAR_simplex, function(y){ y[1]<0.1 })                   # do the prediction falls in the region where B is negligible (B<10% of cases)?
      H1_neg_hat_VAR <- sapply(y_VAR_simplex, function(y){ y[2]<0.1 })
      H3_neg_hat_VAR <- sapply(y_VAR_simplex, function(y){ y[3]<0.1 })
      
      ### prediction uncertainties
      if (coef_corr_sig<99999){   # performance scores are defined only if the covariance matrix exists
        # scores to evaluate the forecast of the composition
        samples_VAR <- lapply(1:N, function(n){ t(mvrnorm(1000, mu=yVAR.test[[n]], Sigma=sig_yVAR[[n]])) })         # samples from the forecast distribution
        ES_VAR <- sapply(1:N, function(n){ es_sample(y=unlist(Y.test[[n]]), dat=samples_VAR[[n]]) })                # energy score
        DSS_VAR <- sapply(1:N, function(n){ tryCatch(my_dss(y=unlist(Y.test[[n]]), mu=yVAR.test[[n]], sig=sig_yVAR[[n]]), error=function(e) NA) })   # Dawid-Sebastiani score
        #DSS_VAR <- sapply(1:N, function(n){ tryCatch(dss(y=unlist(Y.test[[n]]), mu=yVAR.test[[n]], Sigma=sig_yVAR[[n]]), error=function(e) NA) })   # Dawid-Sebastiani score. If the cov. matrix is singular, the score can not be computed
        #DSS_VAR <- sapply(1:N, function(n){ dss_sample(y=unlist(Y.test[[n]]), dat=samples_VAR[[n]]) })   # Dawid-Sebastiani score. If the cov. matrix is singular, the score can not be computed
        VS05_VAR <- sapply(1:N, function(n){ vs_sample(y=unlist(Y.test[[n]]), dat=samples_VAR[[n]], p=0.5) })       # variogram score of order p=0.5
        VS1_VAR <- sapply(1:N, function(n){ vs_sample(y=unlist(Y.test[[n]]), dat=samples_VAR[[n]], p=1) })          # variogram score of order p=1
        # prediction of the dominance state: prediction uncertainties
        #samples_VAR_simplex <- lapply(1:N, function(n){ pivotCoordInv(t(samples_VAR[[n]])) })
        samples_VAR_simplex <- lapply(1:N, function(n){ composition_R2S(v = t(samples_VAR[[n]]), map = coda_map) })
        prob_s_is_B_VAR <- sapply(samples_VAR_simplex, function(s){ sum(s[,1]>=0.5) / 1000 })                 # % of times the predicted dominance state is B
        prob_s_is_H1_VAR <- sapply(samples_VAR_simplex, function(s){ sum(s[,2]>=0.5) / 1000 })                # ... same for AH1, AH3 and X
        prob_s_is_H3_VAR <- sapply(samples_VAR_simplex, function(s){ sum(s[,3]>=0.5) / 1000 }) 
        prob_s_is_X_VAR <- sapply(samples_VAR_simplex, function(s){ sum(apply(s<0.5,1,prod)) / 1000 })
        # prediction of the dominance (yes or no) of B, H1, and H3: prediction uncertainties
        B_prob_to_be_dom_VAR <- prob_s_is_B_VAR                                                               # prob. that B will be dominant = prob. of the prediction to fall in the dominance region of B
        H1_prob_to_be_dom_VAR <- prob_s_is_H1_VAR
        H3_prob_to_be_dom_VAR <- prob_s_is_H3_VAR
        # prediction of the negligibility (yes or no) of B, H1, and H3: prediction uncertainties
        B_prob_to_be_neg_VAR <- sapply(samples_VAR_simplex, function(s){ sum(s[,1]<0.1) / 1000 })             # % of times the prediction falls in the region of negligibility of B
        H1_prob_to_be_neg_VAR <- sapply(samples_VAR_simplex, function(s){ sum(s[,2]<0.1) / 1000 })
        H3_prob_to_be_neg_VAR <- sapply(samples_VAR_simplex, function(s){ sum(s[,3]<0.1) / 1000 })
      }
      
      # fill matrix of results
      mat_VAR[crows_ar:(crows_ar+N-1),1] <- labels(Y)                           # country
      mat_VAR[crows_ar:(crows_ar+N-1),2] <- rep(gr, N)                          # group
      mat_VAR[crows_ar:(crows_ar+N-1),3] <- rep(end_year, N)                    # predicted year
      mat_VAR[crows_ar:(crows_ar+N-1),4] <- rep(l,N)                            # lag
      mat_VAR[crows_ar:(crows_ar+N-1),c(5,6)] <- t(matrix(unlist(Y.test), 2))           # p_obs, q_obs
      mat_VAR[crows_ar:(crows_ar+N-1),7] <- t(matrix(unlist(s_obs), 1))                 # s_obs
      mat_VAR[crows_ar:(crows_ar+N-1),c(8,9)] <- t(matrix(unlist(yVAR.test), 2))        # p_hat, q_hat
      mat_VAR[crows_ar:(crows_ar+N-1),c(10,11,12,13)] <- t(matrix(unlist(sig_yVAR), 4)) # S_pp, S_pq, ...
      if (coef_corr_sig<99999){ mat_VAR[crows_ar:(crows_ar+N-1),c(14,15,16,17)] <- cbind(ES_VAR, DSS_VAR, VS05_VAR, VS1_VAR)}                                   # scores for prediction performance
      mat_VAR[crows_ar:(crows_ar+N-1),18] <- s_hat_VAR                                                                                                          # s_hat
      if (coef_corr_sig<99999){ mat_VAR[crows_ar:(crows_ar+N-1),c(19,20,21,22)] <- cbind(prob_s_is_B_VAR, prob_s_is_H1_VAR, prob_s_is_H3_VAR, prob_s_is_X_VAR)} # predicted probability for each dominance state (B, AH1, AH3, X)
      mat_VAR[crows_ar:(crows_ar+N-1),c(23,24,25)] <- cbind(B_dom_hat_VAR, H1_dom_hat_VAR, H3_dom_hat_VAR)                                                      # binary variables, predicted dominance (yes/no) of each subtype
      if (coef_corr_sig<99999){ mat_VAR[crows_ar:(crows_ar+N-1),c(26,27,28)] <- cbind(B_prob_to_be_dom_VAR, H1_prob_to_be_dom_VAR, H3_prob_to_be_dom_VAR)}      # probability of each subtype to be dominant
      mat_VAR[crows_ar:(crows_ar+N-1),c(29,30,31)] <- cbind(B_neg_hat_VAR, H1_neg_hat_VAR, H3_neg_hat_VAR)                                                      # binary variables, predicted begligibility (yes/no) of each subtype
      if (coef_corr_sig<99999){ mat_VAR[crows_ar:(crows_ar+N-1),c(32,33,34)] <- cbind(B_prob_to_be_neg_VAR, H1_prob_to_be_neg_VAR, H3_prob_to_be_neg_VAR)}      # probability of each subtype to be negligible
      mat_VAR[crows_ar:(crows_ar+N-1),c(35:(35+B*(B*l+1)-1))] <- t(matrix(unlist(B_hat), B*(B*l+1))) # VAR coefs: nu_p, nu_q, A1_pp, A1_qp, A1_pp, A1_qq, A2_pp, A2_qp, A2_pq, A2_qq
      
      ### hVAR(l)     --------------------------------------------------------------------------------------------------------------------------
      #set.seed(1)
      
      ### initial values
      B_hat.0 <- Reduce("+", B_hat)/N  
      #    initial value of w, computed as the mean of the B_hat matrices estimated for the N series. Single vector containing (l*B+1)*B coefficients
      vn.0 <- vector("list", N)                         
      for(n in 1:N) vn.0[[n]] <- B_hat[[n]] - B_hat.0 
      #    initial values of vn. It is a list of N vectors, each one containing l*B² elements
      # initial values for 4 different chains:
      initialA = list(w=B_hat.0, 
                     v_n=vn.0,            
                     Lambda=Posdef(n=B, ev=runif(B, 0, 0.001)),  # Lambda is the 'precision matrix' of the process
                     lambda1sq = runif(B^2*l+B, 0.1, 1000),      # hyperparameters of the L1 and L2 regularization
                     lambda2=runif(B^2*l+B, 0.1, 10),              
                     tausq =runif(B^2*l+B, 0.1, 10),             # tau², together with lambda2 they define the precision matrix of vector w
                     U_v = runif(B^2*l+B, 0.1, 1),               # theta_v, the diagonal elements of the 'precision matrix' Theta_v,
                     alpha = runif(B^2*l+B, 0.1, 1)              # additional term to make the MC chain effective in the research
      )
      initialA$v_n.star <- lapply(1:N, function(n){initialA$alpha*initialA$v_n[[n]]}) # add the initial condition for v*_n = alpha*v_n
      initialA$omega_v <- initialA$U_v/(initialA$alpha^2)                             # add the initial condition for omega_v, that is the precision matrix of v*_n
      
      initialB = list(w=B_hat.0, v_n=vn.0, Lambda=Posdef(n=B, ev=runif(B, 0, 0.001)), lambda1sq = runif(B^2*l+B, 0.1, 1000),
                      lambda2=runif(B^2*l+B, 0.1, 10), tausq =runif(B^2*l+B, 0.1, 10), U_v = runif(B^2*l+B, 0.1, 1), alpha = runif(B^2*l+B, 0.1, 1))
      initialB$v_n.star <- lapply(1:N, function(n){initialB$alpha*initialB$v_n[[n]]})
      initialB$omega_v <- initialB$U_v/(initialB$alpha^2)

      initialC = list(w=B_hat.0, v_n=vn.0, Lambda=Posdef(n=B, ev=runif(B, 0, 0.001)), lambda1sq = runif(B^2*l+B, 0.1, 1000),
                      lambda2=runif(B^2*l+B, 0.1, 10), tausq =runif(B^2*l+B, 0.1, 10), U_v = runif(B^2*l+B, 0.1, 1), alpha = runif(B^2*l+B, 0.1, 1))
      initialC$v_n.star <- lapply(1:N, function(n){initialC$alpha*initialC$v_n[[n]]})
      initialC$omega_v <- initialC$U_v/(initialC$alpha^2) 
      
      initialD = list(w=B_hat.0, v_n=vn.0, Lambda=Posdef(n=B, ev=runif(B, 0, 0.001)), lambda1sq = runif(B^2*l+B, 0.1, 1000),
                      lambda2=runif(B^2*l+B, 0.1, 10), tausq =runif(B^2*l+B, 0.1, 10), U_v = runif(B^2*l+B, 0.1, 1), alpha = runif(B^2*l+B, 0.1, 1))
      initialD$v_n.star <- lapply(1:N, function(n){initialD$alpha*initialD$v_n[[n]]})
      initialD$omega_v <- initialD$U_v/(initialD$alpha^2) 
      
      ### run the MCMC sampling (Gibbs sampler with parameter expansion)
      simsA <- MSGibbs(H=Z, HH=ZZ, eta0=Y0, w_MLE=B_hat, TT=TT, B=B, l=l, N=N, seed=1, numIter=n.iter, thin=n.thin, warmup=n.warmup, initial=initialA)
      simsB <- MSGibbs(H=Z, HH=ZZ, eta0=Y0, w_MLE=B_hat, TT=TT, B=B, l=l, N=N, seed=50001, numIter=n.iter, thin=n.thin, warmup=n.warmup, initial=initialB)
      simsC <- MSGibbs(H=Z, HH=ZZ, eta0=Y0, w_MLE=B_hat, TT=TT, B=B, l=l, N=N, seed=100001, numIter=n.iter, thin=n.thin, warmup=n.warmup, initial=initialC)
      simsD <- MSGibbs(H=Z, HH=ZZ, eta0=Y0, w_MLE=B_hat, TT=TT, B=B, l=l, N=N, seed=150001, numIter=n.iter, thin=n.thin, warmup=n.warmup, initial=initialD)
      #image_name <- sprintf('%s_VAR-hVAR_no_preprocessing_fit%d-%d_%s_lag%d.RData', nb_code, start_year, end_year, group_name, l) 
      #save.image(image_name)
      
      ### posterior inference
      #load(image_name)
      # (1) w
      filename <- sprintf("%s/MC_diagnostic_plots_th%dcases_%s-map/w_trace_hVAR(%d)_%d-%d_%s.pdf", outdir_path, th_annual_cases, coda_map, l, start_year, end_year-2000, group_name)
      sw <- w_inference(chains=list(simsA, simsB, simsC, simsD), B=B, l=l, n.sim=n.sim, warmup=0, thin=1, filename=filename, plot=save_diagnostic_plots, d=3) # d=nb of decimals to use in calculations
      w_hat_Bayes <- sw$w.mode
      # (2) v
      filename <- sprintf("%s/MC_diagnostic_plots_th%dcases_%s-map/v_trace_hVAR(%d)_%d-%d_%s.pdf", outdir_path, th_annual_cases, coda_map, l, start_year, end_year-2000, group_name)
      sv <- v_inference(chains=list(simsA, simsB, simsC, simsD), B=B, l=l, N=N, n.sim=n.sim, warmup=0, thin=1, filename=filename, plot=save_diagnostic_plots, d=3)
      v_hat_Bayes <- sv$v.mode
      # (3) Lambda
      filename <- sprintf("%s/MC_diagnostic_plots_th%dcases_%s-map/Lambda_trace_hVAR(%d)_%d-%d_%s.pdf", outdir_path, th_annual_cases, coda_map, l, start_year, end_year-2000, group_name)
      sL <- L_inference(chains=list(simsA, simsB, simsC, simsD), B=B, l=l, n.sim=n.sim, warmup=0, thin=1, filename=filename, plot=save_diagnostic_plots, d=3)  
      L_hat_Bayes <- sL$L.mode
      # (4) omega_v
      filename <- sprintf("%s/MC_diagnostic_plots_th%dcases_%s-map/Omega_v_trace_hVAR(%d)_%d-%d_%s.pdf", outdir_path, th_annual_cases, coda_map, l, start_year, end_year-2000, group_name)
      somega_v <- omega_v_inference(chains=list(simsA, simsB, simsC, simsD), B=B, l=l, n.sim=n.sim, warmup=0, thin=1, filename=filename, plot=save_diagnostic_plots, d=3) 
      sdv <- sqrt(1/somega_v$somega_v)
      #round(matrix(apply(sdv, 3, function(x){find.mode(x, d=3)}), B),3)   # posterior mode of SD of v_n
      # (5) R_hat - goodness of MCMC exploration
      mon_w <- sw$mon.w                      
      mon_v <- sv$mon.v                      
      mon_L <- sL$mon.L                      
      mon_somega_v <- somega_v$mon.omega_v   
      R_hat <- c(mon_w$Rhat, mon_v$Rhat, mon_L$Rhat, mon_somega_v$Rhat)
      R_hat_max <- max(R_hat)
      if (save_diagnostic_plots){
        filename <- sprintf("%s/MC_diagnostic_plots_th%dcases_%s-map/Rhat_hVAR(%d)_%d-%d_%s.pdf", outdir_path, th_annual_cases, coda_map, l, start_year, end_year-2000, group_name)
        pdf(filename)
        hist(R_hat, breaks=30)
        dev.off()
      }

      ### prediction
      # (1). Bayesian predictive distribution
      nsamples <- length(sw$sw)/(B^2*l+B)                    # nb samples = nb chains*nb samples per chain = total nb of parameter estimated / nb params of the model
      y_hat_post_array <- array(NA, dim=c(nsamples, B, N))
      y_hat_post <- matrix(NA, B, N)                         # posterior mean of y_hat
      for(i in 1:N){   
        temp <- sw$sw + sv$sv[,,((i-1)*(B^2*l+B)+1):(i*(B^2*l+B))]    
        # computation of w_n = w + v_n (we are calculating nchain*nsamples=250*4 values)
        linear <- apply(temp, c(1,2), function(x){matrix(x, B)%*% Z.test[[i]]})
        # computation of y_T+1 = w * y_T (250*4 values)
        set.seed(123) # set seed in order to add the same error to the different trajectories (to compare trajectories)                                 
        error <- apply(sL$sL, c(1,2), function(x){Sigma <- chol2inv(chol(matrix(x, B))); rmvn(n=1, mu=rep(0, B), sigma=chol(Sigma), isChol=TRUE, ncores = 8)})
        # generate random errors ~ MVN(0, Lambda^-1) 
        pred <- linear + error                        
        # 250*4 predictions = expected value + err
        
        y_hat_post_array[,1,i] <- pred[1,,]   # 250*4 estimations for the first coord
        y_hat_post_array[,2,i] <- pred[2,,]   # 250*4 estimations for the second coord
        y_hat_post[,i] <- apply(pred, 1, mean)                               # posterior mean of y_hat (mean over the rows, for each coord. separately)
        #y_hat_post[,i] <- apply(pred, 1, function(x){find.mode(x, d=3)})    # posterior mode of y_hat
      }
      # covariance matrix of posterior predictions
      sig_yhVAR <- sapply(1:N, function(n){c(cov(y_hat_post_array[,,n]))})
      
      # scores to evaluate the forecast of the composition
      ES_hVAR <- sapply(1:N, function(n){ es_sample(y=unlist(Y.test[[n]]), dat=t(y_hat_post_array[,,n])) })                # energy score
      DSS_hVAR <- sapply(1:N, function(n){ tryCatch(my_dss(y=unlist(Y.test[[n]]), mu=matrix(y_hat_post[,n], nrow=B), sig=matrix(sig_yhVAR[,n], nrow=B)), error=function(e) NA) })   # Dawid-Sebastiani score
      #DSS_hVAR <- sapply(1:N, function(n){ dss(y=unlist(Y.test[[n]]), mu=matrix(y_hat_post[,n], nrow=B), Sigma=matrix(sig_yhVAR[,n], nrow=B)) })   # Dawid-Sebastiani score
      #DSS_hVAR <- sapply(1:N, function(n){ dss_sample(y=unlist(Y.test[[n]]), dat=t(y_hat_post_array[,,n])) })   # Dawid-Sebastiani score
      VS05_hVAR <- sapply(1:N, function(n){ vs_sample(y=unlist(Y.test[[n]]), dat=t(y_hat_post_array[,,n]), p=0.5) })       # variogram score of order p=0.5
      VS1_hVAR <- sapply(1:N, function(n){ vs_sample(y=unlist(Y.test[[n]]), dat=t(y_hat_post_array[,,n]), p=1) })          # variogram score of order p=1
      # prediction of the dominance state: point forecast and prediction uncertainties
      #y_hVAR_simplex <- lapply(1:N, function(n){ pivotCoordInv(t(y_hat_post[[n]])) })
      y_hVAR_simplex <- lapply(1:N, function(n){ composition_R2S(v = t(y_hat_post[,n]), map = coda_map) })
      s_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ ifelse(y[1]>=0.5, 'B', ifelse(y[2]>=0.5, 'AH1', ifelse(y[3]>=0.5, 'AH3', 'X'))) })
      #samples_hVAR_simplex <- lapply(1:N, function(n){ pivotCoordInv(y_hat_post_array[,,n]) })
      samples_hVAR_simplex <- lapply(1:N, function(n){ composition_R2S(v = y_hat_post_array[,,n], map = coda_map) })
      prob_s_is_B_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(s[,1]>=0.5) / 1000 })                 # % of times the predicted dominance state is B
      prob_s_is_H1_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(s[,2]>=0.5) / 1000 })                # ... same for AH1, AH3 and X
      prob_s_is_H3_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(s[,3]>=0.5) / 1000 }) 
      prob_s_is_X_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(apply(s<0.5,1,prod)) / 1000 })
      # prediction of the dominance (yes or no) of B, H1, and H3: point forecast and prediction uncertainties
      B_dom_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ y[1]>=0.5 })                                      # do the prediction falls in the dominance region of B?
      H1_dom_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ y[2]>=0.5 })
      H3_dom_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ y[3]>=0.5 })
      B_prob_to_be_dom_hVAR <- prob_s_is_B_hVAR                                                               # prob. that B will be dominant = prob. of the prediction to fall in the dominance region of B
      H1_prob_to_be_dom_hVAR <- prob_s_is_H1_hVAR
      H3_prob_to_be_dom_hVAR <- prob_s_is_H3_hVAR
      # prediction of the negligibility (yes or no) of B, H1, and H3: point forecast and prediction uncertainties
      B_neg_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ y[1]<0.1 })                                       # do the prediction falls in the region where B is negligible (B<10% of cases)?
      H1_neg_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ y[2]<0.1 })
      H3_neg_hat_hVAR <- sapply(y_hVAR_simplex, function(y){ y[3]<0.1 })
      B_prob_to_be_neg_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(s[,1]<0.1) / 1000 })             # % of times the prediction falls in the region of negligibility of B
      H1_prob_to_be_neg_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(s[,2]<0.1) / 1000 })
      H3_prob_to_be_neg_hVAR <- sapply(samples_hVAR_simplex, function(s){ sum(s[,3]<0.1) / 1000 })

      # fill matrix of results
      mat_hVAR[crows_ar:(crows_ar+N-1),1] <- labels(Y)                          # country
      mat_hVAR[crows_ar:(crows_ar+N-1),2] <- rep(gr, N)                         # group
      mat_hVAR[crows_ar:(crows_ar+N-1),3] <- rep(end_year, N)                   # predicted year
      mat_hVAR[crows_ar:(crows_ar+N-1),4] <- rep(l,N)                           # lag
      mat_hVAR[crows_ar:(crows_ar+N-1),c(5,6)] <- t(matrix(unlist(Y.test), 2))  # p_obs, q_obs
      mat_hVAR[crows_ar:(crows_ar+N-1),7] <- t(matrix(unlist(s_obs), 1))        # s_obs
      mat_hVAR[crows_ar:(crows_ar+N-1),c(8,9)] <- t(y_hat_post)                 # p_hat, q_hat
      mat_hVAR[crows_ar:(crows_ar+N-1),c(10,11,12,13)] <- t(sig_yhVAR)          # S_pp, S_pq, ...
      mat_hVAR[crows_ar:(crows_ar+N-1),c(14,15,16,17)] <- cbind(ES_hVAR, DSS_hVAR, VS05_hVAR, VS1_hVAR)  # scores for prediction performance
      mat_hVAR[crows_ar:(crows_ar+N-1),18] <- s_hat_hVAR                        # s_hat
      mat_hVAR[crows_ar:(crows_ar+N-1),c(19,20,21,22)] <- cbind(prob_s_is_B_hVAR, prob_s_is_H1_hVAR, prob_s_is_H3_hVAR, prob_s_is_X_hVAR) # predicted probability for each dominance state (B, AH1, AH3, X)
      mat_hVAR[crows_ar:(crows_ar+N-1),c(23,24,25)] <- cbind(B_dom_hat_hVAR, H1_dom_hat_hVAR, H3_dom_hat_hVAR)                            # binary variable, for the predicted dominance (yes/no) of each subtype
      mat_hVAR[crows_ar:(crows_ar+N-1),c(26,27,28)] <- cbind(B_prob_to_be_dom_hVAR, H1_prob_to_be_dom_hVAR, H3_prob_to_be_dom_hVAR)       # probability of each subtype to be negligible
      mat_hVAR[crows_ar:(crows_ar+N-1),c(29,30,31)] <- cbind(B_neg_hat_hVAR, H1_neg_hat_hVAR, H3_neg_hat_hVAR)                            # binary variable, for the predicted negligible (yes/no) of each subtype
      mat_hVAR[crows_ar:(crows_ar+N-1),c(32,33,34)] <- cbind(B_prob_to_be_neg_hVAR, H1_prob_to_be_neg_hVAR, H3_prob_to_be_neg_hVAR)       # probability of each subtype to be negligible
      mat_hVAR[crows_ar:(crows_ar+N-1),c(35:(35+length(w_hat_Bayes)-1))] <- t(replicate(N, w_hat_Bayes))                # group level coefs: W
      mat_hVAR[crows_ar:(crows_ar+N-1),c(45:(45+length(w_hat_Bayes)-1))] <- t(matrix(v_hat_Bayes, length(w_hat_Bayes))) # country level coefs: V_c
      mat_hVAR[crows_ar:(crows_ar+N-1),c(55,56,57,58)] <- t(replicate(N, L_hat_Bayes))          # Lambda, the cov. matrix of the process. (Lambda⁻1 is the precision matrix)
      mat_hVAR[crows_ar:(crows_ar+N-1),59] <- rep(R_hat_max, N)                                 # Rmax (diagnosis of MC convergence)
      crows_ar <- crows_ar + N
    }
  }
}
end_time <- Sys.time()
end_time - start_time

### save dataframes with results for each prediction method
# AVERAGE
df_average<- data.frame(mat_av)
colnames(df_average) <- cols_av
name_dfaverage <- sprintf('%s/4b_df_th%dcases_%s-map_pred-year-%s_AVERAGE_%diter.csv', outdir_path, th_annual_cases, coda_map, paste(years_to_predict, collapse='-'), n.iter)
write.table(df_average, name_dfaverage, sep = ',', dec = '.', row.names = F, col.names = T)
# VAR
df_VAR <- data.frame(mat_VAR)
colnames(df_VAR) <- cols_VAR
name_dfVAR <- sprintf('%s/4b_df_th%dcases_%s-map_pred-year-%s_VAR_%diter.csv', outdir_path, th_annual_cases, coda_map, paste(years_to_predict, collapse='-'), n.iter)
write.table(df_VAR, name_dfVAR, sep = ',', dec = '.', row.names = F, col.names = T)
# hVAR
df_hVAR <- data.frame(mat_hVAR)
colnames(df_hVAR) <- cols_hVAR
name_dfhVAR <- sprintf('%s/4b_df_th%dcases_%s-map_pred-year-%s_hVAR_%diter.csv', outdir_path, th_annual_cases, coda_map, paste(years_to_predict, collapse='-'), n.iter)
write.table(df_hVAR, name_dfhVAR, sep = ',', dec = '.', row.names = F, col.names = T)
