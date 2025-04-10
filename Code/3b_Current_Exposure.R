##### Fitting ML Models to Data
##### Daniel Dempsey

### Load in ML functions and data, and set working directory
source( 'Code/2_ML_Functions.R' )
hd_bsi <- read_csv( 'Data/hd_bsi.csv' )

home_dir <- 'Output/Model_Fits/'
mkdir( home_dir )
setwd( home_dir )

### Set data
X <- select( hd_bsi, -Exposure, -status ) %>% as.matrix
y <- hd_bsi$Exposure

### Parameter Grids
# XGBoost
param_grid_xgb <- list( max_depth = c(1, 3), eta = c(0.001, 0.01),
                        gamma = c(5, 10), subsample = c(0.5, 0.75), 
                        colsample_bytree = c(0.5, 0.75),
                        lambda = c(1, 5), objective = 'binary:logistic', 
                        nrounds = c(10, 20) )

# Random forest
param_grid_rf <- list( ntree = c(10, 50, 100), mtry = c(4, 6, 8), nodesize = c(1, 3, 5) )

### Run model fit functions
# XGBoost
tic()
fit_fun( X, y, model_mode = 1, param_grid = param_grid_xgb, outer_folds = 10, 
         inner_folds = 10, directory = 'Active_Analysis', seed = 42 )
xgb_time <- toc() # 4879.537 sec elapsed

# Random Forest
tic()
fit_fun( X, y, model_mode = 2, param_grid = param_grid_rf, outer_folds = 10, 
         inner_folds = 10, directory = 'Active_Analysis', seed = 42 )
rf_time <- toc() # 31.337 sec elapsed

# LASSO (SEEMS TO RUN INTO CONVERGENCE ISSUES!!!)
#tic()
#fit_fun( X, y, model_mode = 3, outer_folds = 10, inner_folds = 10, 
#         directory = 'Active_Analysis', seed = 42 )
#lasso_time <- toc()

