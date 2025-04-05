##### Fitting ML Models to Data
##### Daniel Dempsey

### Load in ML functions and data, and set working directory
source( 'Code/2a_ML_Functions.R' )
hd_bsi <- read_csv( 'Data/hd_bsi.csv' )

home_dir <- 'Output/Model_Fits/'
mkdir( home_dir )
setwd( home_dir )

##### Analysis of Prior Exposure vs Never Exposed
n_keep <- !as.logical( hd_bsi$Exposure )
X1 <- select( hd_bsi, -Exposure, -status ) %>% filter( n_keep ) %>% as.matrix
y1 <- ifelse( hd_bsi$status == 'Prior Exposure', 1, 0 )[n_keep]

### Parameter Grids
# XGBoost
param_grid_xgb1 <- list( max_depth = c(1, 3, 5), eta = c(0.005, 0.01, 0.05, 0.1),
                         gamma = c(5, 10, 15), subsample = c(0.5, 0.75, 1), 
                         colsample_bytree = c(0.5, 0.75, 1),
                         lambda = c(1, 5), objective = 'binary:logistic', 
                         nrounds = c(10, 20, 50, 100) )

# Random forest
param_grid_rf1 <- list( ntree = c(10, 50, 100), mtry = c(4, 6, 8), nodesize = c(1, 3, 5) )

### Run model fit functions
# XGBoost
tic()
fit_fun( X1, y1, model_mode = 1, param_grid = param_grid_xgb1, outer_folds = 10, 
         inner_folds = 10, directory = 'Prior_Exposure_Analysis', seed = 42 )
xgb_time1 <- toc() # 4879.537 sec elapsed

# Random Forest
tic()
fit_fun( X1, y1, model_mode = 2, param_grid = param_grid_rf1, outer_folds = 10, 
         inner_folds = 10, directory = 'Prior_Exposure_Analysis', seed = 42 )
rf_time1 <- toc() # 40.663 sec elapsed

##### Analysis of Active Infection vs Non-Active Infection
X2 <- select( hd_bsi, -Exposure, -status ) %>% as.matrix
y2 <- hd_bsi$Exposure

### Parameter Grids
# XGBoost
param_grid_xgb2 <- list( max_depth = c(1, 3), eta = c(0.001, 0.01),
                         gamma = c(5, 10), subsample = c(0.5, 0.75), 
                         colsample_bytree = c(0.5, 0.75),
                         lambda = c(1, 5), objective = 'binary:logistic', 
                         nrounds = c(10, 20) )

# Random forest
param_grid_rf2 <- list( ntree = c(10, 50, 100), mtry = c(4, 6, 8), nodesize = c(1, 3, 5) )

### Run model fit functions
# XGBoost
tic()
fit_fun( X2, y2, model_mode = 1, param_grid = param_grid_xgb2, outer_folds = 10, 
         inner_folds = 10, directory = 'Active_Analysis', seed = 42 )
xgb_time2 <- toc() # 4879.537 sec elapsed

# Random Forest
tic()
fit_fun( X2, y2, model_mode = 2, param_grid = param_grid_rf2, outer_folds = 10, 
         inner_folds = 10, directory = 'Active_Analysis', seed = 42 )
rf_time2 <- toc() # 31.337 sec elapsed

