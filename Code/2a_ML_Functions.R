##### XGBoost and Random Forest Fitting Functions
##### Daniel Dempsey

### Read in libraries
library( xgboost ) # Fit xgboost model
library( randomForest ) # Fit random forest model
library( SHAPforxgboost ) # For shapley value plots
library( treeshap ) # To compute shapley values
library( readr ) # Load in csv files
library( dplyr ) # Data frame manipulation tools
library( tidyr ) # Data frame manipulation tools
library( pROC ) # Computing ROC curves and AUC
library( Gmedian ) # Computing the geometric median 
library( tictoc ) # For timing code
library( caret ) # Cross validation and random forest tools
library( mice ) # Data imputation tools
library( data.table ) # Used when computing Shapley values
library( ggplot2 ) # Plotting functionality
mycol <- c( '#56B4E9', '#E69F00' )

### Function to make directories if they don't already exist
mkdir <- function( x, safe = TRUE, safe_dir = 'HD_BSI' ) {
  if ( safe ) {
    where_now <- sub( ".*/", "", getwd() )
    if ( where_now != safe_dir ) {
      stop( paste0("You are trying to create a folder where you don't want to! Move to the ", safe_dir, " directory.") )
    }
  }
  dir_split <- strsplit( x, '/' )[[1]]
  current_dir <- ''
  for ( i in 1:length(dir_split) ) {
    current_dir <- paste0( current_dir, dir_split[i], sep = '/' )
    if( !dir.exists( current_dir ) ) {
      dir.create( current_dir )
    }
  }
}

### Loss function for cross validation (with hacky fix to avoid bugs)
log_loss <- function( pred, true ) {
  pred <- ifelse( pred == 0, pred + .Machine$double.eps, pred )
  pred <- ifelse( pred == 1, pred - .Machine$double.eps, pred )
  -sum(true*log(pred) + (1-true)*log(1-pred))   
}

### Alternate version of randomForest.unify that fixes a bug that occurs for classification problems
randomForest2.unify <- function (rf_model, data) {
  if (!inherits(rf_model, "randomForest")) {
    stop("Object rf_model was not of class \"randomForest\"")
  }
  if (any(attr(rf_model$terms, "dataClasses") != "numeric")) {
    stop("Models built on data with categorical features are not supported - please encode them before training.")
  }
  n <- rf_model$ntree
  ret <- data.table()
  x <- lapply(1:n, function(tree) {
    tree_data <- as.data.table(randomForest::getTree(rf_model, 
                                                     k = tree, labelVar = TRUE))
    tree_data[, c("left daughter", "right daughter", "split var", 
                  "split point", "prediction")]
  })
  times_vec <- sapply(x, nrow)
  y <- rbindlist(x)
  y$prediction <- as.numeric( y$prediction ) # <-------- this is the fix
  y[, `:=`(Tree, rep(0:(n - 1), times = times_vec))]
  y[, `:=`(Node, unlist(lapply(times_vec, function(x) 0:(x - 1))))]
  setnames(y, c("Yes", "No", "Feature", "Split", "Prediction", 
                "Tree", "Node"))
  y[, `:=`(Feature, as.character(Feature))]
  y[, `:=`(Yes, Yes - 1)]
  y[, `:=`(No, No - 1)]
  y[y$Yes < 0, "Yes"] <- NA
  y[y$No < 0, "No"] <- NA
  y[, `:=`(Missing, NA)]
  y[, `:=`(Missing, as.integer(Missing))]
  ID <- paste0(y$Node, "-", y$Tree)
  y$Yes <- match(paste0(y$Yes, "-", y$Tree), ID)
  y$No <- match(paste0(y$No, "-", y$Tree), ID)
  y$Cover <- 0
  y$Decision.type <- factor(x = rep("<=", times = nrow(y)), 
                            levels = c("<=", "<"))
  y[is.na(Feature), `:=`(Decision.type, NA)]
  y[!is.na(Feature), `:=`(Prediction, NA)]
  y[is.na(Feature), `:=`(Prediction, Prediction/n)]
  setcolorder(y, c("Tree", "Node", "Feature", "Decision.type", 
                   "Split", "Yes", "No", "Missing", "Prediction", "Cover"))
  ret <- list(model = as.data.frame(y), data = as.data.frame(data))
  class(ret) <- "model_unified"
  attr(ret, "missing_support") <- FALSE
  attr(ret, "model") <- "randomForest"
  return(set_reference_dataset(ret, as.data.frame(data)))
}

### Function for missing data imputation via partial mean matching
# A lot of this is taken from the mice package internals
dynamic_pmm_within <- function( ci, dat, ti, donors = 5, ridge = 1e-05 , ... ) {
  
  # Prepare data
  y <- dat[, ci]
  y_obs <- !is.na( y )
  x_all <- cbind( 1, dat[, -ci] )
  x <- x_all[, colSums(is.na(x_all)) == 0]
  
  # Ridge regression on training set only
  x_train <- x[-ti, , drop = FALSE]
  y_train <- y[-ti]
  y_obs_train <- y_obs[-ti]
  
  
  const_flag <- apply( x_train, 2, function(x) !(length(unique(x)) == 1) )
  const_flag[1] <- TRUE
  p <- estimice( x_train[y_obs_train, const_flag, drop = FALSE], y_train[y_obs_train], ... )
  sigma.star <- sqrt( sum((p$r)^2)/rchisq(1, p$df) )
  beta.star <- p$c + ( t(chol(p$v) ) %*% rnorm( ncol(x[, const_flag])) ) * sigma.star
  parm <- list( p$c, beta.star )
  
  # Matching on both sets
  yhatobs <- x[y_obs, const_flag, drop = FALSE] %*% parm[[1]]
  yhatmis <- x[!y_obs, const_flag, drop = FALSE] %*% parm[[2]]
  idx <- matchindex( yhatobs, yhatmis, donors )
  
  # Return imputed values
  y[!y_obs] <- y[y_obs][idx]
  y
  
}

### Function that applies imputation across entire dataset automatically
dynamic_pmm <- function( X, ti, ... ) {
  
  X_imp <- X
  mis_cols <- which( apply(X, 2, function(x) { any(is.na(x)) }) ) %>% sample
  
  for ( i in 1:length(mis_cols) ) {
    X_imp[, mis_cols[i]] <- dynamic_pmm_within( mis_cols[i], X_imp, ti, ... )
  }
  
  X_imp
  
}

### Function to parse the shapley value results into an easier format
# This was mostly written by ChatGPT-4o
shap_convert <- function( shap_res ) {
  
  # Extract variable names and row indices
  varnames <- colnames(shap_res[[1]])
  row_indices <- seq_len(nrow(shap_res[[1]]))
  
  # Create a data frame with all combinations of varname and row_index
  base_df <- expand.grid(
    row_index = row_indices,
    varname = varnames
  )
  
  # Reorder to match the original layout (varname first, row_index second)
  base_df <- base_df[, c("varname", "row_index")]
  
  # For each dataframe in the list, extract values in the correct order
  get_values <- function(df) {
    unlist(df[, varnames], use.names = FALSE)
  }
  
  # Combine values from all dataframes
  value_columns <- lapply(shap_res, get_values)
  
  # Bind them into the final dataframe
  cbind(base_df, setNames(as.data.frame(value_columns), paste0("df", seq_along(shap_res))))
  
}


### Fitting function
# Model mode = 1 corresponds to XGBoost, anything else is Random Forest
fit_fun <- function( X, y, param_grid, model_mode = 1, outer_folds = nrow(X), 
                     inner_folds = 10, directory = 'Unlabelled', seed = 1 ) {
  
  set.seed( seed )
  
  # set working directory
  mkdir( directory, safe_dir = 'Model_Fits' )
  setwd( directory )
  sub_directory <- ifelse( model_mode == 1, 'XGBoost', 'RandomForest' )
  mkdir( sub_directory, safe_dir = directory )
  setwd( sub_directory )
  
  # Prepare parameter grid for cross-validation
  nrounds_ind <- which( names(param_grid) == 'nrounds' )
  all_param_grid <- expand.grid( param_grid )
  all_param_grid_split <- split( all_param_grid, 1:nrow(all_param_grid) )
  
  # Create a fully imputed set purely for Shapley calculation
  # Necessary to produce a consistent set of values (not used in training/validation so no leakage)
  X_all_imp <- complete( mice(X, m = 1, printFlag = FALSE) )
  
  # Begin outer loop
  outer_fold_inds <- createFolds( as.factor(y), k = outer_folds )
  outer_preds_train_all <- outer_preds_test_all <- outer_shaps_all <- best_params_all <- list()
  all_log_loss <- rep( 0, outer_folds )
  for ( i in 1:outer_folds ) {
    
    # Create outer sets
    X_outer_train <- X[-outer_fold_inds[[i]], , drop = FALSE]
    y_outer_train <- y[-outer_fold_inds[[i]]]
    
    X_outer_test <- X[outer_fold_inds[[i]], , drop = FALSE]
    y_outer_test <- y[outer_fold_inds[[i]]]
    
    # Begin inner loop
    inner_fold_inds <- createFolds( as.factor(y_outer_train), k = inner_folds )
    all_test_preds <- matrix( 0, nrow = inner_folds, ncol = nrow(all_param_grid) )
    for ( j in 1:inner_folds ) {
      
      # Dynamic imputation
      X_imp_inner <- dynamic_pmm( X_outer_train, inner_fold_inds[[j]] )
      
      # Create inner sets
      X_inner_train <- X_imp_inner[-inner_fold_inds[[j]], , drop = FALSE]
      y_inner_train <- y_outer_train[-inner_fold_inds[[j]]]
      
      X_inner_test <- X_imp_inner[inner_fold_inds[[j]], , drop = FALSE]
      y_inner_test <- y_outer_train[inner_fold_inds[[j]]]
      
      # Observation weighting based on imbalance
      weight_cv <- rep( 1, length(y_inner_train) )
      weight_cv[y_inner_train == 1] <- sum( 1 - y_inner_train ) / sum( y_inner_train )
      
      # Fit the desired model
      if (model_mode == 1) { # XGBoost
        inner_train_xgb <- xgb.DMatrix( data = X_inner_train, label = y_inner_train, 
                                        weight = weight_cv )
        inner_test_xgb <- xgb.DMatrix( data = X_inner_test, label = y_inner_test )
        cv_fit <- lapply( all_param_grid_split, function(x) {
          xgb.train( params = as.list(x[-nrounds_ind]), data = inner_train_xgb, 
                     nrounds = x[[nrounds_ind]], verbose = 0 )
        } )
        test_preds <- lapply( cv_fit, predict, newdata = inner_test_xgb )
      }
      else { # Random Forest
        cv_fit <- lapply( all_param_grid_split, function(z) {
          randomForest( x = X_inner_train, y = factor(y_inner_train),
                        ntree = z$ntree, nodesize = z$nodesize, mtry = z$mtry ) 
        } )
        test_preds <- lapply( cv_fit, function(x) { predict( x, newdata = X_inner_test, type = 'prob' )[, 2] } )
      }
      all_test_preds[j, ] <- sapply( test_preds, function(x){ log_loss(x, y_inner_test) } )
    }
    
    # Compile results of inner loop
    ce_mean <- apply( all_test_preds, 2, mean )
    best_ind <- which.min( ce_mean )
    best_params_all[[i]] <- best_params <- all_param_grid_split[[best_ind]]
    all_log_loss[i] <- ce_mean[best_ind]
    
    # Remake outer set with imputed data
    X_imp_outer <- dynamic_pmm( X, outer_fold_inds[[i]] )
    X_outer_train <- X_imp_outer[-outer_fold_inds[[i]], ]
    X_outer_test <- X_imp_outer[outer_fold_inds[[i]], , drop = FALSE]
    
    # Observation weighting based on imbalance
    weight_cv <- rep( 1, length(y_outer_train) )
    weight_cv[y_outer_train == 1] <- sum( 1 - y_outer_train ) / sum( y_outer_train )
    
    # Fit the desired model
    if (model_mode == 1) { # XGBoost
      outer_train_xgb <- xgb.DMatrix( data = X_outer_train, label = y_outer_train, 
                                      weight = weight_cv )
      outer_test_xgb <- xgb.DMatrix( data = X_outer_test, label = y_outer_test )
      outer_fit <- xgb.train( params = as.list(best_params[-nrounds_ind]), 
                              data = outer_train_xgb, nrounds = best_params[[nrounds_ind]] )
      outer_preds_train <- predict( outer_fit, newdata = outer_train_xgb )
      outer_preds_test <- predict( outer_fit, newdata = outer_test_xgb )
      fit_unify <- xgboost.unify( outer_fit, X_outer_train )
    }
    else { # Random Forest
      outer_fit <- randomForest( x = X_outer_train, y = factor(y_outer_train),
                                 ntree = best_params$ntree, nodesize = best_params$nodesize, 
                                 mtry = best_params$mtry )
      outer_preds_train <- predict( outer_fit, newdata = X_outer_train, type = 'prob' )[, 2]
      outer_preds_test <- predict( outer_fit, newdata = X_outer_test, type = 'prob' )[, 2]
      fit_unify <- randomForest2.unify( outer_fit, X_outer_train )
    }
    
    # Accuracy and Shapley values
    outer_preds_train_all[[i]] <- data.frame( preds = outer_preds_train, true = y_outer_train )
    outer_preds_test_all[[i]] <- data.frame( preds = outer_preds_test, true = y_outer_test )
    outer_shaps_all[[i]] <- treeshap( fit_unify, X_all_imp, verbose = FALSE )$shaps
    
  }
  
  ### Visualise results
  outer_preds_train_dat <- do.call( 'rbind', outer_preds_train_all )
  outer_preds_test_dat <- do.call( 'rbind', outer_preds_test_all )
  pred_label <- c( "None", "Exposure" )
  outer_preds_train_dat$label <- factor( pred_label[outer_preds_train_dat$true + 1], pred_label )
  outer_preds_test_dat$label <- factor( pred_label[outer_preds_test_dat$true + 1], pred_label )
  outer_preds_train_dat$true <- NULL
  outer_preds_test_dat$true <- NULL
  
  # Prediction vs Response
  pdf( "Predicted_Probability_Train.pdf", height = 10, width = 10 )
  plot( preds ~ label, data = outer_preds_train_dat, col = mycol, pch = 20, ylim = c(0, 1),
        xlab = "", ylab = "", main = "Predicted Probability of Exposure (Training Set)" )
  dev.off()
  
  pdf( "Predicted_Probability_Test.pdf", height = 10, width = 10 )
  plot( preds ~ label, data = outer_preds_test_dat, col = mycol, pch = 20, ylim = c(0, 1),
        xlab = "", ylab = "", main = "Predicted Probability of Exposure (Test Set)" )
  dev.off()
  
  # ROC curve
  pdf( "ROC.pdf", height = 10, width = 10 )
  roc_dat_train <- roc( label ~ preds, data = outer_preds_train_dat, direction = "<", quiet = TRUE )
  roc_dat_test <- roc( label ~ preds, data = outer_preds_test_dat, direction = "<", quiet = TRUE )
  plot( 1 - roc_dat_train$specificities, roc_dat_train$sensitivities, type = 'l', 
        ylab = "Sensitivity", xlab = "False Positive Rate", main = "ROC Curve",
        col = mycol[1], lwd = 2 )
  lines( 1 - roc_dat_test$specificities, roc_dat_test$sensitivities, col = mycol[2], lwd = 2 )
  abline( a = 0, b = 1, lty = 2, col = 'grey' )
  text( 0.8, 0.22, paste0("Train AUC = ", round(roc_dat_train$auc, 2)), cex = 2, col = mycol[1] )
  text( 0.8, 0.08, paste0("Test AUC = ", round(roc_dat_test$auc, 2)), cex = 2, col = mycol[2] )
  dev.off()
  
  # Shapley geometric median
  shap_dat <- shap_convert( outer_shaps_all )
  shap_gmedian <- Gmedian( t(shap_dat[, 3:(outer_folds+2)]) ) %>% t
  med_shap_list <- split( t(shap_gmedian), rep( 1:ncol(X), each = nrow(X) ) )
  names( med_shap_list ) <- colnames( X )
  shap_contrib <- do.call( 'cbind', med_shap_list ) %>% as.data.table
  shap_long <- shap.prep( shap_contrib = shap_contrib, X_train = X_all_imp )
  shap_long_split <- split( shap_long, shap_long$variable )
  shap_long$stdfvalue <- unlist( lapply( shap_long_split, function(x) { (rank(x$rfvalue)-1) / (nrow(x)-1) } ) )
  s_plot <- shap.plot.summary( shap_long )
  
  pdf( "Shapley_All.pdf", height = 10, width = 10 )
  print( s_plot )
  dev.off()
  
  # More condensed version
  cutoff <- min( which( shap_long$variable == 'random_num' ) ) - 1
  shap_long_trunc <- shap_long[1:cutoff, ]
  shap_long_trunc$variable <- factor( shap_long_trunc$variable )
  s_plot_trunc <- shap.plot.summary( shap_long_trunc )
  
  pdf( "Shapley.pdf", height = 10, width = 10 )
  print( s_plot_trunc )
  dev.off()
  
  ### Save chosen parameter data
  best_params_dat <- do.call( 'rbind', best_params_all )
  best_params_dat$param_grid_ind <- rownames( best_params_dat )
  best_params_dat$cv_error <- all_log_loss
  write_csv( best_params_dat, 'Chosen_Parameters.csv' )
  
  ### Save results
  res <- list( accuracy_train = outer_preds_train_dat, accuracy_test = outer_preds_test_dat, 
               shaps = outer_shaps_all, chosen_parameters = best_params_dat )
  save( res, file = 'res.Rdata' )
  setwd( '../..' )
  
}

