##### Logistic Regression LASSO fit for establishing baseline
##### Daniel Dempsey

### Read in libraries
library( readr ) # Read in csv files
library( dplyr ) # Dataframe manipulation tools
library( caret ) # For stratified sampling of training/test splits
library( glmnet ) # Fit a LASSO model
library( pROC ) # For computing ROC curve
mycol <- c( '#56B4E9', '#E69F00' )

### Function to make directories if they don't already exist
mkdir <- function( x ) {
  dir_split <- strsplit( x, '/' )[[1]]
  current_dir <- ''
  for ( i in 1:length(dir_split) ) {
    current_dir <- paste0(current_dir, dir_split[i], sep = '/')
    if( !dir.exists( current_dir ) ) {
      dir.create( current_dir )
    }
  }
}

### Load and prepare data
hd_bsi <- read_csv( 'Data/hd_bsi.csv' )
X <- select( hd_bsi, -Exposure ) %>% as.matrix
y <- hd_bsi$Exposure

### Use CV to find penalty parameter
home_dir <- 'Output/Model_Fits/LASSO'
mkdir( home_dir )
setwd( home_dir )

set.seed( 1001 )
cvfit <- cv.glmnet( X, y, nfolds = nrow(X), type.measure = "deviance", family = "binomial" )
lambda_min <- cvfit$lambda.min

png( 'CV_Deviance.png', width = 700, height = 600 )
plot( cvfit )
dev.off()

### Bootstrap training/test split
B <- 1000
LASSO_res <- matrix( FALSE, ncol = ncol(X), nrow = B )
colnames( LASSO_res ) <- colnames( X )
AUC_res <- matrix( 0, ncol = 2, nrow = B )
colnames( AUC_res ) <- c( 'Train', 'Test' )

set.seed( 10101 )
for ( i in 1:B ) {
  
  train_inds <- createDataPartition( as.factor(y), p = 0.8, list = FALSE )
  test_inds <- setdiff( 1:nrow(X), train_inds )
  
  X_train <- X[train_inds, ]
  X_test <- X[test_inds, ]
  
  y_train <- y[train_inds]
  y_test <- y[test_inds]
  
  fit <- glmnet( X_train, y_train, lambda = lambda_min, family = "binomial" )
  nonshrunk <- rownames( fit$beta )[ which( fit$beta != 0 ) ]
  LASSO_res[i, nonshrunk] <- TRUE
  
  preds_train <- predict( fit, newx = X_train, type = 'response' ) %>% as.vector
  preds_test <- predict( fit, newx = X_test, type = 'response' ) %>% as.vector
  
  AUC_res[i, ] <- c( roc( y_train, preds_train, direction = '<', quiet = TRUE )$auc, 
                     roc( y_test, preds_test, direction = '<', quiet = TRUE )$auc )
  
}

### Plot results
# LASSO-chosen predictors
mean_LASSO <- apply( LASSO_res, 2, mean ) %>% sort %>% rev
mean_LASSO <- mean_LASSO[mean_LASSO > 0]
use_col <- rep( mycol[1], length(mean_LASSO) )
use_col[which(names(mean_LASSO) == 'random_num')] <- mycol[2]

png( 'LASSO_Vars.png', width = 700, height = 600 )
par( mar = c(12, 4, 4, 2) + 0.1 )
barplot( mean_LASSO, col = use_col, las = 2, main = 'LASSO Variables' )
par( mar = c(5, 4, 4, 2) + 0.1 )
dev.off()

# AUC
png( 'AUC.png', width = 700, height = 600 )
boxplot( AUC_res, col = mycol, pch = 20, yaxt = 'n', ylim = c(0, 1), main = 'LASSO AUC' )
y_axt <- seq( 0, 1, 0.25 )
axis( 2, y_axt, y_axt, las = 1 )
abline( h = seq(0, 1, 0.25), lty = 3, col = 'grey' )
abline( h = 0.5, lty = 2 )
dev.off()

