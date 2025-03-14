##### Data cleaning and visualisation of hemodialysis (HD) and bloodstream infection (BSI) data
##### Daniel Dempsey

### Setting up
# Load in libraries
library( readxl ) # Reading in excel files
library( readr ) # Reading and writing csv files
library( dplyr ) # Data frame manipulation tools
library( caret ) # Used for easy one-hot encoding
library( stringi ) # Used when cleaning column names
library( corrplot ) # For visualising correlations
library( mice ) # Missing value imputation
library( ggforce ) # Sina plots
library( ggpubr ) # For placing multiple graphs on the same plot
library( naniar ) # Missing data plots
library( ggcorrplot ) # Correlation matrix plots
mycol <- c( '#56B4E9', '#E69F00' )

### Utility function: make directories if they don't already exist
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
mkdir( 'Output/Data_Visualisations' )

### Read in data, perform initial cleaning and make column names consistent
### Clinic data
hd_raw <- read_xlsx( 'Data/HD-180 - FACS and excel matching -Feb 2023.xlsx', skip = 4, sheet = 'US & 3 antigens pooled' )

# Two columns list study numbers; check that they are consistent
all( hd_raw$StudyNo...1 == hd_raw$StudyNo...34 ) # All TRUE

# Delete one of the redundant StudyNo columns
studynum_inds <- which( substr(colnames(hd_raw), 1, 7) == 'StudyNo' ) # Two columns 
hd <- select( hd_raw, -studynum_inds[2] )
colnames( hd )[1] <- 'StudyNumber'

# Check that the two Exposure variables are consistent
table( hd$HospitaldocumentedSAExposure, hd$Exposure ) # Consistent, remove the redundant column
hd$HospitaldocumentedSAExposure <- NULL

# Change column names to be consistent with BSI set
colnames( hd )[3] <- 'Age'
hd$Sex <- hd_raw$Gender - 1
hd$Gender <- NULL
colnames( hd )[17] <- 'NonClassicalMonocytes'

# Create response variable
hd$status <- ifelse( hd_raw$Exposure == 1, 'Prior Exposure', 'None' )

### BSI data
# Easier if data and column names are read seperately due to formatting in original file
bsi_names_raw <- read_xlsx( 'Data/BSI Pooled Data - Mem Prolif Cyto - 3 Strains.xlsx',
                            range = 'A3:CU3', col_names = FALSE ) %>% unlist

bsi_dat_raw <- read_xlsx( 'Data/BSI Pooled Data - Mem Prolif Cyto - 3 Strains.xlsx',
                          range = 'A4:CU28', col_names = FALSE )

# Check if duplicated names have the same data
dup_fun <- function( x ) {
  dup_flags <- duplicated( x )
  which( x %in% x[ dup_flags ] )
}

bsi_dup_inds <- dup_fun( bsi_names_raw )
bsi_names_raw[ bsi_dup_inds ]

# Check the two CD3 features
all( bsi_dat_raw[, bsi_dup_inds[1]] == bsi_dat_raw[, bsi_dup_inds[2]] ) # Different
bsi_names_raw[ bsi_dup_inds[1] ] <- 'CD3_2'

# Clean the names and attach to dataset
clean_cols <- function( x ) {
  stri_replace_all_regex( x, c('\\?', ' ', '\\)'), '', vectorize = FALSE ) %>%
    stri_replace_all_regex( '@', '_at_', vectorize = FALSE ) %>%
    stri_replace_all_regex( '/', '_or_', vectorize = FALSE ) %>%
    stri_replace_all_regex( '\\+ve', '_positive', vectorize = FALSE ) %>%
    stri_replace_all_regex( c('-', '\\('), '', vectorize = FALSE )
}

bsi_names <- clean_cols( bsi_names_raw )
colnames( bsi_dat_raw ) <- bsi_names
bsi <- bsi_dat_raw

# Fix mistake in the dataset with one of the Transplant Ever values
bsi$TransplantEver[bsi$StudyNo == "CKD032-BSI(1)"] <- 0

# Fix study ID
colnames( bsi )[2] <- 'StudyNumber'
bsi$Infection <- regmatches( bsi$StudyNumber, gregexpr("(?<=\\().*?(?=\\))", bsi$StudyNumber, perl = T) ) %>% 
  unlist %>% as.numeric
bsi$StudyNumber <- sub( '-.*', '', bsi$StudyNumber )

# Fix sex
bsi$Sex <- bsi_dat_raw$Gender - 1
bsi$Gender <- NULL

# Create response variable
bsi$status <- 'Current Exposure'

# Alter cell names for consistency
name_swap <- function( old_name, new_name ) {
  nm <- colnames( bsi )
  inds <- grep( old_name, nm )
  colnames( bsi )[inds] <<- gsub( old_name, new_name, nm[inds] )
}
name_swap( 'MemoryCD4', 'MemCD4' )
name_swap( 'Dividing', 'Div' )
name_swap( 'corIL10', 'IL10cor' )
name_swap( 'TransitioningTh17', 'TransTh17' )

# Manually change Th17 and Th1 to avoid changing things I don't want to
colnames( bsi )[84:86] <- c( "Th17corHKPS80", "Th17corHKMRSA", "Th17corHKMSSA" )
colnames( bsi )[88:90] <- c( "TransTh17corHKPS80", "TransTh17corHKMRSA", "TransTh17corHKMSSA" )
colnames( bsi )[92:94] <- c( "Th1corHKPS80", "Th1corHKMRSA", "Th1corHKMSSA" )
colnames( bsi )[96:98] <- c( "exTh17corHKPS80", "exTh17corHKMRSA", "exTh17corHKMSSA")

### Merge data
keep_cols <- intersect( names(hd), names(bsi) )
hd_bsi_raw <- rbind( hd[, keep_cols], bsi[, keep_cols] )
hd_bsi <- hd_bsi_raw

# Drop unstimulated and non-corrected variables
drop_names <- c( "CD4Unstim", "CD4HKPS80", "CD4HKMRSA", "CD4HKMSSA", 
                 "MemCD4Unstim", "MemCD4HKPS80", "MemCD4HKMRSA", "MemCD4HKMSSA",
                 "IL10Unstim", "Th17Unstim", "TransTh17Unstim", "Th1Unstim", "exTh17Unstim" )
hd_bsi <- select( hd_bsi_raw, -all_of(drop_names) )

### Fix categorical data
# CurrentImmunosuppressiveRx
hd_bsi$CurrentImmunosuppressiveRx <- ifelse( toupper(hd_bsi_raw$CurrentImmunosuppressiveRx) == "Y", 1, 0 )

# Aetiology of CKD
hd_bsi$AetiologyofCKD[which(hd_bsi_raw$AetiologyofCKD == '1&2')] <- '12'
hd_bsi$AetiologyofCKD[which(hd_bsi_raw$AetiologyofCKD == '2&4')] <- '24'
hd_bsi$AetiologyofCKD[which(hd_bsi_raw$AetiologyofCKD == '2&5')] <- '25'
hd_bsi$AetiologyofCKD[which(hd_bsi_raw$AetiologyofCKD == '2 &7')] <- '27'
hd_bsi$AetiologyofCKD <- as.numeric( hd_bsi$AetiologyofCKD )
hd_bsi$AetiologyofCKD <- as.factor( hd_bsi$AetiologyofCKD )

# Vascular Access
hd_bsi$VascularAccess[which(hd_bsi_raw$VascularAccess == '1&2')] <- '12'
hd_bsi$VascularAccess <- as.numeric( hd_bsi$VascularAccess )
hd_bsi$VascularAccess <- as.factor( hd_bsi$VascularAccess )

# Diabetes
hd_bsi$Diabetes <- ifelse( tolower(hd_bsi_raw$Diabetes) %in% c('1', 'yes'), 1, 0 )

### One-hot encoding of non-binary categorical variables
dv <- dummyVars( ~ VascularAccess + AetiologyofCKD, hd_bsi )
new_facs <- predict( dv, hd_bsi ) %>% as.data.frame

# Resolve categorical variables where multiple classes are listed
va12_inds <- which( new_facs$VascularAccess.12 == 1 )
new_facs$VascularAccess.1[va12_inds] <- new_facs$VascularAccess.2[va12_inds] <- 1
new_facs$VascularAccess.12 <- NULL

ackd12_inds <- which( new_facs$AetiologyofCKD.12 == 1 )
new_facs$AetiologyofCKD.1[ackd12_inds] <- new_facs$AetiologyofCKD.2[ackd12_inds] <- 1
new_facs$AetiologyofCKD.12 <- NULL

ackd24_inds <- which( new_facs$AetiologyofCKD.24 == 1 )
new_facs$AetiologyofCKD.2[ackd24_inds] <- new_facs$AetiologyofCKD.4[ackd24_inds] <- 1
new_facs$AetiologyofCKD.24 <- NULL

ackd25_inds <- which( new_facs$AetiologyofCKD.25 == 1 )
new_facs$AetiologyofCKD.2[ackd25_inds] <- new_facs$AetiologyofCKD.5[ackd25_inds] <- 1
new_facs$AetiologyofCKD.25 <- NULL

ackd27_inds <- which( new_facs$AetiologyofCKD.27 == 1 )
new_facs$AetiologyofCKD.2[ackd27_inds] <- new_facs$AetiologyofCKD.7[ackd27_inds] <- 1
new_facs$AetiologyofCKD.27 <- NULL

### Assign more informative names to the categorical variables
read_clean <- function( x ) {
  read_xlsx( 'Data/Clinical Data Codes.xlsx', range = x ) %>% unlist %>% 
    gsub( pattern = '[^[:alnum:]]', replacement = '_' )
}

all_codes <- c( read_clean('A10:A13'), read_clean('A18:A25') )
colnames( new_facs ) <- mapply( function(x, y) { sub('\\..*', paste0('.', y), x) }, 
                                x = colnames(new_facs), y = all_codes )

### Combine one-hot-encoded variables and avoid dummy variable trap
hd_bsi <- cbind( hd_bsi, new_facs )
hd_bsi$VascularAccess.No_access <- hd_bsi$AetiologyofCKD.Unknown <- NULL

### Missing value visualisation and exceptions
setwd( 'Output/Data_Visualisations' )
png( 'Missing_Values_rawData.png', width = 700, height = 600 )
vis_miss( hd_bsi )
dev.off()

# Drop rows with >= 10 missing values
mis_tab <- table( which(is.na(hd_bsi), arr.ind = TRUE)[, 1] )
drop_rows <- names( mis_tab[mis_tab >= 10] ) %>% as.numeric
hd_bsi_trim <- filter( hd_bsi, !(1:nrow(hd_bsi) %in% drop_rows) )

png( 'Missing_Values.png', width = 700, height = 600 )
vis_miss( hd_bsi_trim )
dev.off()

### Correlation Matrix
drop_cor <- c( "StudyNumber", "Diabetes", "TransplantEver", "CurrentImmunosuppressiveRx",
               "Sex", "status", "VascularAccess.T_CVC", "VascularAccess.AVF",                    
               "AetiologyofCKD.Chronic_Glomerulonephritis", "AetiologyofCKD.Ischaemic_Nephrology",      
               "AetiologyofCKD.Polycystic_Kidney_Disease", "AetiologyofCKD.Diabetes",                  
               "AetiologyofCKD.Congenital", "AetiologyofCKD.Other", "VascularAccess",
               "AetiologyofCKD" )

png( 'Correlation_Matrix.png', width = 800, height = 800 )
select( hd_bsi_trim, -all_of(drop_cor) ) %>% cor( use = "complete.obs" ) %>% corrplot
dev.off()

# Drop variables due to multicollinearity (the others will be handled via PCA), and also others for relevancy reasons
hd_bsi_trim2 <- select( hd_bsi_trim, -all_of(c('CD8', 'EffectorTcells', 'CD3', 'Diabetes')) )

### Dimension reduction of different Day 8 strains using PCA
PCA_vars <- c( "MemCD4corHKPS80", "MemCD4corHKMRSA", "MemCD4corHKMSSA",                          
               "DivMemCD4corHKPS80", "DivMemCD4corHKMRSA", "DivMemCD4corHKMSSA",                        
               "IL10corHKPS80", "IL10corHKMRSA", "IL10corHKMSSA", "Th17corHKPS80", 
               "Th17corHKMRSA", "Th17corHKMSSA", "TransTh17corHKPS80", 
               "TransTh17corHKMRSA", "TransTh17corHKMSSA", "Th1corHKPS80", 
               "Th1corHKMRSA", "Th1corHKMSSA", "exTh17corHKPS80",                          
               "exTh17corHKMRSA", "exTh17corHKMSSA" )

PCA_var_split <- split( PCA_vars, c(rep(1, 6), rep(2:6, each = 3)) )
names( PCA_var_split ) <- c( 'memCD4', 'IL10', 'Th17', 'transTh17', 'Th1', 'exTh17' )

PCA_fun <- function( x ) {
  dat <- hd_bsi_trim2[, x] %>% na.omit
  prcomp( dat, scale. = TRUE )
} 

PCA_list_raw <- lapply( PCA_var_split, PCA_fun )
PCA_list <- lapply( PCA_list_raw, function(x) { 
  if ( all(x$rotation[, 1] < 0) ) {
    x$rotation[, 1] <- -x$rotation[, 1]
  }
  
  x } ) # This is so the first principal component is easier to interpret

# Highlight principal components that explain > 80% of the variance
plot_PCA <- function( x, title, ... ) {
  
  var_ex <- x$sdev^2 / sum( x$sdev^2 ) 
  var_ex_cumsum <- cumsum( var_ex )
  last_pca <- which( var_ex_cumsum > 0.8 )[1]
  col_ind <- c( rep( 1, last_pca ), rep( 2, ncol(x$x) - last_pca ) )
  
  var_ex_dat <- data.frame( Var1 = paste0('PC', 1:length(col_ind)), Freq = var_ex )
  var_ex_cumsum_dat <- data.frame( Var1 = paste0('PC', 1:length(col_ind)), Freq = var_ex_cumsum )
  
  ggplot( var_ex_dat, aes( x = Var1, y = Freq ) ) +
    geom_bar(stat="identity", fill=c('cyan', 'grey')[col_ind], col = 'black') + ylim(0, 1) +
    ylab('Proportion of Variance Explained') + xlab('') + ggtitle(paste0(title, ' Principal Components')) +
    theme( text = element_text(size = 15) )
  
}

PCA_plots <- Map( plot_PCA, x = PCA_list, title = names( PCA_list ) )

png( 'PCA_Variance_Explained.png', height = 700, width = 1000 )
ggarrange( plotlist = PCA_plots, nrow = 2, ncol = 3 )
dev.off()

add_cols <- function( var_name, i ) {
  
  x <- PCA_list[[var_name]]
  var_ex <- x$sdev^2 / sum( x$sdev^2 ) 
  var_ex_cumsum <- cumsum( var_ex )
  var_ind <- 1:which( var_ex_cumsum > 0.8 )[1]
  
  PCs <- (as.matrix(hd_bsi_trim2[, PCA_var_split[[i]]]) %*% x$rotation) %>% apply( 2, scale )
  var_pca <- PCs[, var_ind, drop = FALSE]
  colnames( var_pca ) <- paste0( paste0(var_name, '_PC'), var_ind )
  cbind( hd_bsi_pca, var_pca )
  
}

hd_bsi_pca <- hd_bsi_trim2

for ( i in 1:length(PCA_list) ) {
  hd_bsi_pca <- add_cols( names(PCA_list)[i], i )
}

final_dat <- select( hd_bsi_pca, -all_of(PCA_vars) )

### Add a noise variable
set.seed( 1234 )
final_dat$random_num <- abs( rnorm(nrow(final_dat)) )

###### Data Visualisation
final_dat$Exposure <- ifelse( final_dat$status == 'Current Exposure', 1, 0 )

### Plot response data
expose_dat <- as.data.frame( table( final_dat$Exposure ) )

exposure_plot <- ggplot( expose_dat, aes( x = Var1, y = Freq ) ) +
  geom_bar(stat="identity", fill=mycol)+
  ylab('') + xlab('Exposure Status') + ggtitle('Proportion of Exposure Status')

png( 'Exposure_Barplot.png', width = 700, height = 600 )
exposure_plot
dev.off()

### Plot categorical covariates
tab_fun <- function( x ) {
  tab <- as.data.frame( table(final_dat[[x]], final_dat$Exposure) )
  colnames(tab)[1:2] <- c( x, 'Exposure' )
  tab
}

cat_list <- list()

cat_list[[1]] <- tab_fun( 'Sex' )
g_nms <- c('Male', 'Female')
cat_list[[1]]$Sex <- rep( g_nms, 2 )
cat_list[[1]]$Sex <- factor( cat_list[[1]]$Sex, levels = g_nms )
cat_list[[2]] <- tab_fun( 'CurrentImmunosuppressiveRx' )
cat_list[[3]] <- tab_fun( 'TransplantEver' )
cat_list[[4]] <- tab_fun( 'VascularAccess' )
v_nms <- c('T-CVC', 'AVF', 'AVG', 'T-CVC & AVF')
cat_list[[4]]$VascularAccess <- rep( v_nms, 2 )
cat_list[[4]]$VascularAccess <- factor( cat_list[[4]]$VascularAccess, levels = v_nms )
cat_list[[5]] <- tab_fun( 'AetiologyofCKD' )
a_nms <- c('Chronic Glomerulonephritis', 'Ischaemic Nephrology', 
           'Polycystic Kidney Disease', 'Diabetes', 'Congenital',
           'Other', 'Unknown', 'Chronic Glomerulonephritis & Ischaemic Nephrology',
           'Ischaemic Nephrology & Diabetes', 'Ischaemic Nephrology & Congenital',
           'Ischaemic Nephrology & Unknown')
cat_list[[5]]$AetiologyofCKD <- rep( a_nms, 2 )
cat_list[[5]]$AetiologyofCKD <- factor( cat_list[[5]]$AetiologyofCKD, levels = a_nms )

mkdir( 'Categorical_Barplots' )
for( i in 1:length(cat_list) ) {
  
  nm <- colnames( cat_list[[i]] )[1]
  png( paste0('Categorical_Barplots/', nm, '_Barplot.png'), width = 700, height = 600 )
  cat_plot <- ggplot(data=cat_list[[i]], aes(x=cat_list[[i]][, 1], y=Freq, fill=Exposure)) +
    geom_bar(stat="identity", position=position_dodge()) + ylab('') + xlab( nm ) +
    scale_fill_manual(values=mycol) + ggtitle( paste0('Distribution of ', nm) ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print( cat_plot )
  dev.off()
  
}

### Plot numerical covariates
cat_vars <- c( 'StudyNumber', 'Sex', 'CurrentImmunosuppressiveRx', 'VascularAccess',
               'AetiologyofCKD', 'TransplantEver', 'VascularAccess.T_CVC', 'VascularAccess.AVF', 
               'AetiologyofCKD.Chronic_Glomerulonephritis', 'AetiologyofCKD.Ischaemic_Nephrology', 
               'AetiologyofCKD.Polycystic_Kidney_Disease', 'AetiologyofCKD.Diabetes', 
               'AetiologyofCKD.Congenital', 'AetiologyofCKD.Other', 'status' )

num_dat <- select( final_dat, -all_of(cat_vars) ) 
num_dat$Exposure <- factor( num_dat$Exposure )

group_cols <- mycol
use_cols <- group_cols[ match( num_dat$Exposure, c( 1, 0 ) ) ]

mkdir( 'Numerical_Violins' )
set.seed( 1010 )
for ( i in 1:ncol( num_dat )  ) {
  
  nm <- colnames( num_dat )[i]
  if (nm == 'Exposure') { next }
  png( paste0('Numerical_Violins/', nm, '_Violin.png'), width = 700, height = 600 )
  vplot <- ggplot(num_dat, aes(x = Exposure, y = num_dat[[i]])) +
    geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
    geom_sina( aes(x = Exposure, y = num_dat[[i]], color = use_cols), method = "counts", maxwidth = 0.7, 
               alpha = 0.7 ) +
    theme(legend.position = "none") +
    ggtitle( paste0('Distribution of ', nm, ' by Exposure Status') ) + ylab( nm ) + xlab( 'Exposure Status' )
  print( vplot )
  dev.off()
  
}

### Remove redundant columns
drop_final <- c( "StudyNumber", "VascularAccess", "AetiologyofCKD", "status" )
final_dat_trim <- select( final_dat, -all_of(drop_final) )

### Missing value imputation
set.seed( 987 )
final_dat_complete <- complete( mice(final_dat_trim) )

### Write final data to file
write_csv( final_dat_complete, '../../Data/hd_bsi.csv' )

