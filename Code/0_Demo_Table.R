library( readxl )
library( dplyr )

fn <- "Data/IRSABI - Final Set for Pateint Demographics - 180 patients.xlsx"

# Time on Dialysis
tod <- read_xlsx( fn, sheet = "Time on Dialysis" ) %>% select( "Months since starting Haemodialysis", "Exposure" ) %>%
  na.omit
tod$range <- ifelse( tod$`Months since starting Haemodialysis` > 48, 
                     ifelse(tod$`Months since starting Haemodialysis` > 96, '>96', '49-96'), '1-48' )
table( tod$Exposure, tod$range )

# Age at Consent
aoc <- read_xlsx( fn, sheet = "Age at consent" ) %>% select( "Age at Consent", "Exposure" ) %>%
  na.omit
aoc$range <- ifelse( aoc$`Age at Consent` > 45, 
                     ifelse(aoc$`Age at Consent` > 70, '>70', '46-70'), '20-45' )
table( aoc$Exposure, aoc$range )
table( aoc$Exposure, aoc$range ) %>% rowSums

# Sex
sex <- read_xlsx( fn, sheet = "Sex" ) %>% select( "Gender", "Exposure" ) %>%
  na.omit
table( sex$Exposure, sex$Gender )

# Aetiology
aockd <- read_xlsx( fn, sheet = "Aetiology of CKD" ) %>% select( "Aetiology of CKD", "Exposure" ) %>%
  na.omit
aockd$cd1 <- substr( aockd$`Aetiology of CKD`, 1, 1 ) %>% as.numeric
aockd$cd2 <- substr( aockd$`Aetiology of CKD`, nchar(aockd$`Aetiology of CKD`), nchar(aockd$`Aetiology of CKD`) ) %>% as.numeric

double_ind <- setdiff( 1:nrow(aockd), grep('&', aockd$`Aetiology of CKD`) )
aockd$cd2[double_ind] <- 8

causes <- c( "Chronic Glomerulonephritis", "Ischaemic Nephrology", "Polycystic Kidney Disease",
             "Diabetes", "Congenital", "Other", "Unknown", "None" )

aockd$cd1 <- factor( causes[aockd$cd1], levels = causes )
aockd$cd2 <- factor( causes[aockd$cd2], levels = causes )

table( aockd$Exposure, aockd$cd1 ) + table( aockd$Exposure, aockd$cd2 )

