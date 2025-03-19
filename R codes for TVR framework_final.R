# Set working directory
setwd("F:/Co-Work/Peipeixi Indulgence/R codes")

# Load necessary libraries
library(readr)   # For reading CSV files
library(lavaan)  # For structural equation modeling (SEM)
library(psych)   # For psychometric analyses

###################################################
### Data Processing: US Survey ###################

# Import US dataset with specified date formats
US <- read_csv("Dissertation survey- U.S_December 23, 2022_09.21.csv",
               col_types = cols(
                 StartDate = col_date(format = "%m/%d/%Y %H:%M"),
                 EndDate = col_date(format = "%m/%d/%Y %H:%M"),
                 RecordedDate = col_date(format = "%m/%d/%Y %H:%M")
               )
)

# Remove unfinished responses (where Finished column is 0)
US <- US[US$Finished != 0, ]

# Demographic Descriptions: US Survey

# Gender distribution
table(US$Q2)
round(table(US$Q2) / nrow(US) * 100, 2)

# Age group distribution
table(US$Q3)
round(table(US$Q3) / nrow(US) * 100, 2)

# Race distribution
round(table(US$Q16) / nrow(US) * 100, 2)

# Education level distribution
table(US$Q17)
round(table(US$Q17) / nrow(US) * 100, 2)

# Income distribution
round(table(US$Q18) / nrow(US) * 100, 2)

# Employment status distribution
table(US$Q19)
round(table(US$Q19) / nrow(US) * 100, 2)

# Political party affiliation
round(table(US$Q20) / nrow(US) * 100, 2)

# Political orientation
round(table(US$Q21) / nrow(US) * 100, 2)

# Remove unnecessary columns (columns 1-6 and 8-17)
US <- US[, -c(1:6, 8:17)]

# Select columns related to Indulgence and IBB
US <- US[, c(33:38, 39:45, 47)]  # Indulgence (33-38) and IBB (39-45, 47)

###################################################
### Data Processing: India Survey ################

# Import India dataset with specified date formats
Ind <- read_csv("Dissertation survey- India_December 23, 2022_09.19.csv",
                col_types = cols(
                  StartDate = col_date(format = "%m/%d/%Y %H:%M"),
                  EndDate = col_date(format = "%m/%d/%Y %H:%M"),
                  RecordedDate = col_date(format = "%m/%d/%Y %H:%M")
                )
)

# Remove unfinished responses
Ind <- Ind[-which(Ind$Finished == 0), ]

# Remove responses not meeting specific conditions
Ind <- Ind[Ind$`Q8?Masculinity_3` == 3, ]

# Demographic Descriptions: India Survey

# Gender distribution
table(Ind$Q2)
round(table(Ind$Q2) / nrow(Ind) * 100, 2)

# Age group distribution
table(Ind$Q3)
round(table(Ind$Q3) / nrow(Ind) * 100, 2)

# Education level distribution
table(Ind$Q17)
round(table(Ind$Q17) / nrow(Ind) * 100, 2)

# Income distribution
round(table(Ind$Q18) / nrow(Ind) * 100, 2)

# Employment status distribution
table(Ind$Q19)
round(table(Ind$Q19) / nrow(Ind) * 100, 2)

# Remove unnecessary columns
Ind <- Ind[, -c(1:6, 8:17)]

# Select specific columns for analysis
Ind <- Ind[, c(33:38, 39:45, 47)]  # Indulgence (33-38) and IBB (39-45, 47)

###################################################
### Data Processing: China Survey #################

# Load necessary library for reading Excel files
library(readxl)

# Import China dataset
CN <- read_excel("204198085_2_NewChina_259_259.xlsx")

# Remove responses that do not meet the required conditions
CN <- CN[CN$Q1 != 2, ]

# Demographic Descriptions: China Survey

# Gender distribution
table(CN$Q2)
round(table(CN$Q2) / nrow(CN) * 100, 2)

# Age group distribution
table(CN$Q3)
round(table(CN$Q3) / nrow(CN) * 100, 2)

# Education level distribution
table(CN$Q16)
round(table(CN$Q16) / nrow(CN) * 100, 2)

# Income distribution
table(CN$Q17)
round(table(CN$Q17) / nrow(CN) * 100, 2)

# Employment status distribution
table(CN$Q18)
round(table(CN$Q18) / nrow(CN) * 100, 2)

# Race distribution
table(CN$Q19)
round(table(CN$Q19) / nrow(CN) * 100, 2)

# Area distribution
table(CN$Q20)
round(table(CN$Q20) / nrow(CN) * 100, 2)

# Remove unnecessary columns
CN <- CN[, -c(1:6)]

# Select specific columns for analysis
CN <- CN[, c(31:36, 37:43, 45)]  # Indulgence (31-36) and IBB (37-43, 45)

###################################################
### Renaming Variables ###########################

# Create standardized column names for Indulgence (IDG) and IBB
namesIDG_IBB <- c(
  paste(rep("IL_IDG", 6), 1:6, sep = ""),
  paste(rep("IBB", 8), 1:8, sep = "")
)

namesIDG <- namesIDG_IBB[1:6]

# Apply the new column names across datasets
names(CN) <- namesIDG_IBB
names(US) <- namesIDG_IBB
names(Ind) <- namesIDG_IBB

###################################################
### Reliability Analysis #########################

# Ensure `psych` package is installed for Cronbach's alpha
if (!requireNamespace("psych", quietly = TRUE)) {
  install.packages("psych")
}
library(psych)

# Calculate Cronbach's alpha for each dataset (Indulgence subscale)

# US dataset
alpha_results_US <- alpha(US[, c("IL_IDG1", "IL_IDG2", "IL_IDG3", "IL_IDG4", "IL_IDG5", "IL_IDG6")])
print(alpha_results_US)

# India dataset
alpha_results_Ind <- alpha(Ind[, c("IL_IDG1", "IL_IDG2", "IL_IDG3", "IL_IDG4", "IL_IDG5", "IL_IDG6")])
print(alpha_results_Ind)

# China dataset
alpha_results_CN <- alpha(CN[, c("IL_IDG1", "IL_IDG2", "IL_IDG3", "IL_IDG4", "IL_IDG5", "IL_IDG6")])
print(alpha_results_CN)

###################################################
### Confirmatory Factor Analysis (CFA) ###########

# Ensure the `lavaan` package is installed
if (!requireNamespace("lavaan", quietly = TRUE)) {
  install.packages("lavaan")  # Install lavaan for Structural Equation Modeling (SEM)
}
library(lavaan)  # Load the `lavaan` package for CFA and SEM analyses

# Define the baseline measurement model
# IDG is a latent variable measured by six observed indicators (IL_IDG1 - IL_IDG6)
original_structure <- '
  # Measurement model
  IDG =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4 + IL_IDG5 + IL_IDG6
'

############## US Data CFA #######################

# Fit the CFA model for the US dataset
TestingIDG_US <- cfa(
  original_structure, 
  data = US, 
  std.lv = TRUE, 
  estimator = "MLR"
)

# Display model fit and parameter estimates
summary(
  TestingIDG_US,
  fit.measures = TRUE,  # Include model fit indices
  standardized = TRUE   # Show standardized estimates for interpretation
)

############## India Data CFA ###################

# Fit the CFA model for the India dataset
TestingIDG_Ind <- cfa(
  original_structure, 
  data = Ind, 
  std.lv = TRUE, 
  estimator = "MLR"
)

summary(
  TestingIDG_Ind,
  fit.measures = TRUE, 
  standardized = TRUE
)

############## China Data CFA ####################

# Fit the CFA model for the China dataset
TestingIDG_CN <- cfa(
  original_structure, 
  data = CN, 
  std.lv = TRUE, 
  estimator = "MLR"
)

summary(
  TestingIDG_CN,
  fit.measures = TRUE,  
  standardized = TRUE
)

##########################################################
### Refining CFA Model Across Multiple Groups ###########

# Create copies of datasets and assign country identifiers
CN_multi <- CN; CN_multi$CI <- 1  # China = 1
Ind_multi <- Ind; Ind_multi$CI <- 2  # India = 2
US_multi <- US; US_multi$CI <- 3  # US = 3

# Merge all datasets into one for multi-group CFA
Full_data <- rbind(CN_multi, Ind_multi, US_multi)

# Prepare empty lists for constraints
pre_load <- list()   # Items to be estimated freely
pre_unload <- list() # Items to be constrained to zero

# Define fit indices to be used in model optimization
# For instance: 1 = SRMR, 4 = TLI, 5 = CFI
Fit_Use <- c(5, 4, 1) 

# Run the heuristic Ant Colony Optimization (hACO) algorithm to refine model structure (2-factor solution)
results <- hACO(
  Data = as.data.frame(Full_data),              
  n_factors = 2,  # Specify number of latent factors           
  n_variables = 6, # Number of observed variables       
  loaded = pre_load, 
  unloaded = pre_unload, 
  AccList = NULL, 
  BanList = NULL, 
  Phe_level = NULL, 
  standards = c(0.93, 1-0.07, 0.93),  # Fit thresholds (CFI >= 0.93, SRMR <= 0.07, etc.)
  maxRun = 1000,  
  Punish_rate = 0.25, 
  Fit = Fit_Use                     
)

# Fit the refined model using the optimized solution
RefiningIDG_configural <- cfa(
  model = results$OptModelRun,  # Optimal model output from hACO
  data = Full_data,             
  std.lv = TRUE,                
  estimator = "MLR",            
  group = "CI"                  
)

# Display refined model summary
summary(
  RefiningIDG_configural, 
  fit.measures = TRUE, 
  standardized = TRUE
)

###################################################
### Validating the Refined Model #################

# Define refined measurement model (2-factor solution)
refined_structure <- '
  # Measurement model with two correlated latent factors
  IDGF1 =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4  
  IDGF2 =~ IL_IDG5 + IL_IDG6  

  # Correlation between the two IDG factors
  IDGF1 ~~ IDGF2 
'

# Fit the configural invariance model (allows parameters to vary across groups)
ValidatingIDG_configural <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI"
)

summary(
  ValidatingIDG_configural,
  fit.measures = TRUE,  
  standardized = TRUE
)

# Fit metric invariance model (constrains factor loadings across groups)
ValidatingIDG_metric <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI",  
  group.equal = "loadings"
)

summary(
  ValidatingIDG_metric,
  fit.measures = TRUE,  
  standardized = TRUE
)

# Fit scalar invariance model (constrains loadings and intercepts across groups)
ValidatingIDG_scalar <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI",  
  group.equal = c("loadings", "intercepts")
)

summary(
  ValidatingIDG_scalar,
  fit.measures = TRUE,  
  standardized = TRUE
)

# Fit strict invariance model (also constrains residuals across groups)
ValidatingIDG_strict <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI",  
  group.equal = c("loadings", "intercepts", "residuals")
)

summary(
  ValidatingIDG_strict,
  fit.measures = TRUE,  
  standardized = TRUE
)

######################################################
### Empirical Study: Testing Structural Model #######

# Define structural model with original one-factor IDG
structural_model_Original <- '
  IDG =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4 + IL_IDG5 + IL_IDG6  
  IBB =~ IBB1 + IBB2 + IBB3 + IBB4 + IBB5 + IBB6 + IBB7 + IBB8  
  IBB ~ IDG  
'

# Fit the original structural model
fit_original <- sem(
  model = structural_model_Original,  
  data = Full_data,  
  group = "CI",  
  estimator = "MLR",  
  std.lv = TRUE,  
  group.equal = c("loadings")  
)

summary(
  fit_original,
  fit.measures = TRUE,  
  standardized = TRUE  
)

# Define structural model with refined two-factor IDG
structural_model_Refined <- '
  IDGF1 =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4  
  IDGF2 =~ IL_IDG5 + IL_IDG6  
  IBB =~ IBB1 + IBB2 + IBB3 + IBB4 + IBB5 + IBB6 + IBB7 + IBB8  

  # Correlation between IDG factors
  IDGF1 ~~ IDGF2  

  # Structural relationships
  IBB ~ IDGF1 + IDGF2  
'

# Fit the refined structural model
fit_refined <- sem(
  model = structural_model_Refined,  
  data = Full_data,  
  group = "CI",  
  estimator = "MLR",  
  std.lv = TRUE,  
  group.equal = c("loadings")  
)

summary(
  fit_refined,
  fit.measures = TRUE,  
  standardized = TRUE  
)

