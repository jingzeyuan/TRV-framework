# ==================================================================================
# CROSS-CULTURAL VALIDATION OF THE INDULGENCE SCALE USING THE TVR FRAMEWORK
# Testing-Validating-Refining (TVR) Framework for Addressing Construct Bias
# ==================================================================================

# Configuration 
setwd("H:/Co-Work/Peipeixi Indulgence/R_Indulgence")

# Import required libraries
library(readr)   # Data import utilities
library(lavaan)  # Structural equation modeling framework
library(psych)   # Psychometric analysis tools
source("hACO_MultiGroup_CFA.R")  # Hybrid Ant Colony Optimization algorithm for multi-group CFA

# ==================================================================================
# PHASE 1: DATA ACQUISITION AND PREPROCESSING
# ==================================================================================

cat("\n=========================================================\n")
cat("PHASE 1: DATA ACQUISITION AND PREPROCESSING\n")
cat("=========================================================\n")

# -------------------------------------------------------------------------------
# US DATASET PROCESSING
# -------------------------------------------------------------------------------

cat("\n--- US Dataset Processing ---\n")

# Import US dataset with appropriate date format specifications
US <- read_csv("US_Survey.csv",
               col_types = cols(
                 StartDate = col_date(format = "%m/%d/%Y %H:%M"),
                 EndDate = col_date(format = "%m/%d/%Y %H:%M"),
                 RecordedDate = col_date(format = "%m/%d/%Y %H:%M")
               )
)

# Apply data quality filters: remove incomplete responses
US <- US[US$Finished != 0, ]

cat("US Sample Demographic Characteristics:\n")

# Demographic descriptive statistics for US sample
# Gender distribution
cat("\nGender Distribution:\n")
print(table(US$Q2))
print(round(table(US$Q2) / nrow(US) * 100, 2))

# Age distribution
cat("\nAge Distribution:\n")
print(table(US$Q3))
print(round(table(US$Q3) / nrow(US) * 100, 2))

# Ethnic/racial distribution
cat("\nEthnic/Racial Distribution:\n")
print(round(table(US$Q16) / nrow(US) * 100, 2))

# Educational attainment
cat("\nEducational Attainment:\n")
print(table(US$Q17))
print(round(table(US$Q17) / nrow(US) * 100, 2))

# Income distribution
cat("\nIncome Distribution:\n")
print(round(table(US$Q18) / nrow(US) * 100, 2))

# Employment status
cat("\nEmployment Status:\n")
print(table(US$Q19))
print(round(table(US$Q19) / nrow(US) * 100, 2))

# Political affiliation
cat("\nPolitical Affiliation:\n")
print(round(table(US$Q20) / nrow(US) * 100, 2))

# Political orientation
cat("\nPolitical Orientation:\n")
print(round(table(US$Q21) / nrow(US) * 100, 2))

# Extract relevant variables from dataset
US <- US[, -c(1:6, 8:17)]  # Remove extraneous variables
US <- US[, c(33:38, 39:45, 47)]  # Retain Indulgence (33-38) and IBB (39-45, 47) constructs

cat("US data processing completed. Sample size:", nrow(US), "\n")

# -------------------------------------------------------------------------------
# INDIA DATASET PROCESSING
# -------------------------------------------------------------------------------

cat("\n--- India Dataset Processing ---\n")

# Import India dataset with appropriate date format specifications
Ind <- read_csv("India_Survey.csv",
                col_types = cols(
                  StartDate = col_date(format = "%m/%d/%Y %H:%M"),
                  EndDate = col_date(format = "%m/%d/%Y %H:%M"),
                  RecordedDate = col_date(format = "%m/%d/%Y %H:%M")
                )
)

# Apply data quality filters
Ind <- Ind[-which(Ind$Finished == 0), ]  # Remove incomplete responses
Ind <- Ind[Ind$`Q8?Masculinity_3` == 3, ]  # Apply additional quality criteria

cat("India Sample Demographic Characteristics:\n")

# Demographic descriptive statistics for India sample
# Gender distribution
cat("\nGender Distribution:\n")
print(table(Ind$Q2))
print(round(table(Ind$Q2) / nrow(Ind) * 100, 2))

# Age distribution
cat("\nAge Distribution:\n")
print(table(Ind$Q3))
print(round(table(Ind$Q3) / nrow(Ind) * 100, 2))

# Educational attainment
cat("\nEducational Attainment:\n")
print(table(Ind$Q17))
print(round(table(Ind$Q17) / nrow(Ind) * 100, 2))

# Income distribution
cat("\nIncome Distribution:\n")
print(round(table(Ind$Q18) / nrow(Ind) * 100, 2))

# Employment status
cat("\nEmployment Status:\n")
print(table(Ind$Q19))
print(round(table(Ind$Q19) / nrow(Ind) * 100, 2))

# Extract relevant variables from dataset
Ind <- Ind[, -c(1:6, 8:17)]  # Remove extraneous variables
Ind <- Ind[, c(33:38, 39:45, 47)]  # Retain Indulgence (33-38) and IBB (39-45, 47) constructs

cat("India data processing completed. Sample size:", nrow(Ind), "\n")

# -------------------------------------------------------------------------------
# CHINA DATASET PROCESSING
# -------------------------------------------------------------------------------

cat("\n--- China Dataset Processing ---\n")

# Import required library for Excel file processing
library(readxl)

# Import China dataset
CN <- read_excel("China_Survey.xlsx")

# Apply data quality filters
CN <- CN[CN$Q1 != 2, ]  # Remove responses not meeting inclusion criteria

cat("China Sample Demographic Characteristics:\n")

# Demographic descriptive statistics for China sample
# Gender distribution
cat("\nGender Distribution:\n")
print(table(CN$Q2))
print(round(table(CN$Q2) / nrow(CN) * 100, 2))

# Age distribution
cat("\nAge Distribution:\n")
print(table(CN$Q3))
print(round(table(CN$Q3) / nrow(CN) * 100, 2))

# Educational attainment
cat("\nEducational Attainment:\n")
print(table(CN$Q16))
print(round(table(CN$Q16) / nrow(CN) * 100, 2))

# Income distribution
cat("\nIncome Distribution:\n")
print(table(CN$Q17))
print(round(table(CN$Q17) / nrow(CN) * 100, 2))

# Employment status
cat("\nEmployment Status:\n")
print(table(CN$Q18))
print(round(table(CN$Q18) / nrow(CN) * 100, 2))

# Ethnic/racial distribution
cat("\nEthnic/Racial Distribution:\n")
print(table(CN$Q19))
print(round(table(CN$Q19) / nrow(CN) * 100, 2))

# Geographic distribution
cat("\nGeographic Distribution:\n")
print(table(CN$Q20))
print(round(table(CN$Q20) / nrow(CN) * 100, 2))

# Extract relevant variables from dataset
CN <- CN[, -c(1:6)]  # Remove extraneous variables
CN <- CN[, c(31:36, 37:43, 45)]  # Retain Indulgence (31-36) and IBB (37-43, 45) constructs

cat("China data processing completed. Sample size:", nrow(CN), "\n")

# -------------------------------------------------------------------------------
# VARIABLE STANDARDIZATION
# -------------------------------------------------------------------------------

cat("\n--- Variable Standardization ---\n")

# Create consistent variable naming scheme across datasets
namesIDG_IBB <- c(
  paste(rep("IL_IDG", 6), 1:6, sep = ""),  # Indulgence scale items
  paste(rep("IBB", 8), 1:8, sep = "")      # IBB scale items
)

namesIDG <- namesIDG_IBB[1:6]  # Extract Indulgence scale item names

# Standardize variable names across all datasets
names(CN) <- namesIDG_IBB
names(US) <- namesIDG_IBB
names(Ind) <- namesIDG_IBB

cat("Variable naming standardization completed\n")
cat("Total sample sizes - China:", nrow(CN), "India:", nrow(Ind), "US:", nrow(US), "\n")

# ==================================================================================
# PHASE 2: SCALE RELIABILITY ANALYSIS
# ==================================================================================

cat("\n=========================================================\n")
cat("PHASE 2: SCALE RELIABILITY ANALYSIS\n")
cat("=========================================================\n")

# Ensure required packages are available
if (!requireNamespace("psych", quietly = TRUE)) {
  install.packages("psych")
}
library(psych)

# Calculate Cronbach's alpha reliability coefficient for each dataset

# US sample reliability analysis
cat("\n--- US Sample Reliability Analysis ---\n")
alpha_results_US <- alpha(US[, c("IL_IDG1", "IL_IDG2", "IL_IDG3", "IL_IDG4", "IL_IDG5", "IL_IDG6")])
print(alpha_results_US)  # Display reliability statistics

# India sample reliability analysis
cat("\n--- India Sample Reliability Analysis ---\n")
alpha_results_Ind <- alpha(Ind[, c("IL_IDG1", "IL_IDG2", "IL_IDG3", "IL_IDG4", "IL_IDG5", "IL_IDG6")])
print(alpha_results_Ind)  # Display reliability statistics

# China sample reliability analysis
cat("\n--- China Sample Reliability Analysis ---\n")
alpha_results_CN <- alpha(CN[, c("IL_IDG1", "IL_IDG2", "IL_IDG3", "IL_IDG4", "IL_IDG5", "IL_IDG6")])
print(alpha_results_CN)  # Display reliability statistics

# ==================================================================================
# PHASE 3: TVR FRAMEWORK IMPLEMENTATION
# ==================================================================================

cat("\n=========================================================\n")
cat("PHASE 3: TVR FRAMEWORK IMPLEMENTATION\n")
cat("=========================================================\n")

# Ensure required packages are available
if (!requireNamespace("lavaan", quietly = TRUE)) {
  install.packages("lavaan")
}
library(lavaan)

# -------------------------------------------------------------------------------
# STEP 1: TESTING - Initial Measurement Model Evaluation using Multi-Group CFA
# -------------------------------------------------------------------------------

cat("\n--- Step 1: Testing Phase - Original Unidimensional Model Evaluation ---\n")

# Ensure required packages are available
if (!requireNamespace("lavaan", quietly = TRUE)) {
  install.packages("lavaan")
}

# Load lavaan package for SEM analysis
library(lavaan)

# Prepare multi-group dataset for configural invariance testing
cat("\nPreparing multi-group dataset for cross-cultural analysis...\n")

# Duplicate each dataset and add group identifier column 'CI'
CN_multi <- CN    # Copy Chinese dataset
Ind_multi <- Ind  # Copy Indian dataset  
US_multi <- US    # Copy US dataset

# Assign unique group identifiers to each dataset
CN_multi$CI <- 1    # Chinese dataset identifier (CI = 1)
Ind_multi$CI <- 2   # Indian dataset identifier (CI = 2)
US_multi$CI <- 3    # US dataset identifier (CI = 3)

# Combine all datasets into one data frame for multi-group CFA analysis
Full_data <- rbind(CN_multi, Ind_multi, US_multi)

cat("Multi-group dataset completed. Total sample size:", nrow(Full_data), "\n")
cat("Sample sizes - China:", nrow(CN_multi), "India:", nrow(Ind_multi), "US:", nrow(US_multi), "\n")

# Define the baseline measurement model
# Hypothesized structure: Unidimensional Indulgence construct
original_structure <- '
  # Measurement model
  IDG =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4 + IL_IDG5 + IL_IDG6
'

# Fit the configural model for baseline measurement invariance testing
# Configural model allows free estimation across groups to test basic model structure
cat("\nFitting configural model to test baseline measurement invariance...\n")
TestingIDG_configural <- cfa(
  original_structure,           # Specified original measurement model                              
  data = Full_data,            # Combined multi-group dataset
  std.lv = TRUE,               # Standardize latent variables (mean=0, variance=1)
  estimator = "MLR",           # Maximum Likelihood with robust standard errors     
  group = "CI"                 # Specify 'CI' as the grouping variable
)

# Display comprehensive summary of configural model fit
cat("\nConfigural Model Fit Results (Original Unidimensional Structure):\n")
print(summary(
  TestingIDG_configural,
  fit.measures = TRUE,          # Include variety of goodness-of-fit indices
  standardized = TRUE           # Standardized estimates for easier interpretation
))

# -------------------------------------------------------------------------------
# STEP 2: REFINING - Optimizing Factor Structure
# -------------------------------------------------------------------------------

cat("\n--- Step 2: Refining Phase - Factor Structure Optimization ---\n")

# Extract relevant variables for hACO algorithm analysis
# Use existing Full_data from testing step
Data_search <- as.data.frame(Full_data)[,c(1:6,15)]  # Extract Indulgence items (1-6) and group identifier (15)
Full_data <- Data_search
cat("Using multi-group dataset from testing step. Total sample size:", nrow(Full_data), "\n")

# Initialize constraint parameters for hACO algorithm
pre_loaded <- list(c(1),c(6))    # Parameters to be freely estimated
pre_unloaded <- list(c(),c(2))  # Parameters constrained to zero

# Specify fit indices for model evaluation
# 1 = SRMR (Standardized Root Mean Square Residual)
# 4 = TLI (Tucker-Lewis Index)
# 5 = CFI (Comparative Fit Index)
Fit_Use <- c(1, 4, 5) 

cat("\nInitiating hACO algorithm for optimal factor structure search...\n")

# Apply heuristic Ant Colony Optimization algorithm to discover optimal factor structure
results <- hACO(
  Data = Full_data,              
  n_factors = 2,                   # Specify two-factor solution
  n_variables = 6,                 # Number of observed indicators
  loaded = pre_loaded,               # Prior loading constraints (none)
  unloaded = pre_unloaded,           # Prior non-loading constraints (none)
  standards = c(0.92, 1-0.08, 0.92), # Fit thresholds (CFI≥0.92, SRMR≤0.08, TLI≥0.92)
  maxRun = 100,                    # Maximum iterations
  Punish_rate = 0.25,              # Penalty rate for suboptimal solutions
  Fit = Fit_Use,                   # Selected fit indices
  group = "CI"                     # Group identifier for multi-group analysis
)

cat("\nhACO algorithm completed. Fitting optimized measurement model...\n")

# Fit the optimized measurement model from hACO
RefiningIDG_configural <- cfa(
  model = results$OptModelRun,     # Optimal model specification from hACO
  data = Full_data,                # Combined dataset
  std.lv = TRUE,                   # Standardize latent variables
  estimator = "MLR",               # Robust maximum likelihood
  group = "CI"                     # Group identifier
)

# Examine refined model fit
cat("\nOptimized Measurement Model Fit Results:\n")
print(summary(
  RefiningIDG_configural, 
  fit.measures = TRUE, 
  standardized = TRUE
))

# -------------------------------------------------------------------------------
# STEP 3: VALIDATING - Measurement Invariance Testing
# -------------------------------------------------------------------------------

cat("\n--- Step 3: Validating Phase - Measurement Invariance Testing ---\n")

# Define the refined measurement model based on hACO results
# Two-factor solution identified in the Refining step
refined_structure <- '
  # Two-factor measurement model
  IDGF1 =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4  # Factor 1: Items 1-4
  IDGF2 =~ IL_IDG5 + IL_IDG6                      # Factor 2: Items 5-6

  # Factor correlation
  IDGF1 ~~ IDGF2 
'

# Test configural invariance
cat("\nConfigural Invariance Test Results:\n")
ValidatingIDG_configural <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI"
)

print(summary(
  ValidatingIDG_configural,
  fit.measures = TRUE,  
  standardized = TRUE
))

# Test metric invariance
cat("\nMetric Invariance Test Results:\n")
ValidatingIDG_metric <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI",  
  group.equal = "loadings"  # Constrain loadings to be equal
)

print(summary(
  ValidatingIDG_metric,
  fit.measures = TRUE,  
  standardized = TRUE
))

# Test scalar invariance
cat("\nScalar Invariance Test Results:\n")
ValidatingIDG_scalar <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI",  
  group.equal = c("loadings", "intercepts")  # Constrain loadings and intercepts
)

print(summary(
  ValidatingIDG_scalar,
  fit.measures = TRUE,  
  standardized = TRUE
))

# Test strict invariance
cat("\nStrict Invariance Test Results:\n")
ValidatingIDG_strict <- cfa(
  model = refined_structure,  
  data = Full_data,  
  std.lv = TRUE,  
  estimator = "MLR",  
  group = "CI",  
  group.equal = c("loadings", "intercepts", "residuals")  # Constrain all parameters
)

print(summary(
  ValidatingIDG_strict,
  fit.measures = TRUE,  
  standardized = TRUE
))

# -------------------------------------------------------------------------------
# MEASUREMENT INVARIANCE EVALUATION - ΔCFI Analysis
# -------------------------------------------------------------------------------

# Load Function for measurement invariance evaluation using ΔCFI

source("Evaluate_invariance.R")

# Apply invariance evaluation
evaluate_invariance(
  ValidatingIDG_configural,
  ValidatingIDG_metric,
  ValidatingIDG_scalar,
  ValidatingIDG_strict
)

# ==================================================================================
# PHASE 4: EMPIRICAL APPLICATION - STRUCTURAL MODEL TESTING
# ==================================================================================

cat("\n=========================================================\n")
cat("PHASE 4: EMPIRICAL APPLICATION - STRUCTURAL MODEL TESTING\n")
cat("=========================================================\n")

# Test structural model with original unidimensional Indulgence construct
cat("\n--- Original Unidimensional Indulgence Structural Model Results ---\n")
structural_model_Original <- '
  # Measurement models
  IDG =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4 + IL_IDG5 + IL_IDG6  # Single-factor IDG
  IBB =~ IBB1 + IBB2 + IBB3 + IBB4 + IBB5 + IBB6 + IBB7 + IBB8      # IBB construct
  
  # Structural path
  IBB ~ IDG  # Effect of Indulgence on IBB
'

# Fit original structural model
fit_original <- sem(
  model = structural_model_Original,  
  data = Full_data,  
  group = "CI",                      # Multi-group analysis
  estimator = "MLR",                 # Robust estimation
  std.lv = TRUE,                     # Standardize latent variables
  group.equal = c("loadings")        # Metric invariance constraint
)

print(summary(
  fit_original,
  fit.measures = TRUE,  
  standardized = TRUE  
))

# Test structural model with refined two-factor Indulgence construct
cat("\n--- Refined Two-Factor Indulgence Structural Model Results ---\n")
structural_model_Refined <- '
  # Measurement models
  IDGF1 =~ IL_IDG1 + IL_IDG2 + IL_IDG3 + IL_IDG4  # Indulgence Factor 1
  IDGF2 =~ IL_IDG5 + IL_IDG6                      # Indulgence Factor 2
  IBB =~ IBB1 + IBB2 + IBB3 + IBB4 + IBB5 + IBB6 + IBB7 + IBB8  # IBB construct

  # Factor correlation
  IDGF1 ~~ IDGF2  

  # Structural paths
  IBB ~ IDGF1 + IDGF2  # Differential effects of Indulgence factors on IBB
'

# Fit refined structural model
fit_refined <- sem(
  model = structural_model_Refined,  
  data = Full_data,  
  group = "CI",                      # Multi-group analysis
  estimator = "MLR",                 # Robust estimation
  std.lv = TRUE,                     # Standardize latent variables
  group.equal = c("loadings")        # Metric invariance constraint
)

print(summary(
  fit_refined,
  fit.measures = TRUE,  
  standardized = TRUE  
))

cat("\n=========================================================\n")
cat("TVR FRAMEWORK ANALYSIS COMPLETED!\n")
cat("=========================================================\n")

