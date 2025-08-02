########################################################################
#                          hACO_MultiGroup_CFA.R
########################################################################
#
# "Preparation, Patience from Old Drives!"
#
# A comprehensive implementation of the Hybrid Ant Colony Optimization (hACO) 
# algorithm for multi-group Confirmatory Factor Analysis (CFA) model search.
# 
# This script provides functions to automatically discover optimal factor structures
# using the hACO algorithm, with specific support for multi-group analysis. It can
# be used to:
# 1. Find the optimal factor structure for a set of observed variables
# 2. Support multi-group CFA where measurement invariance can be tested across groups
# 3. Provide detailed analysis of model fit and factor loadings
# 4. Generate reports on the optimization process and pheromone distribution
#
# Algorithm overview:
# - The hACO algorithm uses a pheromone-based search strategy inspired by ant colonies
# - Models are generated using probability matrices that evolve over iterations
# - Models with good fit add more pheromone to their parameter configurations
# - The pheromone distribution gradually converges to an optimal model structure
# - The algorithm includes taboo search mechanisms to avoid revisiting poor solutions
# - Multi-group support allows comparing factor structures across different populations
#
# Technical details:
# - Phases: The algorithm operates in three distinct phases:
#   1. Initial search to find a viable model
#   2. Neighborhood search to refine and optimize the model
#   3. Final model selection based on pheromone distribution analysis
# - Adaptive pheromone update with stronger reinforcement for high-quality models
# - Dynamic exploration/exploitation balance through decay rate adjustment
# - Intelligent stagnation detection and recovery with strategic random restarts
# - Early stopping when high-quality solutions stabilize
# - Clean console output with real-time progress tracking
# - Comprehensive ban list management to prevent revisiting poor models
#
# Author: Original by Peipeixi
# Modified for multi-group support with enhanced search efficiency
#
# Usage:
# 1. Source this script: source("hACO_MultiGroup_CFA.R")
# 2. Run the hACO function with your data and parameters:
#    result <- hACO(
#      Data = your_data,        # Data frame with variables (last column = group)
#      n_factors = 2,           # Number of factors to extract
#      n_variables = 6,         # Number of observed variables
#      loaded = list(),         # Prior knowledge about loadings (optional)
#      unloaded = list(),       # Prior knowledge about non-loadings (optional)
#      group = "GroupVar",      # Name of the group variable (last column)
#      maxRun = 100             # Maximum number of iterations (default: 15000)
#    )
#
# Output:
# - OptModel: The optimal model structure in readable format
# - OptModelResults: Detailed fit statistics for the optimal model
# - PheLevel: Final pheromone levels showing variable-factor relationships
# - Additional diagnostic information and alternative models
#
########################################################################

# Required packages
library(MASS)       # For statistical functions
library(lavaan)     # For CFA modeling
library(ShortForm)  # For scale development
library(matrixcalc) # For matrix operations
library(lubridate)  # For date/time functions
library(svMisc)     # For progress indicators
library(stringr)    # For string manipulation

#==============================================================================
#                           UTILITY FUNCTIONS
#==============================================================================
# Helper functions for managing model lists, calculating entropy, and
# other operations used throughout the optimization process.

#' Update Ban List
#'
#' Adds a model to the ban list to avoid re-evaluating poor-fitting models
#'
#' @param BanList Current ban list (or NULL if empty)
#' @param banInd Model indicators to ban
#' @return Updated ban list
BanModel <- function(BanList, banInd) {
  if (missing(BanList)) {
    BanList <- NULL
  }
  BanList <- c(BanList, paste(banInd, collapse = ""))
  return(BanList)
}

#' Update Acceptance List
#' 
#' Adds a model to the acceptance list for tracking good-fitting models
#'
#' @param AccList Current acceptance list (or NULL if empty)
#' @param AccInd Model indicators to accept
#' @return Updated acceptance list
AccModel <- function(AccList, AccInd) {
  if (missing(AccList)) {
    AccList <- NULL
  }
  AccList <- c(AccList, paste(AccInd, collapse = ""))
  return(AccList)
}

#' Build Start Model from Indicators
#'
#' Converts indicator string to numeric vector for starting a new search
#'
#' @param ModInd Model indicator string
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @return Numeric indicator vector
StartBuild <- function(ModInd, n_factors, n_variables) {
  return(as.numeric(unlist(strsplit(ModInd, split = ""))))
}

#' Calculate Entropy
#'
#' Calculates the information entropy of the pheromone distribution
#'
#' @param Phe_level Current pheromone levels
#' @return List containing entropy value
Entropy <- function(Phe_level) {
  entropy <- -sum(Phe_level/sum(Phe_level) * log2(Phe_level/sum(Phe_level)))
  return(list(Entropy = entropy))
}

#' Display Pheromone Distribution
#'
#' Creates a summary of the pheromone distribution for analysis
#'
#' @param Phe_level Current pheromone levels
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param Data Dataset containing variables (for variable names)
#' @param top_n Number of top combinations to display
#' @return Sorted data frame of pheromone distribution (invisible)
DisplayPheromoneDistribution <- function(Phe_level, n_factors, n_variables, Data, top_n = 10) {
  # Create a data frame with indices
  phe_df <- data.frame(
    FactorVar = paste0("F", rep(1:n_factors, each = n_variables), "_V", rep(1:n_variables, n_factors)),
    Pheromone = Phe_level,
    VarIndex = rep(1:n_variables, n_factors),
    FactorIndex = rep(1:n_factors, each = n_variables)
  )
  
  # Add variable names if data has column names
  if (!is.null(Data) && !is.null(names(Data))) {
    var_names <- names(Data)
    phe_df$VarName <- var_names[phe_df$VarIndex]
  }
  
  # Sort by pheromone value
  phe_df_sorted <- phe_df[order(phe_df$Pheromone, decreasing = TRUE), ]
  
  # Display top N combinations with highest pheromone
  cat("\nTop", min(top_n, nrow(phe_df_sorted)), "highest pheromone combinations:\n")
  print(head(phe_df_sorted, top_n), row.names = FALSE)
  
  # Calculate average pheromone by factor
  cat("\nAverage pheromone by factor:\n")
  for (f in 1:n_factors) {
    avg_phe <- mean(phe_df$Pheromone[phe_df$FactorIndex == f])
    cat("Factor", f, ":", sprintf("%.5f", avg_phe), "\n")
  }
  
  return(invisible(phe_df_sorted))
}

#' Display Progress Information
#'
#' Displays detailed progress information during algorithm execution
#' with clean, non-cluttering output by overwriting the previous line.
#' This function handles the real-time status updates during the search process,
#' ensuring that the console remains clean and only shows the most current information.
#'
#' The progress display includes:
#' - Percentage completion with a visual progress bar
#' - Current iteration number and total iterations
#' - Number of models evaluated so far
#' - Number of accepted models (high quality)
#' - Maximum pheromone value found (quality metric)
#' - Optional status message for additional context
#'
#' @param numTotal Current iteration number
#' @param maxRun Maximum number of iterations
#' @param numRuns Number of models evaluated
#' @param numAcc Number of accepted models
#' @param maxPheromone Maximum pheromone value found
#' @param newline Whether to add a newline after progress information
#' @param status Optional status to display after the progress information
#' @return NULL (invisible)
DisplayProgress <- function(numTotal, maxRun, numRuns = NULL, numAcc = NULL, 
                           maxPheromone = NULL, newline = FALSE, status = NULL) {
  # Minimal clearing, using only carriage return
  cat("\r")
  
  # Basic progress information
  percent <- sprintf("%3d%%", as.integer(100 * numTotal / maxRun))
  width <- 40  # Reduced width to leave space for status messages
  n_chars <- as.integer(width * numTotal / maxRun)
  progress_bar <- paste0(
    "[", 
    paste(rep("=", n_chars), collapse = ""),
    paste(rep(" ", width - n_chars), collapse = ""),
    "]"
  )
  
  # Create basic progress string
  progress_info <- paste0(percent, " ", progress_bar, " ", numTotal, "/", maxRun)
  
  # Add detailed information if requested and available
  if (!is.null(numRuns) && !is.null(numAcc) && !is.null(maxPheromone)) {
    detail_text <- sprintf(" | Models: %d | Accepted: %d | Max Ph: %.4f", 
                           numRuns, numAcc, maxPheromone)
    progress_info <- paste0(progress_info, detail_text)
  }
  
  # Display progress
  cat(progress_info, sep="")
  utils::flush.console() # Ensure immediate display
  
  # Add status information after the progress bar if provided
  if (!is.null(status)) {
    # Add status on same line after some space
    cat("  ", status, sep="")
    utils::flush.console() # Ensure display
  }
  
  # Add final newline if requested
  if (newline) {
    cat("\n")
    utils::flush.console()
  }
  
  return(invisible(NULL))
}

#' Display Warning Without Accumulating Output
#'
#' @param message Warning message to display
#' @return NULL (invisible)
DisplayWarning <- function(message) {
  # Simplified clearing method
  utils::flush.console()
  cat("\r", strrep(" ", 120), "\r", sep="")
  cat(message, "\n", sep="")
  utils::flush.console()
  return(invisible(NULL))
}

#==============================================================================
#                         MODEL BUILDING FUNCTIONS
#==============================================================================
# Functions for creating and manipulating CFA model structures, converting
# between probability matrices, indicator vectors, and lavaan syntax.

#' Build Initial Probability Matrix
#'
#' Creates a probability matrix that incorporates prior knowledge about which variables
#' should or should not load on specific factors. This function sets up the initial
#' conditions for the hACO algorithm's search process.
#' 
#' The probability matrix determines how likely a variable is to load on each factor
#' during the model generation process. Values close to 1 indicate high probability
#' of loading, while values close to 0 indicate low probability.
#' 
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param loaded List of variables that should definitely load on each factor (prior knowledge)
#'        Format: Either empty list() for no constraints, or list with exactly n_factors elements
#'        Each element should contain variable indices or be empty c()
#' @param unloaded List of variables that should definitely not load on each factor
#'        Format: Either empty list() for no constraints, or list with exactly n_factors elements  
#'        Each element should contain variable indices or be empty c()
#' @return Probability matrix (dimensions: n_factors * n_variables) for model building
ProbBuild <- function(n_factors, n_variables, loaded, unloaded) {
  # Create base matrices
  baseMatrix <- matrix(rep(0, n_variables*n_factors), n_factors, n_variables)
  loadMatrix <- baseMatrix
  
  # Handle empty lists and validate input
  # If loaded is completely empty, create empty list for each factor
  if (length(loaded) == 0) {
    loaded <- vector("list", n_factors)
  } else if (length(loaded) != n_factors) {
    stop("loaded must be either empty list() or have exactly n_factors elements", call. = FALSE)
  }
  
  # If unloaded is completely empty, create empty list for each factor  
  if (length(unloaded) == 0) {
    unloaded <- vector("list", n_factors)
  } else if (length(unloaded) != n_factors) {
    stop("unloaded must be either empty list() or have exactly n_factors elements", call. = FALSE)
  }
  
  # Fill in known loadings (setting to 1 indicates a forced loading)
  for (i in 1:n_factors) {
    if (length(loaded[[i]]) > 0) {
      # Validate variable indices
      var_indices <- as.numeric(loaded[[i]])
      if (any(var_indices < 1 | var_indices > n_variables)) {
        stop(paste("Invalid variable index in loaded[[", i, "]]. All indices must be between 1 and", n_variables), call. = FALSE)
      }
      loadMatrix[i, var_indices] <- 1
    }
  }
  
  # Fill in known non-loadings (setting to 1 indicates a prohibited loading)
  unloadMatrix <- baseMatrix
  for (i in 1:n_factors) {
    if (length(unloaded[[i]]) > 0) {
      # Validate variable indices
      var_indices <- as.numeric(unloaded[[i]])
      if (any(var_indices < 1 | var_indices > n_variables)) {
        stop(paste("Invalid variable index in unloaded[[", i, "]]. All indices must be between 1 and", n_variables), call. = FALSE)
      }
      unloadMatrix[i, var_indices] <- 1
    }
  }
  
  # Check for invalid configurations
  if (is.element(n_factors, colSums(unloadMatrix))) {
    stop("External misspecification: A variable cannot be prohibited from loading on all factors", call. = FALSE)
  }
  
  if (is.element(2, loadMatrix+unloadMatrix)) {
    stop("Contradiction in settings: A variable cannot be both forced to load and prohibited", call. = FALSE)
  }  
  
  # Create the probability matrix with default probability of 0.55
  # This balanced value gives each loading an approximately equal chance initially
  pre_loads <- list(load = loadMatrix, unload = unloadMatrix)
  prob <- matrix(rep(0.55, n_variables*n_factors), n_factors, n_variables, byrow = TRUE)
  prob[which(pre_loads$load == 1)] <- 1.0  # Force loadings based on prior knowledge
  prob[which(pre_loads$unload == 1)] <- 0.0  # Prohibit loadings based on prior knowledge
  
  return(prob)
} 

#' Build Model Indicators
#'
#' Creates binary indicator vectors for model construction based on probabilities.
#' Each indicator represents whether a variable loads on a factor (1) or not (0).
#' 
#' This function performs the random sampling that generates potential models,
#' and includes a repair mechanism to ensure that every variable loads on at least
#' one factor, which is necessary for model identification.
#' 
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param probability Matrix of probabilities for variable loadings (n_factors * n_variables)
#' @return List containing:
#'   \item{test}{Flattened indicator vector for the model}
#'   \item{ind_use}{Matrix form of indicators (n_factors * n_variables)}
IndBuild <- function(n_factors, n_variables, probability) {
  # Generate random binary values based on probabilities
  # This is the core randomization step in the hACO algorithm
  ind <- rbinom(n_variables*n_factors, 1, as.vector(t(probability)))
  ind_use <- matrix(ind, n_factors, n_variables, byrow = TRUE)
  
  # Repair model to ensure every variable loads on at least one factor
  # This is necessary for model identification in CFA
  repair_factor <- n_factors
  while (is.element(0, colSums(ind_use))) {
    repair <- which(colSums(ind_use) == 0)  # Find variables with no loadings
    repair_line <- which(probability[repair_factor,] != 0)  # Find eligible positions
    repair_element <- repair_line[which(repair_line %in% repair)]  # Select variables to repair
    ind_use[repair_factor,][repair_element] <- 1  # Add loading
    repair_factor <- repair_factor - 1  # Try different factor if needed next time
  }
  
  # Convert back to vector form
  ind <- c(t(ind_use))
  return(list(test = ind, ind_use = ind_use))
}

#' Build Complete CFA Model
#'
#' Constructs a complete CFA model with proper lavaan syntax based on probability matrices.
#' This function handles the full model generation process, including checking against
#' banned configurations and formatting the model for lavaan.
#' 
#' The function ensures that generated models:
#' 1. Have appropriate factor loadings based on probability matrices
#' 2. Are not in the list of previously banned (poor-fitting) models
#' 3. Are formatted correctly for estimation with lavaan
#' 
#' @param Data Dataset containing observed variables
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param probability Matrix of probabilities for variable loadings
#' @param BanList List of previously banned model configurations
#' @param CutChange Counter for unsuccessful attempts
#' @param ChangeRun Maximum attempts before resetting ban list
#' @return List containing:
#'   \item{ModelInd}{Indicator vector for the model}
#'   \item{CFAModel}{Lavaan syntax for the CFA model}
#'   \item{Change}{Flag indicating if ban list was reset}
ModelBuilder <- function(Data, n_factors, n_variables, probability, BanList, CutChange = 0, ChangeRun = 1000){
  Data <- as.data.frame(Data)
  
  # Build model indicators
  ind <- IndBuild(n_factors, n_variables, probability)
  model_preTest <- paste(ind$test, collapse="")
  Change = FALSE
  
  # Keep generating new models until one not in ban list is found
  # This ensures we don't re-evaluate models that were already found to be poor-fitting
  while ((is.element(model_preTest, BanList))) {
    ind <- IndBuild(n_factors, n_variables, probability)
    model_preTest <- paste(ind$test, collapse="")
    CutChange <- CutChange + 1
    if (CutChange > ChangeRun){
      # If we've tried too many times, reset the ban list to avoid getting stuck
      Change <- TRUE
      BanList <- NULL
    }
  }
  
  # Get variable names from data
  x <- names(Data)
  
  # Create factor names (F1, F2, etc.)
  f <- rep("F", n_factors)
  f <- paste(f, c(1:n_factors), sep = "")
  
  # Build model syntax in lavaan format
  mod <- NULL
  for (i in 1:n_factors) {
    # Format: F1 =~ x1 + x2 + x3
    mod_sub <- paste(f[i], "=~", paste(x[which(ind$ind_use[i,] == 1)], collapse = "+"), sep = "")
    mod <- c(mod, mod_sub)
  }
  model_cfa <- paste(mod, collapse = "\n                     ")
  
  return(list(ModelInd = ind$test, CFAModel = model_cfa, Change = Change))
}

#' Rebuild Model from Indicators
#'
#' Reconstructs a CFA model from a model indicator string. This function is used
#' to rebuild models from stored indicators, typically when recovering the best
#' models from a previous run or analyzing accepted models.
#' 
#' This function is the inverse of the model building process, taking a string
#' of indicators and converting it back to a proper lavaan model syntax.
#' 
#' @param ModInd Model indicator string (binary sequence)
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param Data Dataset containing variables (optional, used for variable names)
#' @param group Group variable name for multi-group analysis (optional)
#' @return List containing:
#'   \item{ModelInd}{Numeric indicator vector for the model}
#'   \item{CFAModel}{Lavaan syntax for the CFA model}
ModelReBuild <- function(ModInd, n_factors, n_variables, Data = NULL, group = NULL) {
  # Convert string to matrix of indicators
  ind <- matrix(as.numeric(unlist(strsplit(ModInd, split = ""))), 
                n_factors, n_variables, byrow = TRUE)
  model_ind <- as.numeric(unlist(strsplit(ModInd, split = "")))
  
  # Generate variable names - either from data or default X1, X2, etc.
  if (is.null(names(Data))) {
    x <- rep("X", n_variables)
    x <- paste(x, c(1:n_variables), sep = "")
  } else {
    x <- names(Data)
    if (!is.null(group) && group %in% names(Data)) {
      x <- x[x != group]  # Exclude group variable from the model variables
    }
  }
  
  # Generate factor names (F1, F2, etc.)
  f <- rep("F", n_factors)
  f <- paste(f, c(1:n_factors), sep = "")
  
  # Build model syntax in lavaan format
  mod <- NULL
  for (i in 1:n_factors) {
    mod_sub <- paste(f[i], "=~", paste(x[which(ind[i,] == 1)], collapse = "+"), sep = "")
    mod <- c(mod, mod_sub)
  }
  model_cfa <- paste(mod, collapse = "\n                     ")
  
  return(list(ModelInd = model_ind, CFAModel = model_cfa))
}

#' Build Optimal Model from Pheromone Levels
#'
#' Creates a model based on the distribution of pheromone levels.
#' This function applies a threshold to the pheromone distribution to determine
#' which variables should load on each factor in the final model.
#' 
#' The function works by:
#' 1. Applying a threshold to the pheromone levels based on the rate parameter
#' 2. Converting thresholded values to binary indicators (1=load, 0=don't load)
#' 3. Building the corresponding lavaan syntax for the CFA model
#' 4. Optionally displaying debug information about the model structure
#' 
#' @param Data Dataset containing variables
#' @param PheLevel Current pheromone levels
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param rate Threshold rate for including a variable (proportion of max pheromone)
#' @param group Group variable name (if multi-group analysis)
#' @param verbose Whether to print debug information about the model structure (default: FALSE)
#' @return List containing the model and its components
#'   \item{CFAModel}{Lavaan syntax for the CFA model}
#'   \item{OptModelShow}{Simplified model display string}
#'   \item{ModelInd}{Binary indicators for the model}
#'   \item{ModelPreTest}{Matrix form of indicators (factors * variables)}
OptModel <- function(Data, PheLevel, n_factors, n_variables, rate, group = NULL, verbose = FALSE) {
  # Threshold the pheromone levels to create binary indicators
  PheLevel <- signif(PheLevel, 3)
  PheLevel[which(PheLevel <= rate * max(PheLevel))] <- 0
  PheLevel[which(PheLevel != 0)] <- 1
  Ind_OptModel <- matrix(PheLevel, n_factors, n_variables, byrow = TRUE)
  
  # Get variable names, excluding group variable if specified
  x <- names(Data)
  if (!is.null(group) && group %in% names(Data)) {
    x <- x[x != group]  # Exclude group variable
  }
  
  # Create factor names
  f <- rep("F", n_factors)
  f <- paste(f, c(1:n_factors), sep = "")
  
  # Build model syntax
  Optmod <- NULL
  for (i in 1:n_factors) {
    Optmod_sub <- paste(f[i], "=~", paste(x[which(Ind_OptModel[i, ] == 1)], collapse = "+"), sep = "")
    Optmod <- c(Optmod, Optmod_sub)
  }
  OptModelRun <- paste(Optmod, collapse = "\n                     ")
  OptModelShow <- paste(Optmod, collapse = "   ")
  
  # Only show debug information when verbose=TRUE
  if (verbose) {
    cat("PheLevel max:", max(PheLevel), "\n")
    cat("Factors with variables:\n")
    for (i in 1:n_factors) {
      cat("F", i, ": ", paste(which(Ind_OptModel[i,] == 1), collapse=", "), "\n", sep="")
    }
  }
  
  return(list(
    CFAModel = OptModelRun, 
    OptModelShow = OptModelShow, 
    ModelInd = PheLevel, 
    ModelPreTest = Ind_OptModel
  ))
}

#' Build Model from Beyond Average Method
#'
#' Alternative model building approach using mean pheromone level as threshold
#' 
#' @param PheLevel Current pheromone levels
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param Data Dataset containing variables
#' @param group Group variable name (if multi-group analysis)
#' @return List with model syntax
BeyondModel <- function(PheLevel, n_factors, n_variables, Data = NULL, group = NULL){
  # Round pheromone levels and apply mean threshold
  PheLevel <- round(PheLevel, 4)
  PheLevel[which(PheLevel < sum(PheLevel)/(n_factors*n_variables + 1))] <- 0
  PheLevel[which(PheLevel != 0)] <- 1
  Ind_BeyondModel <- matrix(PheLevel, n_factors, n_variables, byrow = TRUE)
  
  # Generate variable names
  if (is.null(names(Data))){
    x <- rep("X", n_variables)
    x <- paste(x, c(1:n_variables), sep = "")
  } else {
    x <- names(Data)
    if(!is.null(group) && group %in% names(Data)) {
      x <- x[x != group]  # Exclude group variable
    }
  }
  
  # Generate factor names
  f <- rep("F", n_factors)
  f <- paste(f, c(1:n_factors), sep = "")
  
  # Build model syntax
  Beyondmod <- NULL
  for (i in 1:n_factors) {
    Beyondmod_sub <- paste(f[i], "=~", paste(x[which(Ind_BeyondModel[i,] == 1)], collapse = "+"), sep = "")
    Beyondmod <- c(Beyondmod, Beyondmod_sub)
  }
  BeyondModelRun <- paste(Beyondmod, collapse = "\n                     ")
  BeyondModelShow <- paste(Beyondmod, collapse = "   ")
  
  return(list(BeyondModelRun = BeyondModelRun, BeyondModelShow = BeyondModelShow))
}

#==============================================================================
#                        MODEL EVALUATION FUNCTIONS
#==============================================================================
# Functions for estimating CFA models and calculating fit indices.
# These interface with lavaan to evaluate models and extract fit measures.

#' Model Estimation Function
#'
#' Estimates a CFA model using lavaan and returns specified fit indices.
#' This function handles the actual model fitting and extracts the relevant
#' fit measures used to evaluate model quality.
#' 
#' The function also checks for negative loadings, which are generally undesirable
#' in CFA models and can indicate model misspecification.
#' 
#' @param Model Model structure containing CFAModel field with lavaan syntax
#' @param Data Dataset containing variables for model estimation
#' @param Fit Vector of fit indices to calculate (1=srmr, 2=rmsea, 3=ifi, 4=tli, 5=cfi, 6=mfi)
#' @param group Group variable name for multi-group analysis (optional)
#' @return List with:
#'   \item{Fit1}{First fit measure value}
#'   \item{Fit2}{Second fit measure value} 
#'   \item{Fit3}{Third fit measure value}
#'   \item{Modelres}{Full lavaan model object}
#'   \item{BIC}{Bayesian Information Criterion value}
ModelEst <- function(Model, Data, Fit, group = NULL) {
  options(warn = 2)  # Convert warnings to errors to catch estimation problems
  
  # Fit CFA model with or without group variable
  if (is.null(group)) {
    # Single-group CFA
    fit_cfa <- try(cfa(
      Model$CFAModel,
      data = Data,
      std.lv = TRUE,      # Standardize latent variables
      estimator = "MLR"   # Robust maximum likelihood estimator
    ), silent = TRUE)
  } else {
    # Multi-group CFA
    fit_cfa <- try(cfa(
      Model$CFAModel,
      data = Data,
      std.lv = TRUE,
      estimator = "MLR",
      group = group      # Use group variable for multi-group analysis
    ), silent = TRUE)
  }
  
  # Check for negative loadings - these are usually problematic in CFA
  if (!inherits(fit_cfa, "try-error")) {
    params <- try(parameterEstimates(fit_cfa), silent = TRUE)
    
    if (!inherits(params, "try-error") && is.data.frame(params)) {
      loadings_params <- params[params$op == "=~", ]  # Extract factor loadings
      
      if (any(loadings_params$est < 0)) {
        # If negative loadings exist, return poor fit values
        Fit1 <- 0.70
        Fit2 <- 0.70
        Fit3 <- 0.70
        Fits <- cbind(Fit1, Fit2, Fit3)
        bic <- Inf
        options(warn = 0)  # Reset warning behavior
        return(list(Fit1 = Fits[1], Fit2 = Fits[2], Fit3 = Fits[3], 
                    Modelres = fit_cfa, BIC = bic))
      }
    }
  }
  
  # Extract fit measures based on requested indices
  Fit_In <- c("srmr", "rmsea", "ifi", "tli", "cfi", "mfi")
  Fit1 <- try(fitMeasures(fit_cfa, Fit_In[Fit[1]]), silent = TRUE)
  Fit2 <- try(fitMeasures(fit_cfa, Fit_In[Fit[2]]), silent = TRUE)
  Fit3 <- try(fitMeasures(fit_cfa, Fit_In[Fit[3]]), silent = TRUE)
  bic <- try(BIC(fit_cfa), silent = TRUE)
  Fits <- cbind(Fit1, Fit2, Fit3)  
  
  # Handle errors in fit measures
  if (class(Fit1)[1] == "try-error" | class(Fit2)[1] == "try-error" |
      class(Fit3)[1] == "try-error" | class(bic)[1] == "try-error") {
    
    if (length(which(Fit < 3)) > 0) {
      # For indices where lower values are better (SRMR, RMSEA)
      Fits[which(Fit < 3)] <- 1
      Fits[which(Fit >= 3)] <- 0
      Fits <- as.numeric(Fits)
      bic <- Inf
    } else {
      Fits <- 0
      bic <- -Inf
    }
  }
  
  # Process fit values - clip values above 1
  if (length(which(Fits > 1)) > 0) {
    Fits[which(Fits > 1)] <- 1
  }
  
  # For indices where lower values are better (SRMR, RMSEA), invert the score
  # so that higher values always represent better fit
  if (length(which(Fit < 3)) > 0) {
    Fits[which(Fit < 3)] <- 1 - Fits[which(Fit < 3)]
  }
  
  options(warn = 0)  # Reset warning behavior
  return(list(Fit1 = Fits[1], Fit2 = Fits[2], Fit3 = Fits[3], 
              Modelres = fit_cfa, BIC = bic))
}

#' Optimal Model Estimation
#'
#' Enhanced model estimation function that also calculates a pheromone value
#' based on the fit indices. The pheromone value represents the overall quality
#' of the model and is used to guide the hACO optimization process.
#' 
#' This function uses a more sophisticated approach to evaluate models compared
#' to the basic ModelEst function, including applying discrimination parameters
#' to each fit index and calculating a complexity penalty.
#' 
#' @param Model Model structure with lavaan syntax
#' @param Data Dataset containing variables
#' @param Fit Vector of fit indices to calculate (1=srmr, 2=rmsea, 3=ifi, 4=tli, 5=cfi, 6=mfi)
#' @param dicrimination Vector of discrimination parameters for each fit index
#' @param standards Vector of standard/target values for each fit index
#' @param scopes Vector of scope parameters for sensitivity of each fit index
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param Punish_rate Rate for penalizing model complexity
#' @param group Group variable name for multi-group analysis (optional)
#' @return List with:
#'   \item{Fit1}{First fit measure value}
#'   \item{Fit2}{Second fit measure value}
#'   \item{Fit3}{Third fit measure value}
#'   \item{Modelres}{Full lavaan model object}
#'   \item{Ph}{Calculated pheromone value for the model}
OptModelEst <- function(Model, Data, Fit, dicrimination = c(1, 1, 1),
                        standards = c(0.93, 1-0.07, 0.93),
                        scopes = c(0.015, 0.015, 0.015), 
                        n_factors, n_variables, Punish_rate, group = NULL) {
  options(warn = 2)  # Convert warnings to errors to catch estimation problems
  
  # Fit CFA model with or without group variable
  if(is.null(group)) {
    fit_cfa <- try(cfa(
      Model$CFAModel,
      data = Data,
      std.lv = TRUE,
      estimator = "MLR"
    ), silent = TRUE)
  } else {
    fit_cfa <- try(cfa(
      Model$CFAModel,
      data = Data,
      std.lv = TRUE,
      estimator = "MLR",
      group = group
    ), silent = TRUE)
  }
  
  # Check for negative loadings - these are usually problematic in CFA
  if (!inherits(fit_cfa, "try-error")) {
    params <- try(parameterEstimates(fit_cfa), silent = TRUE)
    
    if (!inherits(params, "try-error") && is.data.frame(params)) {
      loadings_params <- params[params$op == "=~", ]
      
      if (any(loadings_params$est < 0)) {
        # If negative loadings exist, return poor fit values and zero pheromone
        if (length(which(Fit < 3)) > 0){
          Fits <- rep(0.70, 3)
          Fits[which(Fit < 3)] <- 1 - Fits[which(Fit < 3)]
          Fits <- as.numeric(Fits)
          bic <- Inf
          ph <- 0
        } else {
          Fits <- rep(0.70, 3)
          bic <- -Inf
          ph <- 0
        }
        options(warn = 0)
        return(list(
          Fit1 = Fits[1], Fit2 = Fits[2], Fit3 = Fits[3],
          Modelres = fit_cfa,
          Ph = ph
        ))
      }
    }
  }
  
  # Extract fit measures based on requested indices
  Fit_In <- c("srmr", "rmsea", "ifi", "tli", "cfi", "mfi")
  Fit1 <- try(fitMeasures(fit_cfa, Fit_In[Fit[1]]), silent = TRUE)
  Fit2 <- try(fitMeasures(fit_cfa, Fit_In[Fit[2]]), silent = TRUE)
  Fit3 <- try(fitMeasures(fit_cfa, Fit_In[Fit[3]]), silent = TRUE)
  bic <- try(BIC(fit_cfa), silent = TRUE)
  Fits <- cbind(Fit1, Fit2, Fit3)  
  
  # Handle errors in fit measures
  if (class(Fit1)[1] == "try-error" | class(Fit2)[1] == "try-error" | 
      class(Fit3)[1] == "try-error" | class(bic)[1] == "try-error"){
    
    if (length(which(Fit < 3)) > 0){
      Fits[which(Fit < 3)] <- 1
      Fits[which(Fit >= 3)] <- 0
      Fits <- as.numeric(Fits)
      bic <- Inf
      ph <- 0
    } else {
      Fits <- 0
      bic <- -Inf
      ph <- 0
    }
  } else {
    # Calculate pheromone value using logistic function
    fit_temp <- Fits
    
    if (length(which(fit_temp > 1)) > 0) {
      fit_temp[which(fit_temp > 1)] <- 1
    }
    
    if (length(which(Fit < 3)) > 0) {
      fit_temp[which(Fit < 3)] <- 1 - fit_temp[which(Fit < 3)]
    }
    
    # Minimum fit thresholds before considering a model
    if (fit_temp[1] < 0.80 | fit_temp[2] < 0.80 | fit_temp[3] < 0.80) {
      ph <- 0
    } else {
      # Calculate pheromone as product of logistic transformations of each fit index
      # This creates a smooth, non-linear response to fit improvement
      ph <- 1
      for (i in 1:(length(fit_temp))) {
        ph_sub <- (exp((
          dicrimination[i] *
            (as.numeric(fit_temp[i]) - standards[i])
        ) / scopes[i])) /
          (1 + exp((
            dicrimination[i] * (as.numeric(fit_temp[i]) - standards[i])
          ) / scopes[i]))
        ph <- ph * ph_sub
      }
      
      # Apply complexity punishment for models with too many parameters
      if (sum(Model$ModelInd) < n_variables) {
        punish <- 0
      } else {
        punish <- Punish_rate * (log(1 + (sum(Model$ModelInd) - n_variables) / 
                                    ((n_factors - 1) * n_variables)))
      }
      
      ph <- ph - punish
      if (ph < 0) {
        ph <- 0
      }
    }
    
    # Process fit values - clip values above 1
    if (length(which(Fits > 1)) > 0) {
      Fits[which(Fits > 1)] <- 1
    }
    
    # For indices where lower values are better (SRMR, RMSEA), invert the score
    if (length(which(Fit < 3)) > 0) {
      Fits[which(Fit < 3)] <- 1 - Fits[which(Fit < 3)]
    }
  }
  
  options(warn = 0)  # Reset warning behavior
  return(list(
    Fit1 = Fits[1], Fit2 = Fits[2], Fit3 = Fits[3],
    Modelres = fit_cfa,
    Ph = ph
  ))
}

#' Complexity Punishment Function
#'
#' Calculates a punishment value for model complexity
#' 
#' @param Model Model structure
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param Punish_rate Rate for punishing model complexity
#' @return Punishment value
ComPunish <- function(Model, n_factors, n_variables, Punish_rate){
  if (sum(Model$ModelInd) < n_variables){
    punish <- 0
  } else {
    punish <- Punish_rate * (log(1 + (sum(Model$ModelInd) - n_variables) / 
                                ((n_factors - 1) * n_variables)))
  }
  return(punish)
}

#' Pheromone Calculation Function
#'
#' Calculates the pheromone value for a model
#' 
#' @param Model Model structure
#' @param Data Dataset containing variables
#' @param dicrimination Discrimination parameters
#' @param standards Standard values for fit indices
#' @param scopes Scope parameters
#' @param BanList Current ban list
#' @param numRuns Counter for model runs
#' @param LowCut Lower threshold for accepting a model
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param Punish_rate Rate for punishing model complexity
#' @param Fit Vector of fit indices to calculate
#' @param group Group variable name (if multi-group analysis)
#' @return List with pheromone value and updated ban list
Pheromone <- function(Model, Data, dicrimination = c(1, 1, 1),
                      standards = c(0.93, 1-0.07, 0.93),
                      scopes = c(0.015, 0.015, 0.015), 
                      BanList, numRuns = 0, LowCut = 0.001, 
                      n_factors, n_variables, Punish_rate, Fit, group = NULL){
  
  # Estimate model
  model_estimation <- ModelEst(Model, Data, Fit, group = group)
  BanInd <- FALSE
  
  # Check if model fits poorly
  if (model_estimation$Fit1 < 0.80 | model_estimation$Fit2 < 0.80 | model_estimation$Fit3 < 0.80){
    BanList <- BanModel(BanList, Model$ModelInd)
    ph <- 0
    BanInd <- TRUE
  } else {
    # Calculate pheromone value
    ph <- 1
    for (i in 1:(length(model_estimation) - 2)){
      ph_sub <- (exp((dicrimination[i] * 
                      (as.numeric(model_estimation[i]) - standards[i])) / scopes[i])) / 
        (1 + exp((dicrimination[i] * (as.numeric(model_estimation[i]) - standards[i])) / scopes[i]))
      ph <- ph * ph_sub
    }
    
    # Apply complexity punishment
    ph <- ph - ComPunish(Model = Model, n_factors = n_factors, 
                         n_variables = n_variables, Punish_rate = Punish_rate)
    if (ph < 0){
      ph <- 0
    }
    numRuns <- numRuns + 1
  }
  
  # Ban models with low pheromone
  if ((ph < LowCut) & (BanInd == FALSE)){
    BanList <- BanModel(BanList, Model$ModelInd)
  }
  
  return(list(Pheromone = ph, BanList = BanList, numRuns = numRuns, BIC = model_estimation$BIC))
}

#' Estimate Model from AccList
#'
#' Estimates a model from the acceptance list
#' 
#' @param ModInd Model indicator string
#' @param Data Dataset containing variables
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param Fit Vector of fit indices to calculate
#' @param group Group variable name (if multi-group analysis)
#' @return List with fit measures and model object
ModelIndEst <- function(ModInd, Data, n_factors, n_variables, Fit = c(5, 2, 3), group = NULL) {
  model <- ModelReBuild(ModInd = ModInd, n_factors = n_factors, 
                        n_variables = n_variables, Data = Data, group = group)
  resModelIndEst <- ModelEst(Model = model, Data = Data, Fit = Fit, group = group)
  return(list(Fit1 = resModelIndEst$Fit1, Fit2 = resModelIndEst$Fit2, 
              Fit3 = resModelIndEst$Fit3, Modelres = resModelIndEst$Modelres))
}

#' Estimate Beyond Average Model
#'
#' Estimates a model using the Beyond Average method
#' 
#' @param PheLevel Current pheromone levels
#' @param Data Dataset containing variables
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param Fit Vector of fit indices to calculate
#' @param group Group variable name (if multi-group analysis)
#' @return List with fit measures and model object
BeyondModelEst <- function(PheLevel, Data, n_factors, n_variables, Fit = c(5, 2, 3), group = NULL){
  beyondmodel <- BeyondModel(PheLevel = PheLevel, n_factors = n_factors, 
                             n_variables = n_variables, Data = Data, group = group)
  beyondmodelres <- OptModelEst(
    Model = list(CFAModel = beyondmodel$BeyondModelRun, ModelInd = PheLevel), 
    Data = Data, 
    Fit = Fit, 
    n_factors = n_factors, 
    n_variables = n_variables, 
    Punish_rate = 0.1, 
    group = group
  )
  return(list(Fit1 = beyondmodelres$Fit1, Fit2 = beyondmodelres$Fit2, 
              Fit3 = beyondmodelres$Fit3, Modelres = beyondmodelres$Modelres))
}

#==============================================================================
#                       PHEROMONE UPDATE FUNCTIONS
#==============================================================================
# Core mechanisms of the hACO algorithm that guide the search toward optimal
# solutions by updating pheromone levels based on model performance.

#' Update Probability Based on Pheromone
#'
#' Adjusts model building probabilities based on pheromone levels
#' 
#' @param Model Current model
#' @param Pheromone Pheromone value for the model
#' @param AccList Current acceptance list
#' @param probability Current probability matrix
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param P_dicrimination Discrimination parameter for probability
#' @param P_standards Standard value for probability
#' @param P_scopes Scope parameter for probability
#' @param HighCut Upper threshold for accepting a model
#' @return List with updated probability matrix and acceptance list
ProbCompensation <- function(Model, Pheromone, AccList, probability, n_factors, n_variables,
                           P_dicrimination = 0.5, P_standards = 0.95, P_scopes = 0.005, HighCut = 0.75){
  if (Pheromone > HighCut){
    # Add model to acceptance list
    AccList <- AccModel(AccList, Model$ModelInd)
    
    if (length(AccList) > 5){
      # Calculate probability increase
      add_rate <- (exp((P_dicrimination * 
                       (Pheromone - P_standards)) / P_scopes)) / 
        (1 + exp((P_dicrimination * (Pheromone - P_standards)) / P_scopes))
      
      # Apply increase to probability matrix
      prob_add <- matrix((Model$ModelInd * add_rate / 5), n_factors, n_variables) 
      prob <- probability + prob_add^20
      prob[which(prob > 1)] <- 1
      
      # Maintain original constraints from ProbBuild
      # Force loadings (probability=1) must remain 1.0
      prob[which(probability == 1)] <- 1.0
      # Prohibit loadings (probability=0) must remain 0.0  
      prob[which(probability == 0)] <- 0.0
      
      return(list(Prob = prob, AccList = AccList))
    }
    return(list(Prob = probability, AccList = AccList))
  } else {
    return(list(Prob = probability, AccList = AccList))
  }
}

#' Update Pheromone Levels
#'
#' Updates the pheromone levels based on model performance. This is one of the
#' core functions of the hACO algorithm, as it reinforces the pathways (variable-factor
#' connections) that lead to good-fitting models.
#' 
#' The key improvement in this function is that it only adds pheromone to
#' variables with Ind=1, meaning only the active factor loadings get reinforced.
#' This ensures the algorithm properly converges to optimal factor structures.
#' High quality models receive stronger reinforcement to accelerate convergence.
#' 
#' @param Model Current model structure
#' @param Phe_level Current pheromone levels (or NULL for initialization)
#' @param Pheromone Pheromone value for the current model
#' @param n_factors Number of factors in the model
#' @param n_variables Number of observed variables
#' @param probability Current probability matrix
#' @param decay_rate Rate at which existing pheromone decays (default: 0.99)
#' @return List with updated pheromone levels
PheLevel <- function(Model, Phe_level = NULL, Pheromone, n_factors, n_variables, 
                     probability, decay_rate = 0.99) {
  # Initialize pheromone levels if not provided
  if (missing(Phe_level) | is.null(Phe_level)) {
    Phe_level <- rep(0.0001, n_factors*n_variables)
  }
  
  # Apply pheromone decay (evaporation) to all connections
  # This helps prevent premature convergence and encourages exploration
  Phe_level <- Phe_level * decay_rate
  
  Ind <- Model$ModelInd
  
  # Only add pheromone to variables with Ind=1
  # This is the key improvement to properly update only active connections
  # Add explicit debug output to verify pheromone updates
  if (Pheromone > 0) {
    # Use an adaptive increment based on model quality
    if (Pheromone > 0.8) {
      # For high quality models, use stronger reinforcement for faster convergence
      PheIncrement <- rep(0, length(Ind))
      PheIncrement[which(Ind == 1)] <- Pheromone^2
      Phe_level <- Phe_level + PheIncrement
    } else {
      # For regular quality models, use standard reinforcement
      PheIncrement <- rep(0, length(Ind))
      PheIncrement[which(Ind == 1)] <- Pheromone^3
      Phe_level <- Phe_level + PheIncrement
    }
  }
  
  # Maintain constraints from probability matrix
  # This ensures that prohibited loadings (probability=0) don't accumulate pheromone
  Phe_level[which(t(probability) == 0)] <- 0.0001
  
  # Ensure minimum pheromone level to allow exploration
  min_pheromone <- 0.0001
  Phe_level[Phe_level < min_pheromone] <- min_pheromone
  
  # Apply maximum pheromone cap to prevent domination by single paths
  max_pheromone <- 0.9999
  if (any(Phe_level > max_pheromone)) {
    Phe_level[Phe_level > max_pheromone] <- max_pheromone
  }
  
  return(list(PheLevel = Phe_level))
}

#' Model Run Function
#'
#' Runs a single model evaluation and updates all parameters
#' 
#' @param Model Current model
#' @param Data Dataset containing variables
#' @param dicrimination Discrimination parameters for fit indices
#' @param standards Standard values for fit indices
#' @param scopes Scope parameters for fit indices
#' @param P_dicrimination Discrimination parameter for probability
#' @param P_standards Standard value for probability
#' @param P_scopes Scope parameter for probability
#' @param AccList Current acceptance list
#' @param BanList Current ban list
#' @param Phe_level Current pheromone levels
#' @param probability Current probability matrix
#' @param n_factors Number of factors
#' @param n_variables Number of variables
#' @param numRuns Counter for model runs
#' @param LowCut Lower threshold for accepting a model
#' @param HighCut Upper threshold for updating probabilities
#' @param Punish_rate Rate for punishing model complexity
#' @param Fit Vector of fit indices to calculate
#' @param group Group variable name (if multi-group analysis)
#' @param decay_rate Rate at which existing pheromone decays
#' @return List with updated values for all parameters
ModelRun <- function(Model, Data, dicrimination, standards, 
                     scopes = c(0.015, 0.015, 0.015), 
                     P_dicrimination, P_standards, P_scopes,
                     AccList, BanList, Phe_level, probability, n_factors, n_variables, 
                     numRuns = 0, LowCut, HighCut, Punish_rate, Fit, group = NULL, 
                     decay_rate = 0.99){
  
  # Calculate pheromone
  pheromone <- Pheromone(
    Model = Model, 
    Data = Data, 
    dicrimination = dicrimination, 
    standards = standards, 
    scopes = scopes, 
    BanList = BanList, 
    numRuns = numRuns, 
    LowCut = LowCut, 
    n_factors = n_factors, 
    n_variables = n_variables, 
    Punish_rate = Punish_rate, 
    Fit = Fit, 
    group = group
  )
  
  # Update probabilities
  prob <- ProbCompensation(
    Model = Model, 
    Pheromone = pheromone$Pheromone, 
    AccList = AccList, 
    probability = probability, 
    n_factors = n_factors, 
    n_variables = n_variables,
    P_dicrimination = P_dicrimination, 
    P_standards = P_standards, 
    P_scopes = P_scopes, 
    HighCut = HighCut
  )
  
  # Update pheromone levels
  phe_level <- PheLevel(
    Model = Model, 
    Phe_level = Phe_level, 
    Pheromone = pheromone$Pheromone, 
    n_factors = n_factors, 
    n_variables = n_variables, 
    probability = probability,
    decay_rate = decay_rate
  )
  
  # Calculate entropy
  entropy <- Entropy(Phe_level = phe_level$PheLevel)
  
  return(list(
    PheLevel = phe_level$PheLevel, 
    Prob = prob$Prob, 
    Entropy = entropy$Entropy, 
    NumRuns = pheromone$numRuns, 
    BanList = pheromone$BanList, 
    AccList = prob$AccList, 
    Pheromone = pheromone$Pheromone, 
    BIC = pheromone$BIC, 
    ModelInd = paste(Model$ModelInd, collapse = "")
  ))
}

#==============================================================================
#                           MAIN hACO FUNCTION
#==============================================================================

#' Heuristic Ant Colony Optimization for CFA Model Search
#'
#' Main function for running the hACO algorithm to find optimal factor structures
#' in multi-group Confirmatory Factor Analysis (CFA). The algorithm automates the 
#' process of discovering the best factor structure based on model fit indices.
#'
#' The algorithm proceeds through three phases:
#' 1. Initial Search: Find a viable starting model with acceptable fit
#' 2. Neighborhood search to refine and optimize the model
#' 3. Final model selection: Analyze pheromone distribution to select the optimal model
#'
#' The key innovation in this approach is the use of pheromone distributions to guide
#' the search toward optimal factor structures, similar to how ants find efficient paths.
#' Models with good fit add more pheromone to their factor-variable connections,
#' gradually steering the algorithm toward the best solution.
#' 
#' @param Data Dataset containing variables for analysis
#' @param n_factors Number of factors to extract
#' @param n_variables Number of observed variables
#' @param BanList Initial list of banned configurations (NULL for new search)
#' @param AccList Initial list of accepted configurations (NULL for new search)
#' @param Phe_level Initial pheromone levels (NULL for new search)
#' @param loaded List specifying which variables should load on each factor (prior knowledge)
#'        Format: Either empty list() for no constraints, or list with exactly n_factors elements
#'        Each element contains variable indices or empty c() if no constraints for that factor
#'        Example: list(c(1,2,3), c(4,5,6)) forces variables 1-3 to load on factor 1, 4-6 on factor 2
#' @param unloaded List specifying which variables should NOT load on each factor
#'        Format: Either empty list() for no constraints, or list with exactly n_factors elements
#'        Each element contains variable indices or empty c() if no constraints for that factor
#'        Example: list(c(4,5,6), c(1,2,3)) prohibits variables 4-6 from loading on factor 1
#' @param dicrimination Vector of discrimination parameters for fit indices
#'        Controls how strongly the algorithm responds to differences in fit measures
#' @param standards Vector of standard/target values for fit indices
#'        Threshold values that determine when a model is considered good
#' @param scopes Vector of scope parameters for sensitivity of fit indices
#'        Controls the range of influence around standard values
#' @param P_dicrimination Discrimination parameter for probability updates
#' @param P_standards Standard value for probability updates
#' @param P_scopes Scope parameter for probability sensitivity
#' @param maxRun Maximum number of iterations for the neighborhood search
#' @param LowCut Lower threshold for accepting a model (pheromone cutoff)
#' @param HighCut Upper threshold for updating probabilities (pheromone threshold)
#' @param ChangeRun Maximum attempts before resetting ban list when stuck
#' @param Punish_rate Rate for penalizing model complexity 
#'        Higher values favor simpler models
#' @param Ban_neibor_length Length of neighbor ban list for local search
#' @param Fit Vector of fit indices to calculate:
#'        1=SRMR, 2=RMSEA, 3=IFI, 4=TLI, 5=CFI, 6=MFI
#'        Default c(5,2,3) uses CFI, RMSEA, and IFI
#' @param group Group variable name for multi-group analysis (optional)
#'        If provided, performs multi-group CFA using this variable
#'
#' @return List containing:
#'   \item{OptModel}{Readable optimal model structure}
#'   \item{OptModelResults}{Fit results for optimal model}
#'   \item{OptModelRun}{Lavaan syntax for the optimal model}
#'   \item{PheLevel}{Final pheromone levels showing factor-variable relationships}
#'   \item{Prob}{Final probability matrix}
#'   \item{BanList}{Final list of banned configurations}
#'   \item{AccList}{Final list of accepted configurations}
#'   \item{BestSearchModel}{Best model found during search phase}
#'   \item{BestEst}{Fit results for best search model}
#'
#' @examples
#' # Simple example with 2 factors and 6 variables (no constraints):
#' result <- hACO(
#'   Data = your_data,
#'   n_factors = 2,
#'   n_variables = 6,
#'   loaded = list(),      # No loading constraints
#'   unloaded = list(),    # No prohibition constraints
#'   maxRun = 100
#' )
#'
#' # Multi-group example with prior knowledge:
#' result <- hACO(
#'   Data = your_data,
#'   n_factors = 3,
#'   n_variables = 9,
#'   loaded = list(c(1,2,3), c(4,5,6), c(7,8,9)),  # Force variable groupings
#'   unloaded = list(c(4,5,6,7,8,9), c(1,2,3,7,8,9), c(1,2,3,4,5,6)),
#'   group = "group_variable",
#'   maxRun = 200
#' )
hACO <- function(Data,
                n_factors,
                n_variables,
                BanList = NULL,
                AccList = NULL,
                Phe_level = NULL,
                loaded,
                unloaded,
                dicrimination = c(1, 1, 1),
                standards = c(0.93, 1-0.07, 0.93),
                scopes = c(0.015, 0.015, 0.015),
                P_dicrimination = 0.5,
                P_standards = 0.95,
                P_scopes = 0.005,
                maxRun = 15000,
                LowCut = 0.001,
                HighCut = 0.75,
                ChangeRun = 1000,
                Punish_rate = 0.1,
                Ban_neibor_length = 5,
                Fit = c(5, 2, 3),  # Default: CFI, RMSEA, SRMR
                group = NULL) {
  
  # Calculate the initial probability matrix incorporating prior knowledge
  probability <- ProbBuild(
    n_factors = n_factors,
    n_variables = n_variables,
    loaded = loaded,
    unloaded = unloaded
  )
  
  # Initialize tracking variables
  entropy_pre   <- 100
  diff          <- 100
  numRuns       <- 1
  numAcc        <- 0
  numTotal      <- 0
  numRuns_pre   <- 0
  diffRuns      <- 0
  MaxPheromone  <- -Inf
  Ban_neibor    <- NULL
  
  # Ensure data is a data frame
  Data <- as.data.frame(Data)
  
  # Generate variable names if not present
  if (is.null(names(Data))) {
    x <- rep("X", n_variables)
    x <- paste(x, c(1:n_variables), sep = "")
    names(Data) <- x
  }
  
  # Check for invalid thresholds
  if (HighCut < LowCut) {
    stop("HighCut cannot be lower than LowCut", call. = FALSE)
  }
  EntropyStart <- FALSE
  
  # Build initial model
  model <- ModelBuilder(
    Data = Data,
    n_factors = n_factors,
    n_variables = n_variables,
    probability = probability,
    BanList = BanList,
    ChangeRun = ChangeRun
  )
  
  # If all models are banned, reset ban list and lower threshold
  if (model$Change == TRUE) {
    BanList <- NULL
    if (HighCut > LowCut) {
      LowCut <- LowCut
      HighCut <- HighCut - 0.10
    } else {
      return(print("Check number of factors, loaded, or unloaded settings"))
    }
  }
  
  # Run initial model
  res <- ModelRun(
    Model = model,
    Data = Data,
    dicrimination = dicrimination,
    standards = standards,
    scopes = scopes,
    P_dicrimination = P_dicrimination,
    P_standards = P_standards,
    P_scopes = P_scopes,
    AccList = AccList,
    BanList = BanList,
    Phe_level = Phe_level,
    probability = probability,
    n_factors = n_factors,
    n_variables = n_variables,
    numRuns = numRuns,
    LowCut = LowCut,
    HighCut = HighCut,
    Punish_rate = Punish_rate,
    Fit = Fit,
    group = group,
    decay_rate = 0.99
  )
  
  cat("\n========================================================================\n")
  cat("  Starting hACO Multi-Group CFA Algorithm\n")
  cat("  Factors:", n_factors, "| Variables:", n_variables, "| Max iterations:", maxRun, "\n")
  if (!is.null(group)) cat("  Multi-group analysis with group variable:", group, "\n")
  cat("========================================================================\n")
  
  # Display motivational phrase before Phase 1
  cat("\n  \"Preparation, Patience from Old Drives!\"\n\n")
  
  # Show Phase 1 start message
  cat("PHASE 1: Searching for initial viable model...\n")
  
  Start <- model$ModelInd
  BestSearch <- res$ModelInd
  
  patience_limit <- 0
  patience_Nolimit <- 10 * n_factors*n_variables
  
  # PHASE 1: INITIAL SEARCH
  # Keep searching until a model with non-zero pheromone is found
  # This establishes a baseline good model before the main search
  phase1_complete <- FALSE
  
  while (res$Pheromone == 0 && !phase1_complete) {
    model <- ModelBuilder(
      Data = Data,
      n_factors = n_factors,
      n_variables = n_variables,
      probability = probability,
      BanList = BanList,
      ChangeRun = ChangeRun
    )
    
    if (model$Change == TRUE) {
      BanList <- NULL
      if (HighCut > LowCut) {
        LowCut <- LowCut
        HighCut <- HighCut - 0.10
      } else {
        return(print("Check n_factor, loaded, or unloaded settings"))
      }
    }
    
    res <- ModelRun(
      Model = model,
      Data = Data,
      dicrimination = dicrimination,
      standards = standards,
      scopes = scopes,
      P_dicrimination = P_dicrimination,
      P_standards = P_standards,
      P_scopes = P_scopes,
      AccList = AccList,
      BanList = BanList,
      Phe_level = Phe_level,
      probability = probability,
      n_factors = n_factors,
      n_variables = n_variables,
      numRuns = numRuns,
      LowCut = LowCut,
      HighCut = HighCut,
      Punish_rate = Punish_rate,
      Fit = Fit,
      group = group,
      decay_rate = 0.99
    )
    
    if (res$Pheromone > MaxPheromone) {
      BestSearch <- res$ModelInd
      MaxPheromone <- res$Pheromone
    }
    
    # Update parameters
    Phe_level <- res$PheLevel
    BanList <- res$BanList
    AccList <- res$AccList
    probability <- res$Prob
    Start <- model$ModelInd
    BestSearch <- res$ModelInd
    
    # Check patience limit
    patience_limit <- patience_limit + 1
    if (patience_limit >= patience_Nolimit || res$Pheromone > 0) {
      phase1_complete <- TRUE
    }
    
    # Display progress without accumulating output
    DisplayProgress(patience_limit, patience_Nolimit, numRuns, 
                   length(AccList), MaxPheromone, FALSE)
  }
  
  # Clear the line after Phase 1
  utils::flush.console()
  cat("\r", strrep(" ", 120), "\r", sep="")
  
  # Phase 1 completion message
  if (res$Pheromone > 0) {
    DisplayWarning(paste("PHASE 1 complete: Found viable model with pheromone =", 
        sprintf("%.4f", res$Pheromone)))
  } else {
    DisplayWarning("PHASE 1 complete: Reached patience limit without finding viable model")
  }
  
  # Clearly indicate transition to Phase 2 using DisplayWarning function
  DisplayWarning("PHASE 2: Starting neighborhood search optimization...")
  
  # Clear the console buffer before starting Phase 2 progress display
  # No additional clearing needed as DisplayWarning already clears the screen
  
  neiber_res <- res
  Good_neibor <- NULL
  
  # PHASE 2: NEIGHBORHOOD SEARCH
  # Main loop for neighborhood search - systematically explore models that
  # are one-change away from current good models to refine the solution
  decay_rate <- 0.99  # Initial pheromone decay rate
  stagnation_counter <- 0  # Counter for detecting search stagnation
  models_evaluated_last_iter <- 0  # Track models evaluated in last iteration
  consecutive_iterations_with_no_models <- 0 # Count iterations with no models evaluated
  
  # Reset numTotal to ensure we display all iterations
  numTotal <- 0
  
  # Variable to track current status message
  current_status <- NULL
  
  # Display Phase 2 start information
  cat("\n===== PHASE 2: Starting Neighborhood Search =====\n")
  utils::flush.console()
  
  # Display initial progress information
  cat(sprintf("%3d%%", 0), 
      " [", paste(rep(" ", 40), collapse=""),
      "] ", 0, "/", maxRun, 
      " | Models: ", numRuns, " | Accepted: ", length(AccList), 
      " | Max Ph: ", sprintf("%.4f", MaxPheromone), sep="")
  utils::flush.console()
  
  while ((numRuns < 3000) &
         (numAcc < 200) &
         ((abs(diff) > 0.0001) | (diffRuns < 1000)) & 
         (numTotal < maxRun)) {
    
    # Ensure progress is printed each iteration - using simple clearing
    cat("\r")
    utils::flush.console()
    
    Good_ph <- NULL
    Good_ind <- NULL
    Good_bic <- NULL
    models_evaluated_this_iter <- 0  # Counter for models evaluated in this iteration
    skipped_positions <- 0
    invalid_models <- 0
    banned_models <- 0
    
    # Reset current status
    current_status <- NULL
    
    # Every 50 iterations, update status message
    if (numTotal %% 50 == 0 && numTotal > 0) {
      current_status <- paste("Iteration milestone:", numTotal)
    }
    
    # Every 100 iterations, adjust decay rate to balance exploration and exploitation
    if (numTotal > 0 && numTotal %% 100 == 0) {
      if (length(AccList) < 5) {
        # If few models are accepted, increase exploration
        decay_rate <- max(0.95, decay_rate - 0.01) 
        current_status <- paste("Few accepted models, increased exploration (decay_rate =", decay_rate, ")")
      } else if (numTotal > 500 && abs(diff) < 0.001 && diffRuns > 500) {
        # If search is stagnating, increase exploration
        decay_rate <- max(0.95, decay_rate - 0.02)
        current_status <- paste("Search stagnating, increased exploration (decay_rate =", decay_rate, ")")
      }
      
      # 10% chance of random restart to increase exploration
      if (runif(1) < 0.1) {
        current_status <- "Performing random restart to increase exploration"
        model <- ModelBuilder(
          Data = Data,
          n_factors = n_factors,
          n_variables = n_variables,
          probability = probability,
          BanList = BanList,
          ChangeRun = ChangeRun
        )
        Start <- model$ModelInd
      }
    }
    
    # Check for search stagnation - if previous iteration found no models
    if (models_evaluated_last_iter == 0) {
      stagnation_counter <- stagnation_counter + 1
      consecutive_iterations_with_no_models <- consecutive_iterations_with_no_models + 1
      
      # After 3 iterations with no models evaluated, force random restart
      if (stagnation_counter >= 3 || consecutive_iterations_with_no_models >= 5) {
        current_status <- "Search stagnating, forcing random restart"
        # Clear neighbor ban list to allow revisiting previous neighborhoods
        Ban_neibor <- NULL
        
        # Generate completely random new starting point while maintaining constraints
        random_ind <- rbinom(n_variables*n_factors, 1, as.vector(t(probability)))
        ind_use <- matrix(random_ind, n_factors, n_variables, byrow = TRUE)
        
        # Ensure each variable loads on at least one factor
        while (is.element(0, colSums(ind_use))) {
          repair <- which(colSums(ind_use) == 0)  # Find variables with no loadings
          repair_line <- which(probability[n_factors,] != 0)  # Find eligible positions
          repair_element <- repair_line[which(repair_line %in% repair)]  # Select variables to repair
          if (length(repair_element) > 0) {
            ind_use[n_factors, repair_element] <- 1  # Add loading to last factor
          } else {
            # If no eligible positions, try other factors
            for (f in (n_factors-1):1) {
              repair_line <- which(probability[f,] != 0)
              repair_element <- repair_line[which(repair_line %in% repair)]
              if (length(repair_element) > 0) {
                ind_use[f, repair_element[1]] <- 1
                break
              }
            }
          }
        }
        
        Start <- c(t(ind_use))
        stagnation_counter <- 0
        
        # More aggressive ban list reduction when repeatedly stagnating
        if (consecutive_iterations_with_no_models >= 10 && length(BanList) > 10) {
          keep_count <- min(10, max(5, ceiling(length(BanList) * 0.1)))
          BanList <- sample(BanList, keep_count)
          current_status <- paste("Severely reduced ban list size to", length(BanList), "models")
          consecutive_iterations_with_no_models <- 0
        } else if (length(BanList) > 100) {
          keep_count <- min(75, max(50, ceiling(length(BanList) * 0.75)))
          BanList <- sample(BanList, keep_count)
          current_status <- paste("Reduced ban list size to", length(BanList), "models")
        }
      }
    } else {
      stagnation_counter <- 0
      consecutive_iterations_with_no_models <- 0
    }
    
    # Loop through all positions to find neighbors
    for (i in 1:length(Start)) {
      iter <- Start
      
      # Skip positions with fixed probabilities
      if (is.element(i, which(t(probability) == 1 | t(probability) == 0))) {
        skipped_positions <- skipped_positions + 1
        next
      } else {
        # Generate neighbor by flipping one position
        iter[i] <- ifelse(Start[i] == 1, 0, 1)
        modelind <- paste(iter, collapse = "")
        
        ind <- matrix(iter, n_factors, n_variables, byrow = TRUE)
        
        # Check if model is valid (each variable loads on at least one factor)
        if (is.element(0, colSums(ind))) {
          invalid_models <- invalid_models + 1
          next
        }
        
        # Check if model is in ban list
        if (is.element(paste(iter, collapse = ""), BanList)) {
          banned_models <- banned_models + 1
          next
        }
        
        # Check if model is in neighbor ban list
        if (is.element(modelind, Ban_neibor)) {
          banned_models <- banned_models + 1
          next
        }
          
        # Rebuild model from indicators
        model <- ModelReBuild(
          ModInd = paste(iter, collapse = ""),
          n_factors = n_factors,
          n_variables = n_variables,
          Data = Data,
          group = group
        )
        
        # Run the model
        neiber_res <- ModelRun(
          Model = model,
          Data = Data,
          dicrimination = dicrimination,
          standards = standards,
          scopes = scopes,
          P_dicrimination = P_dicrimination,
          P_standards = P_standards,
          P_scopes = P_scopes,
          AccList = AccList,
          BanList = BanList,
          Phe_level = Phe_level,
          probability = probability,
          n_factors = n_factors,
          n_variables = n_variables,
          numRuns = numRuns,
          LowCut = LowCut,
          HighCut = HighCut,
          Punish_rate = Punish_rate,
          Fit = Fit,
          group = group,
          decay_rate = decay_rate
        )
        
        # Count models evaluated in this iteration
        models_evaluated_this_iter <- models_evaluated_this_iter + 1
        
        # Track best model found
        if (neiber_res$Pheromone >= MaxPheromone) {
          BestSearch <- neiber_res$ModelInd
          MaxPheromone <- neiber_res$Pheromone
        }
        
        # Store results
        Good_ph <- rbind(Good_ph, neiber_res$Pheromone)
        Good_ind <- rbind(Good_ind, neiber_res$ModelInd)
        Good_bic <- rbind(Good_bic, neiber_res$BIC)
        
        # Update parameters
        Phe_level <- neiber_res$PheLevel
        BanList <- neiber_res$BanList
        AccList <- neiber_res$AccList
        probability <- neiber_res$Prob
        entropy <- neiber_res$Entropy
        diff <- entropy_pre - neiber_res$Entropy
        numRuns <- neiber_res$NumRuns  # Update model count immediately
        
        # Check stop criteria based on entropy
        if (((diff != 0) & (numTotal > 0))) {
          EntropyStart <- TRUE
        }
        if ((abs(diff) < 0.0001) & (EntropyStart == TRUE)) {
          diffRuns <- diffRuns + (neiber_res$NumRuns - numRuns_pre)
        }
        entropy_pre <- neiber_res$Entropy
      }      
    }
    
    # If no models were evaluated in this iteration, display debug information as status
    if (models_evaluated_this_iter == 0) {
      current_status <- paste("No models evaluated | Skipped:", skipped_positions, 
                           "| Invalid:", invalid_models, 
                           "| Banned:", banned_models)
      
      # Create a more relaxed probability matrix to encourage exploration
      relaxed_prob <- matrix(0.5, n_factors, n_variables)
      # Preserve only strong constraints (values of 0 or 1)
      relaxed_prob[probability == 0] <- 0
      relaxed_prob[probability == 1] <- 1
      
      # Aggressive ban list reduction if we've had multiple stagnations
      if (consecutive_iterations_with_no_models >= 3 && length(BanList) > 20) {
        # Keep only a small random sample of banned models to allow exploration
        BanList <- sample(BanList, min(20, length(BanList)))
        current_status <- paste(current_status, "| Reduced ban list to", length(BanList), "models")
      }
      
      # Build a new model with reduced constraints
      model <- ModelBuilder(
        Data = Data,
        n_factors = n_factors,
        n_variables = n_variables,
        probability = relaxed_prob,
        BanList = BanList,
        ChangeRun = ChangeRun
      )
      Start <- model$ModelInd
      Ban_neibor <- NULL  # Clear neighbor ban list
    }
    
    # Choose next starting point based on neighbors
    if (!is.null(Good_ph) && length(Good_ph) > 0) {
      Good_neibor <- data.frame(Good_ph, Good_ind, Good_bic, stringsAsFactors = FALSE)
      
      if (max(Good_neibor$Good_ph) == 0) {
        # If all neighbors have zero pheromone, choose based on BIC
        neiber_max <- which(Good_neibor$Good_bic == min(Good_neibor$Good_bic))
        if (length(neiber_max) > 1) {
          neiber_choose <- sample(neiber_max, 1)
        } else {
          neiber_choose <- neiber_max
        }
      } else {
        # Choose neighbor with highest pheromone
        neiber_max <- which(Good_neibor$Good_ph == max(Good_neibor$Good_ph))
        if (length(neiber_max) > 1) {
          neiber_choose <- sample(neiber_max, 1)
        } else {
          neiber_choose <- neiber_max
        }
      }
      
      neiber_start <- neiber_choose
      
      # Taboo search implementation: Add best model to Ban_neibor list
      # Maintain FIFO (First-In-First-Out) queue with maximum length of Ban_neibor_length
      best_model_ind <- as.character(Good_neibor$Good_ind[neiber_start])
      
      # Only add if not already in ban list
      if (!is.element(best_model_ind, Ban_neibor)) {
        # Remove oldest entry if list is at capacity
        if (length(Ban_neibor) >= Ban_neibor_length) {
          Ban_neibor <- Ban_neibor[-1]  # Remove first (oldest) element
        }
        # Add new best model to end of list
        Ban_neibor <- c(Ban_neibor, best_model_ind)
      }
      
      # Set new starting point
      Start <- StartBuild(Good_neibor$Good_ind[neiber_start], 
                          n_factors = n_factors, 
                          n_variables = n_variables)
    } else if (MaxPheromone > -Inf) {
      # No valid neighbors found, use best model so far
      modelbest <- ModelReBuild(BestSearch, n_factors, n_variables, Data = Data, group = group)
      Start <- modelbest$ModelInd
    }
    
    # Update counters
    numRuns_pre <- numRuns
    numAcc <- length(AccList)
    numTotal <- numTotal + 1
    models_evaluated_last_iter <- models_evaluated_this_iter
    
    # Early stopping condition with stability check
    if (numTotal > 100 && MaxPheromone > 0.85) {
      # Track stability of best solution
      if (!exists("best_ph_unchanged_count")) {
        best_ph_unchanged_count <- 1
        last_best_ph <- MaxPheromone
      } else if (abs(MaxPheromone - last_best_ph) < 0.0001) {
        best_ph_unchanged_count <- best_ph_unchanged_count + 1
        
        # Stop if best solution is stable for some iterations
        if (best_ph_unchanged_count >= 20 && numTotal > 150) {
          current_status <- paste("Best model stable for", best_ph_unchanged_count, 
              "iterations with Ph =", sprintf("%.4f", MaxPheromone))
          # Final progress display before stopping
          DisplayProgress(numTotal, maxRun, numRuns, numAcc, MaxPheromone, 
                      newline = TRUE, status = current_status)
          break
        }
      } else {
        best_ph_unchanged_count <- 1
        last_best_ph <- MaxPheromone
      }
      
      # Only display high quality model message periodically to avoid cluttering
      if (numTotal %% 10 == 0) {
        current_status <- paste("High quality model found. Will stop if stable")
      }
    } else if (numTotal > 200) {
      # Fallback stopping condition based on iteration count
      current_status <- paste("Stopping search after", numTotal, "iterations")
      DisplayProgress(numTotal, maxRun, numRuns, numAcc, MaxPheromone, 
                   newline = TRUE, status = current_status)
      break
    }
    
    if (numTotal == maxRun) {
      current_status <- "Maximum iterations reached!"
    }
    
    # Use direct cat statements to display progress, ensuring it's always shown
    cat("\r")  # Ensure cursor returns to start of line
    cat(sprintf("%3d%%", as.integer(100 * numTotal / maxRun)), 
        " [", paste(rep("=", as.integer(40 * numTotal / maxRun)), collapse=""),
        paste(rep(" ", 40 - as.integer(40 * numTotal / maxRun)), collapse=""),
        "] ", numTotal, "/", maxRun, 
        " | Models: ", numRuns, " | Accepted: ", numAcc, 
        " | Max Ph: ", sprintf("%.4f", MaxPheromone), sep="")
    
    # Add status information (if available)
    if (!is.null(current_status)) {
      cat("  ", current_status, sep="")
    }
    
    # Force console refresh
    utils::flush.console()
    
    # Add a newline every 10 iterations to keep display flowing
    if (numTotal %% 10 == 0) {
      cat("\n")
      utils::flush.console()
    }
    
    # Update counters for next iteration
    numRuns_pre <- numRuns
    numAcc <- length(AccList)
    numTotal <- numTotal + 1
    models_evaluated_last_iter <- models_evaluated_this_iter
  }
  
  # PHASE 3: FINAL MODEL SELECTION
  # Find optimal model from pheromone distribution by testing
  # different thresholds to select variables
  rate_total <- NULL
  model_sum <- NULL
  opt_fit1 <- NULL
  opt_fit2 <- NULL
  opt_fit3 <- NULL
  
  # Simplified cleaning at Phase 3 start
  utils::flush.console()
  cat("\r", strrep(" ", 120), "\r", sep="")
  
  DisplayWarning("PHASE 3: Selecting optimal model from pheromone distribution...")
  
  # Try different thresholds for selecting variables
  for (i in seq(0.98, 0.80, by = -0.01)) {
    Optmodel <- OptModel(
      Data = Data,
      PheLevel = Phe_level,
      n_factors = n_factors,
      n_variables = n_variables,
      rate = i,
      group = group,
      verbose = FALSE  # Don't show details during loop
    )
    OptEst <- OptModelEst(
      Model = Optmodel, 
      Data = Data, 
      Fit = Fit, 
      dicrimination = dicrimination, 
      standards = standards, 
      scopes = scopes, 
      n_factors = n_factors, 
      n_variables = n_variables, 
      Punish_rate = Punish_rate, 
      group = group
    )
    
    if (class(OptEst$Modelres)[1] == "try-error") {
      next
    } else {
      if (length(which(colSums(Optmodel$ModelPreTest) == 0)) == 0) {
        rate_total <- c(rate_total, i)
        model_sum <- c(model_sum, sum(Optmodel$ModelPreTest))
        opt_fit1 <- c(opt_fit1, OptEst$Fit1)
        opt_fit2 <- c(opt_fit2, OptEst$Fit2)
        opt_fit3 <- c(opt_fit3, OptEst$Fit3)
      }
    }
  }
  
  # Choose the optimal threshold based on fit improvement patterns
  if (is.null(rate_total)) {
    rate_opt <- 0.8
    cat("No valid models found, using default rate_opt =", rate_opt, "\n")
  } else {
    opt_sel <- as.data.frame(cbind(rate_total, model_sum, opt_fit1, opt_fit2, opt_fit3))
    opt_sel <- opt_sel[-which(duplicated(opt_sel[, 2])), ]
    rate_opt <- 0
    for (i in nrow(opt_sel):1) {
      if (i == 1) {
        rate_opt <- opt_sel[1, 1]
        cat("Using first model rate_opt =", rate_opt, "\n")
        break
      } else {
        # Look for significant improvements in fit indices
        if (((opt_sel[i, 3] - opt_sel[i - 1, 3]) < 0.02) |
            ((opt_sel[i, 4] - opt_sel[i - 1, 4]) < 0.02) |
            ((opt_sel[i, 5] - opt_sel[i - 1, 5]) < 0.02)) {
          next
        } else {
          rate_opt <- opt_sel[i, 1]
          cat("Selected model at index", i, "with rate_opt =", rate_opt, "\n")
          cat("Fit improvements:", 
              sprintf("%.3f", opt_sel[i, 3] - opt_sel[i - 1, 3]), 
              sprintf("%.3f", opt_sel[i, 4] - opt_sel[i - 1, 4]), 
              sprintf("%.3f", opt_sel[i, 5] - opt_sel[i - 1, 5]), "\n")
          break
        }
      }
    }
  }
  
  # Analyze pheromone distribution
  DisplayWarning("Analyzing pheromone distribution:")
  cat("Pheromone min:", min(Phe_level), "max:", max(Phe_level), "mean:", mean(Phe_level), "\n")
  cat("Number of variables with high pheromone (>", rate_opt * max(Phe_level), "):", 
      sum(Phe_level > rate_opt * max(Phe_level)), "out of", length(Phe_level), "\n\n")
  
  # Build final model using optimal threshold
  Optmodel <- OptModel(
    Data = Data,
    PheLevel = Phe_level,
    n_factors = n_factors,
    n_variables = n_variables,
    rate = rate_opt,
    group = group,
    verbose = TRUE  # Show details for final result
  )
  OptEst <- OptModelEst(
    Model = Optmodel, 
    Data = Data, 
    Fit = Fit, 
    dicrimination = dicrimination, 
    standards = standards, 
    scopes = scopes, 
    n_factors = n_factors, 
    n_variables = n_variables, 
    Punish_rate = Punish_rate, 
    group = group
  )
  
  # Also evaluate the best model found during search
  BestModel <- ModelReBuild(
    ModInd = BestSearch,
    n_factors = n_factors,
    n_variables = n_variables,
    Data = Data,
    group = group
  )
  BestEst <- OptModelEst(
    Model = BestModel, 
    Data = Data, 
    Fit = Fit, 
    dicrimination = dicrimination, 
    standards = standards, 
    scopes = scopes, 
    n_factors = n_factors, 
    n_variables = n_variables, 
    Punish_rate = Punish_rate, 
    group = group
  )
  
  # Display final pheromone distribution
  DisplayWarning("======== Final Pheromone Distribution Analysis ========")
  DisplayPheromoneDistribution(Phe_level, n_factors, n_variables, Data)
  
  DisplayWarning("Iteration Done")
  
  # Return results - both the optimal model based on pheromone distribution
  # and the best model found during the search process
  return(list(
    OptModel = Optmodel$OptModelShow,        # Readable optimal model
    OptModelResults = OptEst,                # Fit results for optimal model
    OptModelRun = Optmodel$CFAModel,         # Lavaan syntax for optimal model
    PheLevel = Phe_level,                    # Final pheromone levels
    Prob = probability,                      # Final probability matrix
    BanList = BanList,                       # Final ban list
    AccList = unique(AccList),               # Final acceptance list
    BestSearchModel = BestModel$CFAModel,    # Best model found during search
    BestEst = BestEst                        # Fit results for best model
  ))
} 