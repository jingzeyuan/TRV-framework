# Function for measurement invariance evaluation using ΔCFI
evaluate_invariance <- function(configural, metric, scalar, strict) {
  # Extract robust CFI values
  cfis <- c(
    Configural = fitMeasures(configural, "cfi.robust"),
    Metric = fitMeasures(metric, "cfi.robust"), 
    Scalar = fitMeasures(scalar, "cfi.robust"),
    Strict = fitMeasures(strict, "cfi.robust")
  )
  
  # Calculate ΔCFI
  dcfis <- c(
    `Metric-Configural` = cfis[2] - cfis[1],
    `Scalar-Metric` = cfis[3] - cfis[2], 
    `Strict-Scalar` = cfis[4] - cfis[3]
  )
  
  # Evaluate invariance (ΔCFI > -0.01 = supported)
  support <- dcfis > -0.01
  
  # Determine highest level achieved
  if (all(support)) {
    highest <- "Strict Invariance"
  } else if (support[1] && support[2]) {
    highest <- "Scalar Invariance"  
  } else if (support[1]) {
    highest <- "Metric Invariance"
  } else {
    highest <- "Configural Only"
  }
  
  # Create results table
  results <- data.frame(
    Model = c("configural", "metric", "scalar", "strict"),
    RobustCFI = round(cfis, 4),
    DCFI = c(NA, round(dcfis, 4)),
    Supported = c(NA, ifelse(support, "✓", "✗")),
    stringsAsFactors = FALSE
  )
  
  cat("\n--- Measurement Invariance Evaluation (ΔCFI > -0.01) ---\n")
  print(results, row.names = FALSE)
  cat("\nHighest Level Achieved:", highest, "\n")
  
  return(invisible(list(cfis = cfis, dcfis = dcfis, support = support, highest = highest)))
}
