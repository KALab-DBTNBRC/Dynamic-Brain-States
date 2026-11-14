
# ==============================================================================
# CLASSICAL TWIN MODEL ANALYSIS: ACE / AE / CE / E
# MZ and DZ twins | 23 phenotypes | OpenMx
# ==============================================================================
#
# Usage:
#   Rscript step15_ace_models.R <twin_workbook.xlsx> <output_dir>
#
# <twin_workbook.xlsx> needs two sheets per phenotype, "<phenotype>_MZ" and
# "<phenotype>_DZ", each with columns Twin1 and Twin2.
# Outputs in <output_dir>: ACE_Master_Results.rds, ACE_Master_Results.csv,
# ACE_Model_Diagnostics.csv


# ------------------------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c("OpenMx", "readxl", "knitr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(OpenMx)
library(readxl)
library(knitr)


# ------------------------------------------------------------------------------
# 2. SETUP FILE PATHS
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript step15_ace_models.R <twin_workbook.xlsx> <output_dir>")
}
file_path <- args[1]

if (!file.exists(file_path)) {
  stop("Excel file not found: ", file_path)
}

output_dir <- args[2]

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

output_rds_path <- file.path(
  output_dir,
  "ACE_Master_Results.rds"
)

output_csv_path <- file.path(
  output_dir,
  "ACE_Master_Results.csv"
)

diagnostics_path <- file.path(
  output_dir,
  "ACE_Model_Diagnostics.csv"
)

cat("\nInput workbook:", file_path, "\n")
cat("Output directory:", output_dir, "\n")

if (file.access(output_dir, 2) != 0) {
  stop("Output directory is not writable: ", output_dir)
}


# ------------------------------------------------------------------------------
# 3. DEFINE VARIABLES
# ------------------------------------------------------------------------------

var_names <- c(
  "FO_State 1", "FO_State 2", "FO_State 3", "FO_State 4",
  "DT_State 1", "DT_State 2", "DT_State 3", "DT_State 4",
  "AR_State 1", "AR_State 2", "AR_State 3", "AR_State 4",
  "SSAGA_Alc_D4_Ab_Dx",
  "SSAGA_Alc_D4_Dp_Dx",
  "SSAGA_Alc_12_Drinks_Per_Day",
  "SSAGA_Alc_12_Frq",
  "Total_Drinks_7days",
  "Num_Days_Drank_7days",
  "Risky_12Month",
  "Risky_Drinks_7days",
  "SSAGA_Alc_Age_1st_Use",
  "Odor_AgeAdj",
  "Odor_Unadj"
)

selVars <- c("Twin1", "Twin2")

available_sheets <- excel_sheets(file_path)

master_results <- list()
diagnostic_results <- list()
failed_variables <- character(0)

cat("\nNumber of variables:", length(var_names), "\n")
cat("Number of worksheets:", length(available_sheets), "\n")


# ------------------------------------------------------------------------------
# 4. CHECK WORKSHEET NAMES
# ------------------------------------------------------------------------------

missing_sheets <- character(0)

for (var in var_names) {
  
  mz_sheet <- paste0(var, "_MZ")
  dz_sheet <- paste0(var, "_DZ")
  
  if (!(mz_sheet %in% available_sheets)) {
    missing_sheets <- c(missing_sheets, mz_sheet)
  }
  
  if (!(dz_sheet %in% available_sheets)) {
    missing_sheets <- c(missing_sheets, dz_sheet)
  }
}

if (length(missing_sheets) > 0) {
  cat("\nMissing worksheets:\n")
  print(missing_sheets)
} else {
  cat("\nAll expected MZ and DZ worksheets found.\n")
}


# ------------------------------------------------------------------------------
# 5. FUNCTION: READ AND PREPARE TWIN DATA
# ------------------------------------------------------------------------------

read_twin_data <- function(sheet_name) {
  
  dat <- as.data.frame(
    read_excel(
      file_path,
      sheet = sheet_name
    )
  )
  
  if (!all(selVars %in% names(dat))) {
    stop(
      "Twin1/Twin2 columns missing in worksheet: ",
      sheet_name
    )
  }
  
  dat <- dat[, selVars, drop = FALSE]
  
  # Convert to numeric while detecting unexpected text
  for (j in seq_along(selVars)) {
    
    original <- dat[[j]]
    
    converted <- suppressWarnings(
      as.numeric(as.character(original))
    )
    
    invalid <- !is.na(original) & is.na(converted)
    
    if (any(invalid)) {
      stop(
        "Non-numeric values found in ",
        sheet_name,
        ", column ",
        selVars[j]
      )
    }
    
    dat[[j]] <- converted
  }
  
  # Remove rows where both twins are missing
  dat <- dat[
    !apply(is.na(dat), 1, all),
    ,
    drop = FALSE
  ]
  
  if (nrow(dat) < 10) {
    stop(
      "Insufficient twin pairs in worksheet: ",
      sheet_name
    )
  }
  
  return(dat)
}


# ------------------------------------------------------------------------------
# 6. FUNCTION: EXTRACT A, C, E PROPORTIONS
# ------------------------------------------------------------------------------

get_proportions <- function(model) {
  
  a_sq <- as.numeric(
    mxEval(shared.A, model)
  )[1]
  
  c_sq <- as.numeric(
    mxEval(shared.C, model)
  )[1]
  
  e_sq <- as.numeric(
    mxEval(shared.E, model)
  )[1]
  
  v_tot <- as.numeric(
    mxEval(shared.V, model)
  )[1]
  
  if (!is.finite(v_tot) || v_tot <= 0) {
    return(
      c(A = NA_real_, C = NA_real_, E = NA_real_)
    )
  }
  
  return(
    c(
      A = a_sq / v_tot,
      C = c_sq / v_tot,
      E = e_sq / v_tot
    )
  )
}


# ------------------------------------------------------------------------------
# 7. FUNCTION: EXTRACT MODEL FIT INFORMATION
# ------------------------------------------------------------------------------

get_model_info <- function(model) {
  
  n_parameters <- length(model$output$estimate)
  
  minus2LL <- suppressWarnings(as.numeric(model$output$fit))
  if (length(minus2LL) == 0 || !is.finite(minus2LL[1])) {
    stop("Model fit (-2LL) is missing or non-finite.")
  }
  minus2LL <- minus2LL[1]
  
  # Structural df for two observed variables in each of two groups:
  # 5 unique moments per group (2 means + 3 covariance elements) = 10.
  # The saturated model estimates all 10 moments; model df = 10 - k.
  model_df <- 10 - n_parameters
  
  # AIC = -2LL + 2k
  AIC_value <- minus2LL + 2 * n_parameters
  
  list(
    parameters = n_parameters,
    minus2LL = minus2LL,
    df = model_df,
    AIC = AIC_value
  )
}


# ------------------------------------------------------------------------------
# 8. FUNCTION: LIKELIHOOD RATIO TEST
# ------------------------------------------------------------------------------

compare_models <- function(full_model, restricted_model) {
  
  full_LL <- suppressWarnings(as.numeric(full_model$output$fit))
  restricted_LL <- suppressWarnings(as.numeric(restricted_model$output$fit))
  
  if (length(full_LL) == 0 || length(restricted_LL) == 0 ||
      !is.finite(full_LL[1]) || !is.finite(restricted_LL[1])) {
    return(c(diffLL = NA_real_, diffdf = NA_real_, p_value = NA_real_))
  }
  
  full_LL <- full_LL[1]
  restricted_LL <- restricted_LL[1]
  diff_LL <- restricted_LL - full_LL
  
  # The restricted model should have fewer free parameters.
  # This is equivalent to the df difference when both models use
  # the same observations and the same likelihood specification.
  full_parameters <- length(full_model$output$estimate)
  restricted_parameters <- length(restricted_model$output$estimate)
  diff_df <- full_parameters - restricted_parameters
  
  if (!is.finite(diff_LL) || diff_df <= 0) {
    return(c(diffLL = diff_LL, diffdf = diff_df, p_value = NA_real_))
  }
  
  if (diff_LL < -1e-4) {
    warning(
      "Restricted model has a lower -2LL than the full model; ",
      "check convergence and model nesting."
    )
    return(c(diffLL = diff_LL, diffdf = diff_df, p_value = NA_real_))
  }
  
  # Correct tiny numerical differences only.
  if (diff_LL < 0) diff_LL <- 0
  
  p_value <- pchisq(diff_LL, df = diff_df, lower.tail = FALSE)
  
  c(diffLL = diff_LL, diffdf = diff_df, p_value = p_value)
}


# ------------------------------------------------------------------------------
# 8A. FUNCTION: SAFELY EXTRACT OPTIMIZER STATUS
# ------------------------------------------------------------------------------

get_status <- function(model) {
  
  status_code <- tryCatch(
    model$output$status$code,
    error = function(e) NULL
  )
  
  if (is.null(status_code) || length(status_code) == 0 ||
      is.na(status_code[1])) {
    return(NA_integer_)
  }
  
  as.integer(status_code[1])
}


# ------------------------------------------------------------------------------
# 8B. FUNCTION: BUILD A TWO-GROUP SATURATED MODEL
# ------------------------------------------------------------------------------

build_saturated_model <- function(mzData, dzData) {
  
  make_group_model <- function(group_name, dat, mean_labels, cov_labels) {
    
    # Starting values from the observed data.
    mean_start <- vapply(dat, function(x) mean(x, na.rm = TRUE), numeric(1))
    mean_start[!is.finite(mean_start)] <- 0
    
    cov_start <- suppressWarnings(
      cov(dat, use = "pairwise.complete.obs")
    )
    
    # Ensure a finite, positive-definite starting covariance matrix.
    if (!is.matrix(cov_start) || any(!is.finite(cov_start))) {
      cov_start <- diag(1, 2)
    }
    cov_start <- (cov_start + t(cov_start)) / 2
    diag(cov_start) <- pmax(diag(cov_start), 1e-4)
    min_eigenvalue <- min(eigen(cov_start, symmetric = TRUE,
                                only.values = TRUE)$values)
    if (!is.finite(min_eigenvalue) || min_eigenvalue <= 1e-6) {
      adjustment <- if (is.finite(min_eigenvalue)) {
        abs(min_eigenvalue) + 1e-3
      } else {
        1e-3
      }
      cov_start <- cov_start + diag(adjustment, 2)
    }
    
    mxModel(
      group_name,
      mxData(observed = dat, type = "raw"),
      mxMatrix(
        type = "Full", nrow = 1, ncol = 2,
        free = TRUE, values = as.numeric(mean_start),
        labels = mean_labels, name = "expMean"
      ),
      mxMatrix(
        type = "Symm", nrow = 2, ncol = 2,
        free = TRUE, values = cov_start,
        labels = matrix(c(cov_labels[1], cov_labels[2],
                          cov_labels[2], cov_labels[3]), 2, 2),
        name = "expCov"
      ),
      mxExpectationNormal(
        covariance = "expCov", means = "expMean", dimnames = selVars
      ),
      mxFitFunctionML()
    )
  }
  
  sat_mz <- make_group_model(
    "SatMZ", mzData,
    c("mean_MZ_T1", "mean_MZ_T2"),
    c("var_MZ_T1", "cov_MZ", "var_MZ_T2")
  )
  
  sat_dz <- make_group_model(
    "SatDZ", dzData,
    c("mean_DZ_T1", "mean_DZ_T2"),
    c("var_DZ_T1", "cov_DZ", "var_DZ_T2")
  )
  
  mxModel(
    "SAT",
    sat_mz,
    sat_dz,
    mxFitFunctionMultigroup(c("SatMZ", "SatDZ"))
  )
}


# ------------------------------------------------------------------------------
# 8C. FUNCTION: CALCULATE MODEL CHI-SQUARE AGAINST SATURATED MODEL
# ------------------------------------------------------------------------------

get_chisq_vs_saturated <- function(model, saturated_model) {
  
  model_fit <- as.numeric(model$output$fit)[1]
  saturated_fit <- as.numeric(saturated_model$output$fit)[1]
  # Use the same structural df calculation as get_model_info().
  model_df <- 10 - length(model$output$estimate)
  
  if (!is.finite(model_fit) || !is.finite(saturated_fit) ||
      !is.finite(model_df) || model_df <= 0) {
    return(c(chi_square = NA_real_, df = model_df, p_value = NA_real_))
  }
  
  chi_square <- model_fit - saturated_fit
  
  # A negative difference beyond numerical tolerance indicates that the
  # saturated model or one of the fitted models needs to be checked.
  if (chi_square < -1e-4) {
    warning("Model has lower -2LL than saturated model; check model fits.")
    return(c(chi_square = chi_square, df = model_df, p_value = NA_real_))
  }
  
  chi_square <- max(0, chi_square)
  p_value <- pchisq(chi_square, df = model_df, lower.tail = FALSE)
  
  c(chi_square = chi_square, df = model_df, p_value = p_value)
}

# ------------------------------------------------------------------------------
# 9. MASTER LOOP
# ------------------------------------------------------------------------------

for (var in var_names) {
  
  cat("\n\n")
  cat("==================================================================\n")
  cat("PROCESSING VARIABLE:", var, "\n")
  cat("==================================================================\n")
  
  tryCatch({
    
    # --------------------------------------------------------------------------
    # 9A. IDENTIFY MZ AND DZ SHEETS
    # --------------------------------------------------------------------------
    
    mz_sheet <- paste0(var, "_MZ")
    dz_sheet <- paste0(var, "_DZ")
    
    if (!(mz_sheet %in% available_sheets) ||
        !(dz_sheet %in% available_sheets)) {
      
      stop(
        "Target sheets not found: ",
        mz_sheet, " or ", dz_sheet
      )
    }
    
    
    # --------------------------------------------------------------------------
    # 9B. READ DATA
    # --------------------------------------------------------------------------
    
    mzData <- read_twin_data(mz_sheet)
    dzData <- read_twin_data(dz_sheet)
    
    cat("MZ twin pairs:", nrow(mzData), "\n")
    cat("DZ twin pairs:", nrow(dzData), "\n")
    
    # Count observed individual values
    total_observations <-
      sum(!is.na(as.matrix(mzData))) +
      sum(!is.na(as.matrix(dzData)))
    
    
    # --------------------------------------------------------------------------
    # 9C. DEFINE SHARED ACE COMPONENTS
    # --------------------------------------------------------------------------
    
    shareVars <- mxModel(
      
      "shared",
      
      # Means
      mxMatrix(
        "Full", 1, 2,
        free = TRUE,
        values = c(0, 0),
        labels = c("mean", "mean"),
        name = "expMean"
      ),
      
      # Additive genetic path
      mxMatrix(
        "Full", 1, 1,
        free = TRUE,
        values = 0.6,
        labels = "a_param",
        name = "a"
      ),
      
      # Common environmental path
      mxMatrix(
        "Full", 1, 1,
        free = TRUE,
        values = 0.6,
        labels = "c_param",
        name = "c"
      ),
      
      # Unique environmental path
      mxMatrix(
        "Full", 1, 1,
        free = TRUE,
        values = 0.6,
        labels = "e_param",
        name = "e"
      ),
      
      # Variance components
      mxAlgebra(
        a %*% t(a),
        name = "A"
      ),
      
      mxAlgebra(
        c %*% t(c),
        name = "C"
      ),
      
      mxAlgebra(
        e %*% t(e),
        name = "E"
      ),
      
      # Total phenotypic variance
      mxAlgebra(
        A + C + E,
        name = "V"
      )
    )
    
    
    # --------------------------------------------------------------------------
    # 9D. DEFINE MZ MODEL
    # --------------------------------------------------------------------------
    
    modelMZ <- mxModel(
      
      "MZ",
      
      mxData(
        observed = mzData,
        type = "raw"
      ),
      
      mxAlgebra(
        
        rbind(
          
          cbind(
            shared.V,
            shared.A + shared.C
          ),
          
          cbind(
            shared.A + shared.C,
            shared.V
          )
          
        ),
        
        name = "expCovMZ"
      ),
      
      mxExpectationNormal(
        covariance = "expCovMZ",
        means = "shared.expMean",
        dimnames = selVars
      ),
      
      mxFitFunctionML()
    )
    
    
    # --------------------------------------------------------------------------
    # 9E. DEFINE DZ MODEL
    # --------------------------------------------------------------------------
    
    modelDZ <- mxModel(
      
      "DZ",
      
      mxData(
        observed = dzData,
        type = "raw"
      ),
      
      mxAlgebra(
        
        rbind(
          
          cbind(
            shared.V,
            0.5 * shared.A + shared.C
          ),
          
          cbind(
            0.5 * shared.A + shared.C,
            shared.V
          )
          
        ),
        
        name = "expCovDZ"
      ),
      
      mxExpectationNormal(
        covariance = "expCovDZ",
        means = "shared.expMean",
        dimnames = selVars
      ),
      
      mxFitFunctionML()
    )
    
    
    # --------------------------------------------------------------------------
    # 9F. FIT ACE MODEL
    # --------------------------------------------------------------------------
    
    cat("\nFitting ACE model...\n")
    
    aceModel <- mxModel(
      "ACE",
      shareVars,
      modelMZ,
      modelDZ,
      mxFitFunctionMultigroup(
        c("MZ", "DZ")
      )
    )
    
    aceFit <- mxRun(
      aceModel,
      silent = TRUE
    )
    
    
    # --------------------------------------------------------------------------
    # 9F2. FIT SATURATED MODEL FOR ABSOLUTE MODEL CHI-SQUARE
    # --------------------------------------------------------------------------
    
    cat("Fitting saturated model...\n")
    
    saturatedModel <- build_saturated_model(mzData, dzData)
    saturatedFit <- mxRun(saturatedModel, silent = TRUE)
    
    
    # --------------------------------------------------------------------------
    # 9G. FIT AE MODEL: DROP C
    # --------------------------------------------------------------------------
    
    cat("Fitting AE model...\n")
    
    aeModel <- omxSetParameters(
      aceFit,
      labels = "c_param",
      free = FALSE,
      values = 0
    )
    
    aeFit <- mxRun(
      mxModel(aeModel, name = "AE"),
      silent = TRUE
    )
    
    
    # --------------------------------------------------------------------------
    # 9H. FIT CE MODEL: DROP A
    # --------------------------------------------------------------------------
    
    cat("Fitting CE model...\n")
    
    ceModel <- omxSetParameters(
      aceFit,
      labels = "a_param",
      free = FALSE,
      values = 0
    )
    
    ceFit <- mxRun(
      mxModel(ceModel, name = "CE"),
      silent = TRUE
    )
    
    
    # --------------------------------------------------------------------------
    # 9I. FIT E MODEL: DROP A AND C
    # --------------------------------------------------------------------------
    
    cat("Fitting E model...\n")
    
    eModel <- omxSetParameters(
      aceFit,
      labels = c("a_param", "c_param"),
      free = FALSE,
      values = 0
    )
    
    eFit <- mxRun(
      mxModel(eModel, name = "E"),
      silent = TRUE
    )
    
    
    # --------------------------------------------------------------------------
    # 9J. CHECK OPTIMIZER STATUS
    # --------------------------------------------------------------------------
    
    fitted_models <- list(
      ACE = aceFit,
      AE  = aeFit,
      CE  = ceFit,
      E   = eFit
    )
    
    model_status <- vapply(
      fitted_models,
      get_status,
      integer(1)
    )
    
    print(model_status)
    
    if (any(!is.na(model_status) & model_status != 0)) {
      warning(
        "Non-zero optimizer status for ",
        var,
        ": ",
        paste(
          names(model_status)[
            !is.na(model_status) & model_status != 0
          ],
          collapse = ", "
        )
      )
    }
    
    
    # --------------------------------------------------------------------------
    # 9K. EXTRACT A, C, E PROPORTIONS
    # --------------------------------------------------------------------------
    
    prop_ACE <- get_proportions(aceFit)
    prop_AE <- get_proportions(aeFit)
    prop_CE <- get_proportions(ceFit)
    prop_E <- get_proportions(eFit)
    
    
    # --------------------------------------------------------------------------
    # 9L. EXTRACT MODEL FIT INDICES
    # --------------------------------------------------------------------------
    
    info_ACE <- get_model_info(aceFit)
    info_AE <- get_model_info(aeFit)
    info_CE <- get_model_info(ceFit)
    info_E <- get_model_info(eFit)
    
    # Absolute model fit: compare each model with the saturated model.
    chisq_ACE <- get_chisq_vs_saturated(aceFit, saturatedFit)
    chisq_AE  <- get_chisq_vs_saturated(aeFit, saturatedFit)
    chisq_CE  <- get_chisq_vs_saturated(ceFit, saturatedFit)
    chisq_E   <- get_chisq_vs_saturated(eFit, saturatedFit)
    
    
    # --------------------------------------------------------------------------
    # 9M. LIKELIHOOD RATIO TESTS
    # --------------------------------------------------------------------------
    
    comparison_AE <- compare_models(
      aceFit,
      aeFit
    )
    
    comparison_CE <- compare_models(
      aceFit,
      ceFit
    )
    
    comparison_E <- compare_models(
      aceFit,
      eFit
    )
    
    
    # --------------------------------------------------------------------------
    # 9N. CREATE SUMMARY TABLE (INCLUDING ACE HERITABILITY h2)
    # --------------------------------------------------------------------------
    
    summary_df <- data.frame(
      
      Variable = var,
      
      Model = c(
        "ACE (Base)",
        "AE (Drop C)",
        "CE (Drop A)",
        "E (Drop A & C)"
      ),
      
      Parameters = c(
        info_ACE$parameters,
        info_AE$parameters,
        info_CE$parameters,
        info_E$parameters
      ),
      
      minus2LL = c(
        info_ACE$minus2LL,
        info_AE$minus2LL,
        info_CE$minus2LL,
        info_E$minus2LL
      ),
      
      df = c(
        info_ACE$df,
        info_AE$df,
        info_CE$df,
        info_E$df
      ),
      
      saturated_minus2LL = rep(as.numeric(saturatedFit$output$fit)[1], 4),
      
      chi_square = c(
        chisq_ACE["chi_square"],
        chisq_AE["chi_square"],
        chisq_CE["chi_square"],
        chisq_E["chi_square"]
      ),
      
      chi_square_p = c(
        chisq_ACE["p_value"],
        chisq_AE["p_value"],
        chisq_CE["p_value"],
        chisq_E["p_value"]
      ),
      
      AIC = c(
        info_ACE$AIC,
        info_AE$AIC,
        info_CE$AIC,
        info_E$AIC
      ),
      
      diffLL = c(
        NA,
        comparison_AE["diffLL"],
        comparison_CE["diffLL"],
        comparison_E["diffLL"]
      ),
      
      diffdf = c(
        NA,
        comparison_AE["diffdf"],
        comparison_CE["diffdf"],
        comparison_E["diffdf"]
      ),
      
      p_value = c(
        NA,
        comparison_AE["p_value"],
        comparison_CE["p_value"],
        comparison_E["p_value"]
      ),
      
      prop_A = c(
        prop_ACE["A"],
        prop_AE["A"],
        prop_CE["A"],
        prop_E["A"]
      ),
      
      prop_C = c(
        prop_ACE["C"],
        prop_AE["C"],
        prop_CE["C"],
        prop_E["C"]
      ),
      
      prop_E = c(
        prop_ACE["E"],
        prop_AE["E"],
        prop_CE["E"],
        prop_E["E"]
      ),
      
      # Heritability from the full ACE model.
      # Reduced-model rows are left NA because their A estimates are
      # conditional on fixing one or more variance components to zero.
      h2 = c(
        prop_ACE["A"],
        NA_real_,
        NA_real_,
        NA_real_
      ),
      
      h2_percent = c(
        100 * prop_ACE["A"],
        NA_real_,
        NA_real_,
        NA_real_
      ),
      
      MZ_pairs = nrow(mzData),
      
      DZ_pairs = nrow(dzData),
      
      Observed_values = total_observations,
      
      stringsAsFactors = FALSE
    )
    
    
    # --------------------------------------------------------------------------
    # 9O. PRINT RESULTS
    # --------------------------------------------------------------------------
    
    cat("\nMODEL RESULTS:", var, "\n\n")
    
    print(
      kable(
        summary_df,
        caption = paste(
          "Twin Model Summary:",
          var
        ),
        digits = 4
      )
    )
    
    
    # --------------------------------------------------------------------------
    # 9P. STORE RESULTS
    # --------------------------------------------------------------------------
    
    master_results[[var]] <- summary_df
    
    diagnostic_results[[var]] <- data.frame(
      
      Variable = var,
      
      ACE_status = model_status["ACE"],
      
      AE_status = model_status["AE"],
      
      CE_status = model_status["CE"],
      
      E_status = model_status["E"],
      
      MZ_pairs = nrow(mzData),
      
      DZ_pairs = nrow(dzData),
      
      Saturated_minus2LL = as.numeric(saturatedFit$output$fit)[1],
      
      ACE_minus2LL = info_ACE$minus2LL,
      
      ACE_chi_square = chisq_ACE["chi_square"],
      
      ACE_chi_square_df = chisq_ACE["df"],
      
      ACE_chi_square_p = chisq_ACE["p_value"],
      
      AE_minus2LL = info_AE$minus2LL,
      
      CE_minus2LL = info_CE$minus2LL,
      
      E_minus2LL = info_E$minus2LL,
      
      stringsAsFactors = FALSE
    )
    
    cat("\nSuccessfully completed:", var, "\n")
    
    
  }, error = function(e) {
    
    cat(
      "\nERROR PROCESSING:",
      var,
      "\nError:",
      conditionMessage(e),
      "\n"
    )
    
    cat("Call stack at error:\n")
    calls <- sys.calls()
    print(tail(calls, 8))
    
    failed_variables <<- c(failed_variables, var)
    
  })
  
}


# ------------------------------------------------------------------------------
# 10. COMBINE AND SAVE RESULTS
# ------------------------------------------------------------------------------

if (length(master_results) > 0) {
  
  combined_results <- do.call(
    rbind,
    master_results
  )
  
  rownames(combined_results) <- NULL
  
  # Save RDS
  saveRDS(
    master_results,
    output_rds_path
  )
  
  # Save CSV
  write.csv(
    combined_results,
    output_csv_path,
    row.names = FALSE
  )
  
  cat("\nResults saved:\n")
  cat(output_rds_path, "\n")
  cat(output_csv_path, "\n")
  
} else {
  
  warning(
    "No variables were successfully processed. ",
    "No result files were saved."
  )
}


# ------------------------------------------------------------------------------
# 11. SAVE MODEL DIAGNOSTICS
# ------------------------------------------------------------------------------

if (length(diagnostic_results) > 0) {
  
  combined_diagnostics <- do.call(
    rbind,
    diagnostic_results
  )
  
  rownames(combined_diagnostics) <- NULL
  
  write.csv(
    combined_diagnostics,
    diagnostics_path,
    row.names = FALSE
  )
  
  cat("\nDiagnostics saved:\n")
  cat(diagnostics_path, "\n")
}


# ------------------------------------------------------------------------------
# 12. FINAL EXECUTION REPORT
# ------------------------------------------------------------------------------

cat("\n\n")
cat("==================================================================\n")
cat("ACE TWIN MODEL ANALYSIS COMPLETE\n")
cat("==================================================================\n")

cat(
  "\nSuccessfully processed:",
  length(master_results),
  "of",
  length(var_names),
  "variables\n"
)

if (length(failed_variables) > 0) {
  
  cat("\nVariables requiring attention:\n")
  print(failed_variables)
  
} else {
  
  cat("\nAll variables processed successfully.\n")
}

cat("\nOutput folder:\n")
cat(output_dir, "\n")

cat("\n==================================================================\n")
