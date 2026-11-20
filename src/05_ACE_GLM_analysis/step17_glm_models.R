# ============================================================
# BRAIN DYNAMICS, ALCOHOL USE AND ODOR IDENTIFICATION
# Combined MZ + DZ twin sample
# Primary odor measure: Odor_Unadj
#
# Usage:
#   Rscript step17_glm_models.R <mz_workbook.xlsx> <dz_workbook.xlsx> <output_dir>
#
# Output: <output_dir>/Brain_Dynamics_Alcohol_Odor_Results_Revised.xlsx
#
# Corrected version:
# - Uses lmerTest::lmer() so mixed-model coefficient p-values exist.
# - Converts every results list to a data frame before FDR correction.
# - Uses pair mean sex for pair-level models (handles opposite-sex DZ pairs).
# - Checks input columns, pair structure, model output and export.
# ============================================================

# 1. PACKAGES --------------------------------------------------
packages <- c(
  "readxl", "dplyr", "tibble", "lme4", "lmerTest",
  "sandwich", "lmtest", "openxlsx"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

# 2. FILE LOCATIONS --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript step17_glm_models.R <mz_workbook.xlsx> <dz_workbook.xlsx> <output_dir>")
}
mz_file <- args[1]
dz_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(mz_file)) stop("MZ file not found: ", mz_file)
if (!file.exists(dz_file)) stop("DZ file not found: ", dz_file)

# 3. READ DATA -------------------------------------------------
mz <- readxl::read_excel(mz_file, sheet = 1)
dz <- readxl::read_excel(dz_file, sheet = 1)

# 4. VARIABLE DEFINITIONS --------------------------------------
id_var       <- "Subject"
pair_var     <- "TwinPairID"
zygosity_var <- "ZygosityGT1"
age_var      <- "Age_in_Yrs"
sex_var      <- "Gender"
fd_var       <- "MeanFD_from_HCP_regressors"
smoking_var  <- "SSAGA_TB_Still_Smoking"

alcohol_vars <- c(
  "SSAGA_Alc_12_Drinks_Per_Day",
  "Total_Drinks_7days"
)

odor_var <- "Odor_Unadj"

brain_vars <- c(
  "FO_State 1", "FO_State 2", "FO_State 3", "FO_State 4",
  "DT_State 1", "DT_State 2", "DT_State 3", "DT_State 4",
  "AR_State 1", "AR_State 2", "AR_State 3", "AR_State 4"
)

state_labels <- c(
  "FO_State 1" = "Visual Network-A FO",
  "FO_State 2" = "DMN-C FO",
  "FO_State 3" = "Visual Network-B FO",
  "FO_State 4" = "Limbic-B FO",
  "DT_State 1" = "Visual Network-A DT",
  "DT_State 2" = "DMN-C DT",
  "DT_State 3" = "Visual Network-B DT",
  "DT_State 4" = "Limbic-B DT",
  "AR_State 1" = "Visual Network-A AR",
  "AR_State 2" = "DMN-C AR",
  "AR_State 3" = "Visual Network-B AR",
  "AR_State 4" = "Limbic-B AR"
)

# 5. VALIDATE INPUT COLUMNS ------------------------------------
required_vars <- unique(c(
  id_var, pair_var, zygosity_var, age_var, sex_var,
  fd_var, smoking_var, alcohol_vars, odor_var, brain_vars
))

missing_mz <- setdiff(required_vars, names(mz))
missing_dz <- setdiff(required_vars, names(dz))

if (length(missing_mz) > 0) {
  stop("Missing columns in MZ workbook: ",
       paste(missing_mz, collapse = ", "))
}
if (length(missing_dz) > 0) {
  stop("Missing columns in DZ workbook: ",
       paste(missing_dz, collapse = ", "))
}

mz_zyg <- unique(na.omit(as.character(mz[[zygosity_var]])))
dz_zyg <- unique(na.omit(as.character(dz[[zygosity_var]])))

if (!all(mz_zyg == "MZ")) {
  stop("Unexpected zygosity value in MZ workbook: ",
       paste(mz_zyg, collapse = ", "))
}
if (!all(dz_zyg == "DZ")) {
  stop("Unexpected zygosity value in DZ workbook: ",
       paste(dz_zyg, collapse = ", "))
}

# 6. COMBINE AND CLEAN DATA ------------------------------------
dat <- dplyr::bind_rows(mz, dz)

numeric_vars <- c(
  age_var, sex_var, fd_var, smoking_var,
  alcohol_vars, odor_var, brain_vars
)

for (v in numeric_vars) {
  dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
}

dat$SubjectID <- as.character(dat[[id_var]])
dat$PairID <- paste(
  as.character(dat[[zygosity_var]]),
  as.character(dat[[pair_var]]),
  sep = "_"
)
dat$Zygosity <- factor(
  as.character(dat[[zygosity_var]]),
  levels = c("DZ", "MZ")
)

if (any(is.na(dat$SubjectID) | dat$SubjectID == "")) {
  stop("Missing subject IDs detected.")
}
if (anyDuplicated(dat$SubjectID)) {
  stop("Duplicate subject IDs detected.")
}

# Keep only pairs with exactly two participants.
pair_sizes <- table(dat$PairID)
complete_pair_ids <- names(pair_sizes[pair_sizes == 2])
dat <- dat[dat$PairID %in% complete_pair_ids, , drop = FALSE]

# Order IDs numerically when possible, otherwise alphabetically.
id_numeric <- suppressWarnings(as.numeric(dat$SubjectID))
if (all(!is.na(id_numeric))) {
  dat <- dat[order(dat$PairID, id_numeric), , drop = FALSE]
} else {
  dat <- dat[order(dat$PairID, dat$SubjectID), , drop = FALSE]
}

if (nrow(dat) == 0) stop("No complete twin pairs found.")

cat("\nDATA CHECK\n")
cat("Participants:", nrow(dat), "\n")
cat("Complete pairs:", length(unique(dat$PairID)), "\n")
cat("MZ participants:", sum(dat$Zygosity == "MZ"), "\n")
cat("DZ participants:", sum(dat$Zygosity == "DZ"), "\n")
cat("Sex values:", paste(sort(unique(na.omit(dat[[sex_var]]))),
                         collapse = ", "), "\n")

# 7. HELPER FUNCTIONS ------------------------------------------

make_formula <- function(outcome, predictors) {
  quote_name <- function(x) paste0("`", x, "`")
  stats::as.formula(
    paste(
      quote_name(outcome), "~",
      paste(vapply(predictors, quote_name, character(1)),
            collapse = " + ")
    )
  )
}

empty_result <- function() {
  tibble::tibble(
    beta = NA_real_, SE = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, p = NA_real_
  )
}

# HC3 robust regression coefficient extraction.
get_hc3 <- function(model, predictor) {
  tryCatch({
    robust_vcov <- sandwich::vcovHC(model, type = "HC3")
    ct <- lmtest::coeftest(model, vcov. = robust_vcov)
    coefficient_names <- rownames(ct)
    
    matched <- which(coefficient_names == predictor)
    if (length(matched) == 0) {
      matched <- which(gsub("`", "", coefficient_names) == predictor)
    }
    if (length(matched) != 1) {
      stop("Could not uniquely identify coefficient ", predictor,
           ". Found: ", paste(coefficient_names, collapse = ", "))
    }
    
    beta <- unname(ct[matched, 1])
    se <- unname(ct[matched, 2])
    p_value <- unname(ct[matched, 4])
    
    tibble::tibble(
      beta = beta, SE = se,
      CI_low = beta - 1.96 * se,
      CI_high = beta + 1.96 * se,
      p = p_value
    )
  }, error = function(e) {
    message("HC3 calculation failed: ", conditionMessage(e))
    empty_result()
  })
}

# Create one row per twin pair. Pair mean sex is used rather than
# the first twin's sex, so opposite-sex DZ pairs are represented.
make_pair_data <- function(data, variables) {
  split_data <- split(data, data$PairID)
  
  rows <- lapply(split_data, function(g) {
    id_num <- suppressWarnings(as.numeric(g$SubjectID))
    if (all(!is.na(id_num))) {
      g <- g[order(id_num), , drop = FALSE]
    } else {
      g <- g[order(g$SubjectID), , drop = FALSE]
    }
    
    if (nrow(g) != 2) return(NULL)
    
    z <- list(
      PairID = as.character(g$PairID[1]),
      zygosity = droplevels(g$Zygosity[1]),
      age_mean = mean(g[[age_var]], na.rm = TRUE),
      fd_mean = mean(g[[fd_var]], na.rm = TRUE),
      smoking_mean = mean(g[[smoking_var]], na.rm = TRUE),
      sex_mean = mean(g[[sex_var]], na.rm = TRUE),
      fd_signed = g[[fd_var]][2] - g[[fd_var]][1],
      smoking_signed = g[[smoking_var]][2] - g[[smoking_var]][1]
    )
    
    for (v in variables) {
      z[[paste0(v, "_first")]] <- g[[v]][1]
      z[[paste0(v, "_second")]] <- g[[v]][2]
      z[[paste0(v, "_mean")]] <- mean(g[[v]], na.rm = TRUE)
    }
    
    as.data.frame(z, check.names = FALSE)
  })
  
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) stop("No complete twin pairs available.")
  dplyr::bind_rows(rows)
}

# Standardized absolute and signed twin differences.
add_pair_differences <- function(pair_data, variables, source_data) {
  for (v in variables) {
    sd_v <- stats::sd(source_data[[v]], na.rm = TRUE)
    
    if (!is.finite(sd_v) || sd_v == 0) {
      warning("Invalid SD for variable: ", v)
      pair_data[[paste0(v, "_abs")]] <- NA_real_
      pair_data[[paste0(v, "_signed")]] <- NA_real_
    } else {
      diff_v <- pair_data[[paste0(v, "_second")]] -
        pair_data[[paste0(v, "_first")]]
      pair_data[[paste0(v, "_abs")]] <- abs(diff_v) / sd_v
      pair_data[[paste0(v, "_signed")]] <- diff_v / sd_v
    }
  }
  pair_data
}

# Pair-level linear regression with HC3 standard errors.
fit_pair_lm <- function(d, outcome, predictors, predictor) {
  model_vars <- c(outcome, predictors)
  d <- d[stats::complete.cases(d[, model_vars, drop = FALSE]), ,
         drop = FALSE]
  
  numeric_vars_in_model <- model_vars[
    vapply(d[, model_vars, drop = FALSE], is.numeric, logical(1))
  ]
  
  if (length(numeric_vars_in_model) > 0 && nrow(d) > 0) {
    finite_rows <- Reduce(
      `&`,
      lapply(numeric_vars_in_model, function(v) is.finite(d[[v]]))
    )
    d <- d[finite_rows, , drop = FALSE]
  }
  
  if (nrow(d) <= length(predictors) + 3) {
    return(list(result = empty_result(), n = nrow(d)))
  }
  
  model <- tryCatch(
    stats::lm(make_formula(outcome, predictors), data = d),
    error = function(e) {
      message("Linear model failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(model)) {
    return(list(result = empty_result(), n = nrow(d)))
  }
  
  list(result = get_hc3(model, predictor), n = nrow(d))
}

# Within-between mixed-effects model.
# IMPORTANT: lmerTest::lmer() is used explicitly to obtain p-values.
fit_wb_model <- function(data, predictor, outcome) {
  needed <- c(
    predictor, outcome, "PairID",
    age_var, sex_var, fd_var, smoking_var, "Zygosity"
  )
  
  d <- data[
    complete.cases(data[, needed, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  # Remove rows with non-finite numeric values.
  numeric_needed <- c(predictor, outcome, age_var, sex_var, fd_var,
                      smoking_var)
  finite_rows <- Reduce(
    `&`,
    lapply(numeric_needed, function(v) is.finite(d[[v]]))
  )
  d <- d[finite_rows, , drop = FALSE]
  
  if (nrow(d) < 10 || length(unique(d$PairID)) < 5) {
    stop("Insufficient complete observations or twin pairs.")
  }
  
  d <- d %>%
    dplyr::group_by(PairID) %>%
    dplyr::mutate(
      pair_mean_predictor = mean(.data[[predictor]], na.rm = TRUE),
      within_component = .data[[predictor]] - pair_mean_predictor,
      between_component = pair_mean_predictor
    ) %>%
    dplyr::ungroup()
  
  sx <- stats::sd(d[[predictor]], na.rm = TRUE)
  sy <- stats::sd(d[[outcome]], na.rm = TRUE)
  
  if (!is.finite(sx) || sx == 0 || !is.finite(sy) || sy == 0) {
    stop("Zero or invalid variance for predictor or outcome.")
  }
  
  d$within_component <- d$within_component / sx
  d$between_component <- d$between_component / sx
  d$outcome_scaled <- d[[outcome]] / sy
  
  d$AgeModel <- d[[age_var]]
  d$SexModel <- factor(d[[sex_var]])
  d$FDModel <- d[[fd_var]]
  d$SmokingModel <- factor(d[[smoking_var]])
  d$Zygosity <- droplevels(factor(d$Zygosity))
  
  model_formula <- outcome_scaled ~
    within_component +
    between_component +
    AgeModel +
    SexModel +
    FDModel +
    SmokingModel +
    Zygosity +
    (1 | PairID)
  
  model <- lmerTest::lmer(
    model_formula,
    data = d,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
  
  co <- as.data.frame(summary(model)$coefficients)
  terms <- c("within_component", "between_component")
  
  if (!all(terms %in% rownames(co))) {
    stop("Within/between coefficients not found in model.")
  }
  
  # lmerTest normally returns this exact column. Use a safe fallback
  # in case the installed package version labels it differently.
  p_col <- grep("^Pr\\(", names(co), value = TRUE)
  if (length(p_col) == 0) {
    stop("Mixed-model p-value column is missing. Confirm lmerTest is loaded.")
  }
  
  tibble::tibble(
    effect = c("within", "between"),
    beta = unname(co[terms, "Estimate"]),
    SE = unname(co[terms, "Std. Error"]),
    CI_low = unname(co[terms, "Estimate"] -
                      1.96 * co[terms, "Std. Error"]),
    CI_high = unname(co[terms, "Estimate"] +
                       1.96 * co[terms, "Std. Error"]),
    p = unname(co[terms, p_col[1]]),
    n = nrow(d),
    pairs = length(unique(d$PairID)),
    singular = lme4::isSingular(model, tol = 1e-4)
  )
}

# Convert model-result lists into data frames and add BH FDR values.
add_fdr <- function(x, group_vars) {
  x <- dplyr::bind_rows(x)
  if (!"p" %in% names(x)) {
    stop("Results table has no p-value column; FDR cannot be calculated.")
  }
  if (length(group_vars) > 0) {
    x <- x %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::mutate(q_12 = stats::p.adjust(p, method = "BH")) %>%
      dplyr::ungroup()
  } else {
    x <- x %>% dplyr::mutate(q_12 = stats::p.adjust(p, method = "BH"))
  }
  x
}

# 8. ANALYSIS 1: ALCOHOL PAIR-LEVEL SIMILARITY ------------------
cat("\nAnalysis 1: Alcohol pair-level similarity\n")

pair_brain_alcohol <- make_pair_data(dat, c(brain_vars, alcohol_vars))
pair_brain_alcohol <- add_pair_differences(
  pair_brain_alcohol, c(brain_vars, alcohol_vars), dat
)

pair_similarity_results <- list()

for (outcome in alcohol_vars) {
  for (brain in brain_vars) {
    x <- paste0(brain, "_abs")
    y <- paste0(outcome, "_abs")
    predictors <- c(x, "age_mean", "sex_mean", "fd_mean",
                    "smoking_mean", "zygosity")
    
    ans <- fit_pair_lm(pair_brain_alcohol, y, predictors, x)
    
    pair_similarity_results[[length(pair_similarity_results) + 1]] <-
      dplyr::mutate(
        ans$result,
        brain = brain,
        brain_label = unname(state_labels[brain]),
        outcome = outcome,
        n_pairs = ans$n
      )
  }
}

pair_similarity_results <- add_fdr(pair_similarity_results, "outcome")

# 9. ANALYSIS 2: ALCOHOL SIGNED WITHIN-PAIR DIFFERENCES ---------
cat("\nAnalysis 2: Alcohol signed within-pair differences\n")

signed_results <- list()

for (outcome in alcohol_vars) {
  for (brain in brain_vars) {
    x <- paste0(brain, "_signed")
    y <- paste0(outcome, "_signed")
    predictors <- c(x, "fd_signed", "smoking_signed", "zygosity")
    
    ans <- fit_pair_lm(pair_brain_alcohol, y, predictors, x)
    
    signed_results[[length(signed_results) + 1]] <-
      dplyr::mutate(
        ans$result,
        brain = brain,
        brain_label = unname(state_labels[brain]),
        outcome = outcome,
        n_pairs = ans$n
      )
  }
}

signed_results <- add_fdr(signed_results, "outcome")

# 10. ANALYSIS 3: ALCOHOL INDIVIDUAL-LEVEL MIXED MODELS ---------
cat("\nAnalysis 3: Alcohol individual-level mixed-effects models\n")

mixed_alcohol_results <- list()

for (outcome in alcohol_vars) {
  for (brain in brain_vars) {
    row <- tryCatch({
      dplyr::mutate(
        fit_wb_model(dat, predictor = brain, outcome = outcome),
        brain = brain,
        brain_label = unname(state_labels[brain]),
        outcome = outcome
      )
    }, error = function(e) {
      message("Alcohol mixed model failed: ", brain, " / ", outcome,
              " -- ", conditionMessage(e))
      tibble::tibble(
        effect = c("within", "between"),
        beta = NA_real_, SE = NA_real_,
        CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
        n = NA_integer_, pairs = NA_integer_, singular = NA,
        brain = brain, brain_label = unname(state_labels[brain]),
        outcome = outcome
      )
    })
    mixed_alcohol_results[[length(mixed_alcohol_results) + 1]] <- row
  }
}

# 12 tests per outcome and effect (within and between separately).
mixed_alcohol_results <- add_fdr(
  mixed_alcohol_results, c("outcome", "effect")
)

# 11. ANALYSIS 4: ODOR-BRAIN INDIVIDUAL-LEVEL MIXED MODELS ------
cat("\nAnalysis 4: Odor-brain mixed-effects models\n")

odor_brain_results <- list()

for (brain in brain_vars) {
  row <- tryCatch({
    dplyr::mutate(
      fit_wb_model(dat, predictor = brain, outcome = odor_var),
      brain = brain,
      brain_label = unname(state_labels[brain]),
      odor_measure = odor_var
    )
  }, error = function(e) {
    message("Odor-brain mixed model failed: ", brain,
            " -- ", conditionMessage(e))
    tibble::tibble(
      effect = c("within", "between"),
      beta = NA_real_, SE = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      n = NA_integer_, pairs = NA_integer_, singular = NA,
      brain = brain, brain_label = unname(state_labels[brain]),
      odor_measure = odor_var
    )
  })
  odor_brain_results[[length(odor_brain_results) + 1]] <- row
}

# Separate 12-test FDR families for within and between effects.
odor_brain_results <- add_fdr(odor_brain_results, "effect")

# 12. ANALYSIS 5: ODOR-BRAIN PAIR-LEVEL SIMILARITY -------------
cat("\nAnalysis 5: Odor-brain pair-level similarity\n")

pair_odor <- make_pair_data(dat, c(brain_vars, odor_var))
pair_odor <- add_pair_differences(pair_odor, c(brain_vars, odor_var), dat)

odor_pair_results <- list()

for (brain in brain_vars) {
  x <- paste0(brain, "_abs")
  y <- paste0(odor_var, "_abs")
  predictors <- c(x, "age_mean", "sex_mean", "fd_mean",
                  "smoking_mean", "zygosity")
  
  ans <- fit_pair_lm(pair_odor, y, predictors, x)
  
  odor_pair_results[[length(odor_pair_results) + 1]] <-
    dplyr::mutate(
      ans$result,
      brain = brain,
      brain_label = unname(state_labels[brain]),
      odor_measure = odor_var,
      n_pairs = ans$n
    )
}

odor_pair_results <- add_fdr(odor_pair_results, character(0))

# 13. ANALYSIS 6: ODOR-ALCOHOL INDIVIDUAL-LEVEL MIXED MODELS --
cat("\nAnalysis 6: Odor-alcohol mixed-effects models\n")

odor_alcohol_mixed <- list()

for (outcome in alcohol_vars) {
  row <- tryCatch({
    dplyr::mutate(
      fit_wb_model(dat, predictor = odor_var, outcome = outcome),
      odor_measure = odor_var,
      outcome = outcome
    )
  }, error = function(e) {
    message("Odor-alcohol mixed model failed: ", outcome,
            " -- ", conditionMessage(e))
    tibble::tibble(
      effect = c("within", "between"),
      beta = NA_real_, SE = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      n = NA_integer_, pairs = NA_integer_, singular = NA,
      odor_measure = odor_var, outcome = outcome
    )
  })
  odor_alcohol_mixed[[length(odor_alcohol_mixed) + 1]] <- row
}

odor_alcohol_mixed <- dplyr::bind_rows(odor_alcohol_mixed)

# 14. ANALYSIS 7: ODOR-ALCOHOL PAIR-LEVEL SIMILARITY ----------
cat("\nAnalysis 7: Odor-alcohol pair-level similarity\n")

pair_odor_alcohol <- make_pair_data(dat, c(odor_var, alcohol_vars))
pair_odor_alcohol <- add_pair_differences(
  pair_odor_alcohol, c(odor_var, alcohol_vars), dat
)

odor_alcohol_pair <- list()

for (outcome in alcohol_vars) {
  x <- paste0(odor_var, "_abs")
  y <- paste0(outcome, "_abs")
  predictors <- c(x, "age_mean", "sex_mean", "fd_mean",
                  "smoking_mean", "zygosity")
  
  ans <- fit_pair_lm(pair_odor_alcohol, y, predictors, x)
  
  odor_alcohol_pair[[length(odor_alcohol_pair) + 1]] <-
    dplyr::mutate(
      ans$result,
      odor_measure = odor_var,
      outcome = outcome,
      n_pairs = ans$n
    )
}

odor_alcohol_pair <- dplyr::bind_rows(odor_alcohol_pair)

# 15. FINAL CHECKS ---------------------------------------------
results_list <- list(
  Alcohol_Pair_Similarity = pair_similarity_results,
  Alcohol_Signed_Differences = signed_results,
  Alcohol_Mixed_Effects = mixed_alcohol_results,
  Odor_Brain_Mixed = odor_brain_results,
  Odor_Brain_Pair = odor_pair_results,
  Odor_Alcohol_Mixed = odor_alcohol_mixed,
  Odor_Alcohol_Pair = odor_alcohol_pair
)

cat("\nRESULTS CHECK\n")
for (sheet_name in names(results_list)) {
  result_table <- results_list[[sheet_name]]
  cat(
    sheet_name, ": ", nrow(result_table), " rows; ",
    sum(!is.na(result_table$p)), " valid p-values\n",
    sep = ""
  )
}

expected_rows <- c(
  Alcohol_Pair_Similarity = 24,
  Alcohol_Signed_Differences = 24,
  Alcohol_Mixed_Effects = 48,
  Odor_Brain_Mixed = 24,
  Odor_Brain_Pair = 12,
  Odor_Alcohol_Mixed = 4,
  Odor_Alcohol_Pair = 2
)

for (nm in names(expected_rows)) {
  if (nrow(results_list[[nm]]) != expected_rows[[nm]]) {
    stop(
      "Unexpected number of rows in ", nm, ": ",
      nrow(results_list[[nm]]), " (expected ",
      expected_rows[[nm]], ")."
    )
  }
}

# Check that the mixed-effects tables contain the p-value column.
stopifnot("p" %in% names(mixed_alcohol_results))
stopifnot("q_12" %in% names(mixed_alcohol_results))

# 16. EXPORT ---------------------------------------------------
output_file <- file.path(
  output_dir,
  "Brain_Dynamics_Alcohol_Odor_Results_Revised.xlsx"
)

openxlsx::write.xlsx(
  results_list,
  file = output_file,
  overwrite = TRUE,
  asTable = TRUE
)

cat("\n============================================\n")
cat("ALL ANALYSES COMPLETED\n")
cat("Output workbook:\n")
cat(output_file, "\n")
cat("============================================\n")
