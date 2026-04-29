setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr)
})

# --------------------------------------------------------
# Automatically find the latest Step 1 run_id by modification time
# --------------------------------------------------------
run_dirs <- list.files("models", pattern = "^sdm_knn6_run_", full.names = TRUE)
if (length(run_dirs) == 0) stop("No folder starting with 'sdm_knn6_run_' found in the models directory.")

run_dirs <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
run_dir_latest <- run_dirs[1]
run_id <- basename(run_dir_latest)

cat("Reading SDM run:", run_id, "\n")

# --------------------------------------------------------
# Timestamp for output filenames
# --------------------------------------------------------
time_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat("Current timestamp:", time_tag, "\n")

# --------------------------------------------------------
# Step 1 input data
# --------------------------------------------------------
gdf_model <- readRDS("data/03_model_input_knn6.rds")
df_used   <- st_drop_geometry(gdf_model)

lwK6      <- readRDS("data/02_listw_knn6.rds")
sdm_model <- readRDS(file.path(run_dir_latest, "sdm_model.rds"))

# ===========================================================
# Automatically identify the dependent variable
# ===========================================================
all_cols   <- colnames(df_used)
coef_names <- names(sdm_model$coefficients)

x_candidates <- coef_names[coef_names != "(Intercept)"]
x_vars <- intersect(all_cols, x_candidates)

target_var <- setdiff(all_cols, x_vars)
target_var <- target_var[sapply(df_used[target_var], is.numeric)][1]

cat("Dependent variable identified as:", target_var, "\n")

# ===========================================================
# Automatically construct regression formula
# ===========================================================
form_sdm <- as.formula(
  paste0(target_var, " ~ ", paste(x_vars, collapse = " + "))
)

cat("\nRegression formula used:\n")
print(form_sdm)

# --------------------------------------------------------
# 1) OLS + Moran's I / LM tests
# --------------------------------------------------------
ols <- lm(form_sdm, data = df_used)

mi_ols <- moran.test(residuals(ols), lwK6)

lm_tests <- lm.RStests(
  ols, lwK6,
  test = c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")
)

lm_out <- data.frame(
  Test      = names(lm_tests),
  Statistic = sapply(lm_tests, function(x) x$statistic),
  p_value   = sapply(lm_tests, function(x) x$p.value)
)

write_csv(
  lm_out,
  paste0("models/LM_tests_OLS_", run_id, "_", time_tag, ".csv")
)

# --------------------------------------------------------
# 2) Estimate SAR and SEM models using maximum likelihood
# --------------------------------------------------------
sar_ml <- lagsarlm(form_sdm, data = df_used, listw = lwK6, method = "eigen")
sem_ml <- errorsarlm(form_sdm, data = df_used, listw = lwK6, method = "eigen")

# --------------------------------------------------------
# 3) SDM simplification tests against SAR and SEM
# --------------------------------------------------------

# SDM -> SAR test: theta = 0
test_theta <- LR.Sarlm(sdm_model, sar_ml)

# SDM -> SEM test based on log-likelihood LR test
LL_SDM <- logLik(sdm_model)
LL_SEM <- logLik(sem_ml)

LR_SEM <- 2 * (LL_SDM - LL_SEM)
df_SEM <- attr(LL_SDM, "df") - attr(LL_SEM, "df")
p_SEM  <- 1 - pchisq(LR_SEM, df_SEM)

lr_out <- data.frame(
  Test = c("SDM -> SAR (LR test)",
           "SDM -> SEM (LR test)"),
  Statistic = c(test_theta$LR, as.numeric(LR_SEM)),
  df        = c(test_theta$df, as.numeric(df_SEM)),
  p_value   = c(test_theta$p.value, as.numeric(p_SEM))
)

write_csv(
  lr_out,
  paste0("models/LR_tests_SDM_", run_id, "_", time_tag, ".csv")
)

# --------------------------------------------------------
# 4) Log-likelihood summary
# --------------------------------------------------------
loglik_out <- data.frame(
  Model = c("OLS", "SAR", "SEM", "SDM"),
  logLik = c(
    as.numeric(logLik(ols)),
    as.numeric(logLik(sar_ml)),
    as.numeric(logLik(sem_ml)),
    as.numeric(logLik(sdm_model))
  )
)

write_csv(
  loglik_out,
  paste0("models/logLik_models_", run_id, "_", time_tag, ".csv")
)

# --------------------------------------------------------
# 5) Model comparison table: logLik, AIC, and BIC
# --------------------------------------------------------
compare_tbl <- data.frame(
  Model = c("OLS", "SAR", "SEM", "SDM"),
  logLik = c(
    as.numeric(logLik(ols)),
    as.numeric(logLik(sar_ml)),
    as.numeric(logLik(sem_ml)),
    as.numeric(logLik(sdm_model))
  ),
  AIC = c(AIC(ols), AIC(sar_ml), AIC(sem_ml), AIC(sdm_model)),
  BIC = c(BIC(ols), BIC(sar_ml), BIC(sem_ml), BIC(sdm_model))
)

write_csv(
  compare_tbl,
  paste0("models/model_comparison_", run_id, "_", time_tag, ".csv")
)

cat("\nStep 2 completed. All outputs were saved with run_id and timestamp.\n")