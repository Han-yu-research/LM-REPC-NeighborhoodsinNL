# =========================================================
# SDM (Spatial Durbin Model) workflow with KNN = 6 and eigen method
# Read data -> clean data -> construct KNN = 6 weights -> standardize variables
# -> estimate SDM -> calculate impacts -> diagnostics -> save and plot residuals
# =========================================================

setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(writexl)
})

# ---------------------------
# 0) Paths and variables
# ---------------------------
csv_path <- "model/form_energy_final_cleaned_renamed.csv"
shp_path <- "shp/10580_filtered_neighborhood_boundariesx.shp"
out_dir_models <- "models"
out_dir_tables <- "tables"
out_dir_data   <- "data"
dir.create(out_dir_models, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_data,   showWarnings = FALSE, recursive = TRUE)

target_var <- "ec_inten_percapita"
x_vars <- c(
  'Aplus_pc','pop_density','SF_pc',
  'pre2000_pc','inc_hi_pc','ave_floor', 'hh_size',
  'IJI_water','AI_water','IJI_built','AI_built',
  'AI_Dense_Tall','AI_Dense_Short',
  'PLAND_water','PLAND_Dense_Tall','PLAND_Dense_Short',
  'PD','SHDI'
)

# ---------------------------
# 0.1) Unified run_id and output helper functions
# ---------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0("sdm_knn6_run_", timestamp)

dir_models_run <- file.path(out_dir_models, run_id)
dir_tables_run <- file.path(out_dir_tables, run_id)
dir_logs_run   <- file.path("logs", run_id)
dir.create(dir_models_run, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables_run, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs_run,   recursive = TRUE, showWarnings = FALSE)

write_rds  <- function(obj, name) saveRDS(obj, file = file.path(dir_models_run, paste0(name, ".rds")))
write_csv2 <- function(df, name)  readr::write_csv(df, file.path(dir_tables_run, paste0(name, ".csv")))
write_txt  <- function(name, expr){
  p <- file.path(dir_logs_run, paste0(name, ".txt"))
  zz <- capture.output({ eval.parent(substitute(expr)) })
  writeLines(zz, p); invisible(p)
}

# Convert impacts summary output into long format
as_long_res <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.data.frame(x)) {
    out <- lapply(names(x), function(eff){
      mat <- x[[eff]]; if (is.null(mat)) return(NULL)
      df <- as.data.frame(mat)
      df$variable <- rownames(mat)
      df$effect <- eff
      df
    })
    do.call(rbind, out)
  } else {
    df <- as.data.frame(x); df$variable <- rownames(x); df
  }
}

# ---------------------------
# 1) Read and merge data
# ---------------------------
df  <- read_csv(csv_path, show_col_types = FALSE)
shp <- st_read(shp_path, quiet = TRUE)

shp <- shp %>% mutate(buurtcode = trimws(as.character(buurtcode)))
df  <- df  %>% mutate(buurtcode = trimws(as.character(buurtcode)))

if (nrow(shp) == 0) stop("The shapefile is empty. Please check the input data.")
if (any(duplicated(df$buurtcode))) stop("Duplicate 'buurtcode' values found in the CSV. Please remove duplicates or aggregate the data first.")

gdf <- shp %>% left_join(df, by = "buurtcode")
cat("Total rows after merging:", nrow(gdf), "\n")

missing_cols <- setdiff(c(target_var, x_vars), names(gdf))
if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))

# ---------------------------
# 2) Clean data
# ---------------------------
classes <- sapply(st_drop_geometry(gdf)[, c(target_var, x_vars)], class)
cat("Variable type preview:\n"); print(classes); write_txt("var_classes_print", { print(classes) })

gdf <- gdf %>%
  dplyr::select(buurtcode, all_of(c(target_var, x_vars, "ste_mvs")), geometry) %>%
  mutate(across(all_of(c(target_var, x_vars)), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(c(target_var, x_vars)), ~ is.finite(.) & !is.na(.))) %>%
  filter(is.finite(.data[[target_var]]) & !is.na(.data[[target_var]]))

gdf <- st_make_valid(gdf)
cat("Candidate sample size after cleaning:", nrow(gdf), "rows\n")

# ---------------------------
# 3) Construct KNN = 6 spatial weights and export .gal file
# ---------------------------
gdf_proj <- tryCatch(st_transform(gdf, 28992), error = function(e) gdf)
coords6  <- st_coordinates(st_centroid(st_geometry(gdf_proj)))
nbK6     <- spdep::knn2nb(spdep::knearneigh(coords6, k = 6), sym = TRUE)

stopifnot(all(spdep::card(nbK6) >= 1))

gal_path <- file.path(out_dir_data, "knn6_after_clean.gal")
write.nb.gal(nbK6, file = gal_path)
cat(".gal file saved:", gal_path, "\n")

# ---------------------------
# 4) Create row-standardized spatial weights
# ---------------------------
lwK6 <- nb2listw(nbK6, style = "W", zero.policy = FALSE)
cat("KNN = 6 spatial weights created:", length(nbK6), "regions; average degree:",
    round(mean(spdep::card(nbK6)), 2), "\n")

# ---------------------------
# 5) Standardize X and log-transform Y
# ---------------------------
log_eps <- 1e-6
df_no_geom <- st_drop_geometry(gdf)

X_num    <- as.data.frame(lapply(df_no_geom[, x_vars, drop = FALSE], function(v) as.numeric(v)))
X_scaled <- as.data.frame(scale(X_num, center = TRUE, scale = TRUE))

y_raw <- as.numeric(df_no_geom[[target_var]])
y_log <- ifelse(y_raw <= 0, log(y_raw + log_eps), log(y_raw))

model_gdf_scaled <- gdf
model_gdf_scaled[, x_vars]     <- X_scaled
model_gdf_scaled[[target_var]] <- y_log

if (is.null(rownames(model_gdf_scaled))) {
  rownames(model_gdf_scaled) <- as.character(seq_len(nrow(model_gdf_scaled)))
}

cat("Model data prepared: sample size =", nrow(model_gdf_scaled),
    "; number of X variables =", length(x_vars), "\n")

# ---------------------------
# 6) Diagnostics: near-zero variance and VIF
# ---------------------------
sd_vals <- apply(st_drop_geometry(model_gdf_scaled)[, x_vars, drop = FALSE], 2, sd)
nzv <- sd_vals[sd_vals < 1e-4]

if (length(nzv) > 0) {
  cat("Near-zero variance variables:\n")
  print(nzv)
}

ols_form  <- as.formula(paste(target_var, "~", paste(x_vars, collapse = " + ")))
ols_model <- lm(ols_form, data = st_drop_geometry(model_gdf_scaled))

if (requireNamespace("car", quietly = TRUE)) {
  vif_vals <- car::vif(ols_model)
  cat("\n===== VIF Check =====\n")
  print(vif_vals)
  
  high_vif <- vif_vals[vif_vals > 5]
  if (length(high_vif) > 0) {
    cat("\nHigh VIF variables (> 5):\n")
    print(high_vif)
  } else {
    cat("\nAll variables have VIF <= 5.\n")
  }
} else {
  cat("Package {car} is not installed. VIF check skipped.\n")
}

if (length(nzv) > 0) {
  nzv_tbl <- data.frame(variable = names(nzv), sd = as.numeric(nzv))
  write_csv2(nzv_tbl, "near_zero_variance_vars")
}

if (exists("vif_vals")) {
  vif_tbl <- data.frame(variable = names(vif_vals), vif = as.numeric(vif_vals))
  write_csv2(vif_tbl, "vif_values")
}

# ---------------------------
# 7) Estimate SDM
# ---------------------------
form <- as.formula(paste0(target_var, " ~ ", paste(x_vars, collapse = " + ")))

sdm_model <- spatialreg::lagsarlm(
  formula = form,
  data = model_gdf_scaled,
  listw = lwK6,
  type = "mixed",
  method = "eigen",
  zero.policy = FALSE
)

cat("\n===== SDM Results (KNN = 6 + eigen) =====\n")
print(summary(sdm_model))

write_rds(sdm_model, "sdm_model_knn6_eigen")

sdm_coef <- data.frame(term = names(coef(sdm_model)), estimate = as.numeric(coef(sdm_model)))
write_csv2(sdm_coef, "sdm_knn6_coefficients")

sdm_vcov <- as.data.frame(vcov(sdm_model))
sdm_vcov$..row <- rownames(vcov(sdm_model))
write_csv2(sdm_vcov, "sdm_knn6_vcov")

write_txt("sdm_knn6_summary", { print(summary(sdm_model)) })

# ---------------------------
# 8) Covariance matrix check and impacts
# ---------------------------
vcov_matrix <- vcov(sdm_model)

if (is.matrix(vcov_matrix) && isSymmetric(vcov_matrix)) {
  ev <- eigen(vcov_matrix, symmetric = TRUE)$values
  cat("Is the covariance matrix approximately positive definite?", min(ev) > -1e-8, "\n")
  cat("Minimum eigenvalue:", signif(min(ev), 6), "\n")
}

set.seed(1)
R_draws <- 2000
imp <- impacts(sdm_model, listw = lwK6, R = R_draws, tr = TRUE)
imp_sum <- summary(imp, zstats = TRUE)

cat("\n===== Impacts (Direct / Indirect / Total) =====\n")
print(imp_sum, zstats = TRUE, short = TRUE)

# Export impacts with p-values
export_impacts_with_p <- function(imp_sum, pmat, filename = "impacts_with_p.xlsx") {
  output_dir <- "output"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  out_path <- file.path(output_dir, filename)
  
  clean_p <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("<", "", x)
    as.numeric(x)
  }
  
  p_df <- as.data.frame(pmat)
  p_df$variable <- rownames(p_df)
  p_df$Direct   <- clean_p(p_df$Direct)
  p_df$Indirect <- clean_p(p_df$Indirect)
  p_df$Total    <- clean_p(p_df$Total)
  
  get_mean <- function(x) x$statistics[, "Mean"]
  
  df <- data.frame(
    variable = rownames(imp_sum$direct$statistics),
    direct_mean   = get_mean(imp_sum$direct),
    indirect_mean = get_mean(imp_sum$indirect),
    total_mean    = get_mean(imp_sum$total),
    direct_p    = p_df$Direct,
    indirect_p  = p_df$Indirect,
    total_p     = p_df$Total
  )
  
  stars <- function(p) {
    ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
  }
  
  df$direct_star   <- stars(df$direct_p)
  df$indirect_star <- stars(df$indirect_p)
  df$total_star    <- stars(df$total_p)
  df$direct_effect   <- paste0(round(df$direct_mean, 4),   df$direct_star)
  df$indirect_effect <- paste0(round(df$indirect_mean, 4), df$indirect_star)
  df$total_effect    <- paste0(round(df$total_mean, 4),    df$total_star)
  
  writexl::write_xlsx(df, out_path)
  cat("\nImpacts and p-values saved to:", normalizePath(out_path), "\n\n")
  
  return(df)
}

df_final <- export_impacts_with_p(imp_sum, pmat = imp_sum$pzmat, filename = "sdm_impacts_with_p.xlsx")

# ---------------------------
# 9) Model fit statistics and residual saving
# ---------------------------
y_obs <- as.numeric(sdm_model$y)
resid <- as.numeric(sdm_model$residuals)
n_obs <- length(y_obs)
k_par <- sdm_model$parameters
ll    <- sdm_model$LL

# Pseudo R2
SSE <- sum(resid^2)
SST <- sum((y_obs - mean(y_obs))^2)
R2_pseudo <- 1 - SSE / SST

metrics <- data.frame(
  logLik     = ll,
  AIC        = -2 * ll + 2 * k_par,
  BIC        = -2 * ll + log(n_obs) * k_par,
  R2_pseudo  = R2_pseudo,
  n          = n_obs,
  k          = k_par
)

print(metrics)
write_csv2(metrics, "sdm_model_metrics")
write_rds(metrics, "sdm_model_metrics")

# ---------------------------
# 10) Save residuals to model_gdf_scaled
# ---------------------------
model_gdf_scaled$residual <- resid

write_rds(model_gdf_scaled, file = file.path(out_dir_data, "03_model_input_knn6_with_resid.rds"))
write_csv2(st_drop_geometry(model_gdf_scaled), "03_model_input_knn6_with_resid")

cat("Residuals saved to model_gdf_scaled.\n")

# ---------------------------
# 11) Plot residuals against ste_mvs
# ---------------------------
if (!"ste_mvs" %in% colnames(model_gdf_scaled)) {
  stop("Column 'ste_mvs' is not found in the data.")
}

df_violin <- data.frame(
  ste_mvs = factor(model_gdf_scaled$ste_mvs),
  residual = model_gdf_scaled$residual
)

p <- ggplot(df_violin, aes(x = ste_mvs, y = residual, fill = ste_mvs)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(
    x = "ste_mvs",
    y = "Residual (SDM)",
    title = "Residual vs ste_mvs (Violin Plot)"
  ) +
  scale_fill_brewer(palette = "Set2")

plot(p)

ggsave(
  filename = file.path(out_dir_tables, "residual_vs_ste_mvs_violin.png"),
  plot = p,
  width = 6,
  height = 4
)

cat("Residual vs ste_mvs violin plot saved.\n")

# ---------------------------
# 12) Save input files required by Step 2
# ---------------------------
saveRDS(model_gdf_scaled, file = file.path(out_dir_data, "03_model_input_knn6.rds"))
saveRDS(lwK6, file = file.path(out_dir_data, "02_listw_knn6.rds"))
saveRDS(sdm_model, file = file.path(dir_models_run, "sdm_model.rds"))

cat(
  "Files required by Step 2 saved:\n",
  "   - data/03_model_input_knn6.rds\n",
  "   - data/02_listw_knn6.rds\n",
  "   - models/", run_id, "/sdm_model.rds\n",
  sep = ""
)