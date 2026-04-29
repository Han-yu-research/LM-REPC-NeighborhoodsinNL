# =========================================================
# SDM (Spatial Durbin Model) workflow with Queen contiguity and eigen method
# Read data -> clean data -> construct spatial weights -> standardize variables
# -> estimate SDM -> calculate impacts -> run diagnostics
# =========================================================
# Please check the working directory and input paths before running this script.

setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2)
})

# ---------------------------
# 0) Paths and variables
# ---------------------------
csv_path <- "model/form_energy_final_cleaned_renamed.csv"
shp_path <- "shp/10580_filtered_neighborhood_boundariesx.shp"
out_dir_models <- "models"
out_dir_tables <- "tables"
out_dir_data   <- "data"
dir.create(out_dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_data,   recursive = TRUE, showWarnings = FALSE)

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
# 1) Read and merge data
# ---------------------------
df  <- read_csv(csv_path, show_col_types = FALSE)
shp <- st_read(shp_path, quiet = TRUE)

# Keep original case and remove leading/trailing spaces.
shp <- shp %>% mutate(buurtcode = trimws(as.character(buurtcode)))
df  <- df  %>% mutate(buurtcode = trimws(as.character(buurtcode)))

# If case mismatch exists, use the following lines:
# shp$buurtcode <- toupper(shp$buurtcode)
# df$buurtcode <- toupper(df$buurtcode)

if (nrow(shp) == 0) stop("The shapefile is empty. Please check the input data.")
if (any(duplicated(df$buurtcode))) stop("Duplicate 'buurtcode' values found in the CSV. Please remove duplicates or aggregate the data first.")

gdf <- shp %>% left_join(df, by = "buurtcode")
cat("Total rows after merging:", nrow(gdf), "\n")

# Check whether required columns exist
missing_cols <- setdiff(c(target_var, x_vars), names(gdf))
if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))

# ---------------------------
# 2) Clean data
# ---------------------------
classes <- sapply(st_drop_geometry(gdf)[, c(target_var, x_vars)], class)
cat("Variable type preview:\n"); print(classes)

# Convert model variables to numeric and remove non-finite values.
gdf <- gdf %>%
  dplyr::select(buurtcode, all_of(c(target_var, x_vars)), geometry) %>%
  mutate(across(all_of(c(target_var, x_vars)), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(c(target_var, x_vars)), ~ is.finite(.) & !is.na(.))) %>%
  filter(is.finite(.data[[target_var]]) & !is.na(.data[[target_var]]))

cat("Candidate sample size after cleaning:", nrow(gdf), "rows\n")

# ---------------------------
# 3) Queen contiguity and island removal
# ---------------------------
gdf <- st_make_valid(gdf)
nb0 <- poly2nb(gdf, queen = TRUE, snap = 1e-6)

iso <- which(card(nb0) == 0)
if (length(iso) > 0) {
  cat("Detected", length(iso), "isolated units. They will be removed and the neighborhood structure will be rebuilt.\n")
  gdf_no_islands <- gdf[-iso, , drop = FALSE]
  nb <- poly2nb(gdf_no_islands, queen = TRUE, snap = 1e-6)
} else {
  cat("No isolated units detected.\n")
  gdf_no_islands <- gdf
  nb <- nb0
}

# Optional: keep only the largest connected component if needed.
# library(igraph)
# g <- igraph::graph_from_adj_list(nb, mode = "undirected")
# comp <- igraph::components(g)$membership
# main <- which.max(tabulate(comp))
# keep <- which(comp == main)
# gdf_no_islands <- gdf_no_islands[keep, ]
# nb <- spdep::subset.nb(nb, keep, remap = TRUE)

# Export .gal file for reuse.
write.nb.gal(nb, file = file.path(out_dir_data, "queen_after_clean_no_islands.gal"))
cat(".gal file saved:", file.path(out_dir_data, "queen_after_clean_no_islands.gal"), "\n")

# ---------------------------
# 3.x) Save shapefile after removing islands
# ---------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Output path with timestamp to avoid overwriting files.
out_shp <- file.path(out_dir_data, sprintf("gdf_no_islands_%s.shp", timestamp))

# Write shapefile and keep attributes, including buurtcode.
st_write(gdf_no_islands, out_shp, delete_layer = TRUE)

cat("Shapefile after removing islands saved:", out_shp, "\n")
cat("Attribute columns:", paste(colnames(st_drop_geometry(gdf_no_islands)), collapse = ", "), "\n")

# ---------------------------
# 4) Create row-standardized spatial weights
# ---------------------------
lw <- nb2listw(nb, style = "W", zero.policy = FALSE)
cat("Queen spatial weights created:", length(nb), "regions.\n")

# ---------------------------
# 5) Standardize X and log-transform Y based on gdf_no_islands
# ---------------------------
log_eps <- 1e-6
df_no_geom <- st_drop_geometry(gdf_no_islands)

# Convert X variables to numeric and standardize.
X_num    <- as.data.frame(lapply(df_no_geom[, x_vars, drop = FALSE], function(v) as.numeric(v)))
X_scaled <- as.data.frame(scale(X_num, center = TRUE, scale = TRUE))

# Log-transform Y. For non-positive values, add a small constant before log transformation.
y_raw <- as.numeric(df_no_geom[[target_var]])
y_log <- ifelse(y_raw <= 0, log(y_raw + log_eps), log(y_raw))

# Assemble the sf object using the same sample as the spatial weights.
model_gdf_scaled <- gdf_no_islands
model_gdf_scaled[, x_vars]     <- X_scaled
model_gdf_scaled[[target_var]] <- y_log

cat("Model data prepared: sample size =", nrow(model_gdf_scaled),
    "; number of X variables =", length(x_vars), "\n")

# ---------------------------
# 6) Align data order with nb region.id
# ---------------------------
rid <- attr(nb, "region.id")
if (is.null(rownames(model_gdf_scaled))) {
  rownames(model_gdf_scaled) <- as.character(seq_len(nrow(model_gdf_scaled)))
}
stopifnot(all(rid %in% rownames(model_gdf_scaled)))
model_gdf_scaled <- model_gdf_scaled[rid, , drop = FALSE]

# Final dimension and order checks.
stopifnot(nrow(model_gdf_scaled) == length(nb))
stopifnot(identical(rownames(model_gdf_scaled), rid))

# ---------------------------
# 7) Diagnostics: near-zero variance and VIF
# ---------------------------
sd_vals <- apply(st_drop_geometry(model_gdf_scaled)[, x_vars, drop = FALSE], 2, sd)
nzv <- sd_vals[sd_vals < 1e-4]
if (length(nzv) > 0) { cat("Near-zero variance variables:\n"); print(nzv) }

ols_form  <- as.formula(paste(target_var, "~", paste(x_vars, collapse = " + ")))
ols_model <- lm(ols_form, data = st_drop_geometry(model_gdf_scaled))
if (requireNamespace("car", quietly = TRUE)) {
  vif_vals <- car::vif(ols_model)
  cat("\n===== VIF Check =====\n"); print(vif_vals)
  high_vif <- vif_vals[vif_vals > 5]
  if (length(high_vif) > 0) {
    cat("\nHigh VIF variables (> 5):\n"); print(high_vif)
  } else {
    cat("\nAll variables have VIF <= 5.\n")
  }
} else {
  cat("Package {car} is not installed. VIF check skipped.\n")
}

# ---------------------------
# 8) Estimate SDM
# ---------------------------
form <- as.formula(paste0(target_var, " ~ ", paste(x_vars, collapse = " + ")))
sdm_model <- spatialreg::lagsarlm(
  formula = form,
  data    = model_gdf_scaled,
  listw   = lw,
  type    = "mixed",   # SDM specification
  method  = "eigen",
  zero.policy = FALSE
)
cat("\n===== SDM Results (Queen + eigen) =====\n")
print(summary(sdm_model))

saveRDS(sdm_model, file = file.path(out_dir_models, "sdm_model_full_eigen.rds"))

# ---------------------------
# 9) Covariance matrix check and impacts
# ---------------------------
vcov_matrix <- vcov(sdm_model)
if (is.matrix(vcov_matrix) && isSymmetric(vcov_matrix)) {
  ev <- eigen(vcov_matrix, symmetric = TRUE)$values
  cat("Is the covariance matrix approximately positive definite?", min(ev) > -1e-8, "\n")
  cat("Minimum eigenvalue:", signif(min(ev), 6), "\n")
}

# ---------------------------
# 10) Model fit statistics
# ---------------------------
if (!exists("timestamp")) timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Sample size based on the final modelling data.
n_obs <- nrow(model_gdf_scaled)

# Extract LogLik, AIC, and BIC directly from the model object.
ll  <- as.numeric(logLik(sdm_model))
aic <- AIC(sdm_model)
bic <- BIC(sdm_model)

# Extract Y directly from the final data to avoid formula environment issues.
y_vec <- as.numeric(st_drop_geometry(model_gdf_scaled)[[target_var]])

# Extract residuals from the model object.
resid_vec <- residuals(sdm_model)

rss <- sum(resid_vec^2, na.rm = TRUE)
tss <- sum((y_vec - mean(y_vec, na.rm = TRUE))^2, na.rm = TRUE)
pseudo_r2_rss <- 1 - rss / tss

model_fit <- tibble::tibble(
  model = "SDM (Queen + eigen)",
  n = n_obs,
  logLik = ll,
  AIC = aic,
  BIC = bic,
  pseudo_R2_RSS = pseudo_r2_rss
)

cat("\n===== Model Fit Statistics =====\n")
print(model_fit)

fit_out <- file.path(out_dir_tables, sprintf("model_fit_stats_%s.csv", timestamp))
readr::write_csv(model_fit, fit_out)
cat("Model fit statistics saved:", fit_out, "\n")

# ---------------------------
# 11) Export impacts as a long table
# ---------------------------
set.seed(1)
R_draws <- 2000
imp <- impacts(sdm_model, listw = lw, R = R_draws, tr = TRUE)
imp_sum <- summary(imp, zstats = TRUE)

cat("\n===== Impacts (Direct / Indirect / Total) =====\n")
print(imp_sum, zstats = TRUE, short = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Save RDS files with timestamp.
saveRDS(imp,     file.path(out_dir_models, sprintf("impacts_sdm_eigen_R%d_%s.rds", R_draws, timestamp)))
saveRDS(imp_sum, file.path(out_dir_models, sprintf("impacts_sdm_eigen_R%d_summary_%s.rds", R_draws, timestamp)))

# Helper function to convert impact summary matrices into long format.
mat_to_long <- function(mat, effect_name, stat_name) {
  if (is.null(mat)) return(NULL)
  df <- as.data.frame(mat)
  df$variable <- rownames(mat)
  df$effect <- effect_name
  names(df)[1] <- stat_name
  df[, c("variable", "effect", stat_name)]
}

# imp_sum$res  : estimates for Direct, Indirect, and Total effects
# imp_sum$se   : Monte Carlo standard errors for Direct, Indirect, and Total effects
# imp_sum$zmat : z statistics for Direct, Indirect, and Total effects
res_long <- do.call(rbind, lapply(names(imp_sum$res), function(eff) {
  mat_to_long(imp_sum$res[[eff]], eff, "estimate")
}))

se_long <- NULL
if (!is.null(imp_sum$se)) {
  se_long <- do.call(rbind, lapply(names(imp_sum$se), function(eff) {
    mat_to_long(imp_sum$se[[eff]], eff, "std_error")
  }))
}

z_long <- NULL
if (!is.null(imp_sum$zmat)) {
  z_long <- do.call(rbind, lapply(names(imp_sum$zmat), function(eff) {
    mat_to_long(imp_sum$zmat[[eff]], eff, "z")
  }))
}

# Merge into one table: variable + effect + estimate + std_error + z.
impacts_long <- res_long
if (!is.null(se_long)) impacts_long <- merge(impacts_long, se_long, by = c("variable", "effect"), all.x = TRUE)
if (!is.null(z_long))  impacts_long <- merge(impacts_long, z_long,  by = c("variable", "effect"), all.x = TRUE)

# Add two-sided p-values based on z statistics.
if ("z" %in% names(impacts_long)) {
  impacts_long$p_value <- 2 * pnorm(-abs(impacts_long$z))
}

# Sort by variable and effect.
effect_order <- c("direct", "indirect", "total")
impacts_long$effect <- tolower(impacts_long$effect)
impacts_long$effect <- factor(impacts_long$effect, levels = effect_order)
impacts_long <- impacts_long[order(impacts_long$variable, impacts_long$effect), ]

# Export the final long-format impact table.
out_csv <- file.path(out_dir_tables, sprintf("impacts_sdm_eigen_long_R%d_%s.csv", R_draws, timestamp))
readr::write_csv(impacts_long, out_csv)

cat("Long-format impacts table saved:", out_csv, "\n")