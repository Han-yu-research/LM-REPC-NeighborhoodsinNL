# ============================================
# SDM (Matrix) — Drop-off comparison by variable group
# KNN = 6, row-standardized W, standardized X, Y = log(y)
# ============================================

setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(stringr); library(purrr)
})

# ---------------------------
# 0) Paths, variables, and groups
# ---------------------------
csv_path <- "model/form_energy_final_cleaned_renamed.csv"
shp_path <- "shp/10580_filtered_neighborhood_boundariesx.shp"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id    <- paste0("sdm_matrix_bigcomp_", timestamp)

dir_models <- file.path("models", run_id)
dir_tables <- file.path("tables", run_id)
dir_logs   <- file.path("logs",   run_id)
dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs,   recursive = TRUE, showWarnings = FALSE)

write_rds  <- function(obj, name) saveRDS(obj, file = file.path(dir_models, paste0(name, ".rds")))
write_csv2 <- function(df, name, append = FALSE) {
  f <- file.path(dir_tables, paste0(name, ".csv"))
  if (!append || !file.exists(f)) {
    readr::write_csv(df, f)
  } else {
    readr::write_csv(df, f, append = TRUE)
  }
}
write_txt  <- function(name, expr){
  p <- file.path(dir_logs, paste0(name, ".txt"))
  zz <- capture.output({ eval.parent(substitute(expr)) })
  writeLines(zz, p); invisible(p)
}

target_var <- "ec_inten_percapita"

# Variable groups, consistent with the main modelling workflow
grp_demographic <- c('pop_density','SF_pc','inc_hi_pc','hh_size')
grp_building    <- c('ave_floor','Aplus_pc','pre2000_pc')
grp_landscape   <- c(
  'IJI_built','IJI_water',
  'AI_built','AI_water','AI_Dense_Tall','AI_Dense_Short',
  'PLAND_water','PLAND_Dense_Tall','PLAND_Dense_Short',
  'PD','SHDI'
)
x_vars_full <- c(grp_demographic, grp_building, grp_landscape)

# ---------------------------
# 1) Read, merge, and clean data
# ---------------------------
df  <- read_csv(csv_path, show_col_types = FALSE)
shp <- st_read(shp_path, quiet = TRUE)

shp <- shp %>% mutate(buurtcode = trimws(as.character(buurtcode)))
df  <- df  %>% mutate(buurtcode = trimws(as.character(buurtcode)))

stopifnot(nrow(shp) > 0)
if (any(duplicated(df$buurtcode))) stop("Duplicate 'buurtcode' values found in the CSV. Please remove duplicates first.")

gdf <- shp %>% left_join(df, by = "buurtcode")
missing_cols <- setdiff(c(target_var, x_vars_full), names(gdf))
if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))

gdf <- gdf %>%
  dplyr::select(buurtcode, all_of(c(target_var, x_vars_full)), geometry) %>%
  mutate(across(all_of(c(target_var, x_vars_full)), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(c(target_var, x_vars_full)), ~ is.finite(.) & !is.na(.))) %>%
  st_make_valid()

cat("Candidate sample size after cleaning:", nrow(gdf), "\n")

# ---------------------------
# 2) Construct KNN = 6 spatial weights
# ---------------------------
gdf_proj <- tryCatch(st_transform(gdf, 28992), error = function(e) gdf)
coords6  <- st_coordinates(st_centroid(st_geometry(gdf_proj)))

nbK6 <- spdep::knn2nb(spdep::knearneigh(coords6, k = 6), sym = TRUE)
stopifnot(all(spdep::card(nbK6) >= 1))

lwK6 <- nb2listw(nbK6, style = "W", zero.policy = FALSE)

cat("Number of connected components for KNN = 6:", spdep::n.comp.nb(nbK6)$nc, 
    "(expected: 1)\n")

# Export .gal file for reuse
gal_path <- file.path("data", sprintf("knn6_bigcomp_%s.gal", timestamp))
write.nb.gal(nbK6, file = gal_path)
cat(".gal file saved:", gal_path, "\n")

# Check connectivity
comp_big <- spdep::n.comp.nb(nbK6)
cat("Number of connected components after reconstruction:", comp_big$nc, "(expected: 1)\n")

# ---------------------------
# 3) Standardize X and log-transform Y
# ---------------------------
log_eps <- 1e-6
df_ng   <- st_drop_geometry(gdf)

X_num    <- as.data.frame(lapply(df_ng[, x_vars_full, drop = FALSE], function(v) as.numeric(v)))
X_scaled <- as.data.frame(scale(X_num, center = TRUE, scale = TRUE))

y_raw <- as.numeric(df_ng[[target_var]])
y_log <- ifelse(y_raw <= 0, log(y_raw + log_eps), log(y_raw))

model_gdf <- gdf
model_gdf[, x_vars_full] <- X_scaled
model_gdf[[target_var]]  <- y_log
if (is.null(rownames(model_gdf))) rownames(model_gdf) <- as.character(seq_len(nrow(model_gdf)))

cat("Model data prepared: sample size =", nrow(model_gdf), "; number of X variables =", length(x_vars_full), "\n")

# ---------------------------
# 4) General SDM fitting function using Matrix method
# ---------------------------
fit_sdm_matrix <- function(vars) {
  stopifnot(length(vars) >= 1)
  form <- as.formula(paste(target_var, "~", paste(vars, collapse = " + ")))
  spatialreg::lagsarlm(
    formula = form,
    data    = model_gdf,
    listw   = lwK6,
    type    = "mixed",
    method  = "Matrix",
    zero.policy = FALSE
  )
}

metrics_row <- function(model_name, vars, fit) {
  tibble::tibble(
    model = model_name,
    k     = length(vars),
    logLik = as.numeric(logLik(fit)),
    AIC    = AIC(fit),
    BIC    = BIC(fit)
  )
}

# Save metrics incrementally: create the file if it does not exist, otherwise append
metrics_file <- "sdm_matrix_bigcomp_metrics"
save_metrics <- function(row_df) {
  file_path <- file.path(dir_tables, paste0(metrics_file, ".csv"))
  if (!file.exists(file_path)) {
    readr::write_csv(row_df, file_path)
  } else {
    readr::write_csv(row_df, file_path, append = TRUE)
  }
}

# ---------------------------
# 5) Estimate four models: full model and three drop-off models
# ---------------------------
cat("\n===== Fitting Full model =====\n")
m_full <- fit_sdm_matrix(x_vars_full)
write_rds(m_full, "sdm_full_matrix_bigcomp")
mr_full <- metrics_row("Full", x_vars_full, m_full); save_metrics(mr_full)

x_drop_demographic <- setdiff(x_vars_full, grp_demographic)
x_drop_building    <- setdiff(x_vars_full, grp_building)
x_drop_landscape   <- setdiff(x_vars_full, grp_landscape)

cat("\n===== Fitting Drop Demographic model =====\n")
m_drop_demographic <- fit_sdm_matrix(x_drop_demographic)
write_rds(m_drop_demographic, "sdm_drop_demographic_matrix_bigcomp")
mr_dem <- metrics_row("Drop Demographic", x_drop_demographic, m_drop_demographic); save_metrics(mr_dem)

cat("\n===== Fitting Drop Building model =====\n")
m_drop_building <- fit_sdm_matrix(x_drop_building)
write_rds(m_drop_building, "sdm_drop_building_matrix_bigcomp")
mr_bld <- metrics_row("Drop Building", x_drop_building, m_drop_building); save_metrics(mr_bld)

cat("\n===== Fitting Drop Landscape model =====\n")
m_drop_landscape <- fit_sdm_matrix(x_drop_landscape)
write_rds(m_drop_landscape, "sdm_drop_landscape_matrix_bigcomp")
mr_lm <- metrics_row("Drop Landscape", x_drop_landscape, m_drop_landscape); save_metrics(mr_lm)

# ---------------------------
# 6) Summarize model comparison and calculate delta AIC/BIC relative to the full model
# ---------------------------
summ_tbl <- bind_rows(mr_full, mr_dem, mr_bld, mr_lm) %>%
  mutate(model = factor(model, levels = c("Full","Drop Demographic","Drop Building","Drop Landscape")))

full_AIC <- summ_tbl$AIC[summ_tbl$model == "Full"][1]
full_BIC <- summ_tbl$BIC[summ_tbl$model == "Full"][1]

summ_tbl <- summ_tbl %>%
  mutate(dAIC = AIC - full_AIC,
         dBIC = BIC - full_BIC)

print(summ_tbl)
write_csv2(summ_tbl, "sdm_matrix_bigcomp_dropoff_compare")

# ---------------------------
# 7) Plot delta AIC/BIC bar chart
# ---------------------------
library(ggplot2); library(dplyr); library(tidyr); library(forcats)

plot_df <- summ_tbl %>%
  filter(model != "Full") %>%
  mutate(group = case_when(
    model == "Drop Demographic" ~ "Demographic\nindicators",
    model == "Drop Building"    ~ "Building\ncharacteristics",
    model == "Drop Landscape"   ~ "Landscape\nmetrics"
  )) %>%
  select(group, dAIC, dBIC) %>%
  pivot_longer(cols = c(dAIC, dBIC), names_to = "metric", values_to = "delta") %>%
  mutate(metric = ifelse(metric == "dAIC", "AIC", "BIC"))

cols <- c(AIC = "#A7CAE5", BIC = "#E8A089")

plot_df <- plot_df %>%
  group_by(group) %>% mutate(delta_max = max(delta)) %>% ungroup() %>%
  mutate(group = fct_reorder(group, delta_max))

# Automatically calculate the upper limit rounded to the nearest 500
max_val <- ceiling(max(plot_df$delta) / 500) * 500

p <- ggplot(plot_df, aes(x = group, y = delta, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = NULL, y = NULL,
       title = "Delta AIC / BIC relative to the full model") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.06)),
    breaks = seq(0, max_val, by = 500),
    limits = c(0, max_val)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.4),
    axis.text.y = element_text(size = 14, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0, face = "bold", color = "black"),
    axis.title = element_text(color = "black")
  )

styled_plot_path <- file.path(dir_tables, "sdm_dropoff_delta_barplot_styled_eigen_all_500grid.png")
ggsave(styled_plot_path, p, width = 10, height = 4, dpi = 300)
cat("Plot saved:", styled_plot_path, "\n")

# Save reproducibility bundle
bundle <- list(
  models = list(
    full = m_full,
    drop_demographic = m_drop_demographic,
    drop_building    = m_drop_building,
    drop_landscape   = m_drop_landscape
  ),
  lw   = lwK6,
  nb   = nbK6,
  data = model_gdf,
  groups = list(
    demographic = grp_demographic,
    building    = grp_building,
    landscape   = grp_landscape
  ),
  compare_table = summ_tbl,
  meta = list(timestamp = timestamp, gal = gal_path)
)

# Save bundle to data/<run_id>/
dir_data_run <- file.path(out_dir_data, run_id)
dir.create(dir_data_run, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  bundle,
  file.path(dir_data_run, sprintf("bundle_matrix_bigcomp_%s.rds", timestamp))
)

cat("\nRun completed. Main output directories:\n")
cat("   models: ", normalizePath(dir_models, winslash = "/"), "\n")
cat("   tables: ", normalizePath(dir_tables, winslash = "/"), "\n")
cat("   logs  : ", normalizePath(dir_logs,   winslash = "/"), "\n")
cat("   data  : ", normalizePath(dir_data_run, winslash = "/"), "\n")