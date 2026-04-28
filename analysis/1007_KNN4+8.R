# =========================================================
# 批量运行 SDM（Spatial Durbin Model）
# KNN = 4, 6, 8 版本，自动比较并保存结果
# =========================================================

setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2)
})


# ---------------------------
# 0) 路径与变量设定
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
# 1) 读取 & 合并
# ---------------------------
df  <- read_csv(csv_path, show_col_types = FALSE)
shp <- st_read(shp_path, quiet = TRUE)

shp <- shp %>% mutate(buurtcode = trimws(as.character(buurtcode)))
df  <- df  %>% mutate(buurtcode = trimws(as.character(buurtcode)))
gdf <- shp %>% left_join(df, by = "buurtcode") %>%
  dplyr::select(buurtcode, all_of(c(target_var, x_vars)), geometry) %>%
  mutate(across(all_of(c(target_var, x_vars)), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(c(target_var, x_vars)), ~ is.finite(.) & !is.na(.)))
gdf <- st_make_valid(gdf)

# ---------------------------
# ---------------------------
# 2) 函数定义（替换为此版本）
# ---------------------------
run_sdm_for_k <- function(k, gdf, target_var, x_vars) {
  cat("\n============================\n")
  cat("🔹 开始运行 SDM（KNN =", k, "）\n")
  cat("============================\n")
  
  # ---- 创建权重矩阵 ----
  gdf_proj <- tryCatch(st_transform(gdf, 28992), error = function(e) gdf)
  coords <- st_coordinates(st_centroid(st_geometry(gdf_proj)))
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k), sym = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = FALSE)
  gal_path <- file.path(out_dir_data, paste0("knn", k, "_after_clean.gal"))
  write.nb.gal(nb, gal_path)
  cat("✅ 权重矩阵创建完成，保存至：", gal_path, "\n")
  
  # ---- 数据标准化 ----
  df_no_geom <- st_drop_geometry(gdf)
  X_num    <- as.data.frame(lapply(df_no_geom[, x_vars, drop = FALSE], as.numeric))
  X_scaled <- as.data.frame(scale(X_num, center = TRUE, scale = TRUE))
  y_raw <- as.numeric(df_no_geom[[target_var]])
  y_log <- ifelse(y_raw <= 0, log(y_raw + 1e-6), log(y_raw))
  
  gdf_scaled <- gdf
  gdf_scaled[, x_vars]     <- X_scaled
  gdf_scaled[[target_var]] <- y_log
  
  # ---- 模型拟合 ----
  form <- as.formula(paste(target_var, "~", paste(x_vars, collapse = " + ")))
  sdm_model <- spatialreg::lagsarlm(
    formula = form,
    data    = gdf_scaled,
    listw   = lw,
    type    = "mixed",
    method  = "eigen",
    zero.policy = FALSE
  )
  
  cat("\n📈 SDM 模型（KNN =", k, "）结果：\n")
  print(summary(sdm_model))
  
  # ---- impacts（如需保留就留着）----
  set.seed(1)
  R_draws <- 2000
  imp <- impacts(sdm_model, listw = lw, R = R_draws, tr = TRUE)
  imp_sum <- summary(imp, zstats = TRUE)
  cat("\n📊 Impacts 汇总（KNN =", k, "）\n")
  print(imp_sum, zstats = TRUE, short = TRUE)
  
  # ----【新增】关键指标：使用模型内置 y/fitted/residuals ----
  y_obs  <- as.numeric(sdm_model$y)
  y_hat  <- as.numeric(sdm_model$fitted.values)
  resid  <- as.numeric(sdm_model$residuals)
  n      <- length(y_obs)
  
  SSE <- sum(resid^2)
  SST <- sum((y_obs - mean(y_obs))^2)
  R2_corr <- as.numeric(cor(y_hat, y_obs)^2)
  R2_sse  <- 1 - SSE / SST
  
  px <- length(x_vars)
  k_eff <- 2*px + 3  # 近似：截距 + rho + sigma^2 + X 和 WX
  Adj_R2 <- 1 - (1 - R2_sse) * (n - 1) / (n - k_eff - 1)
  
  ll_sdm  <- as.numeric(logLik(sdm_model))
  aic_sdm <- AIC(sdm_model)
  bic_sdm <- BIC(sdm_model)
  
  cat("\n===== 📌 关键指标（KNN =", k, ") =====\n")
  cat(sprintf("R² (corr^2)      : %.5f\n", R2_corr))
  cat(sprintf("R² (1 - SSE/SST) : %.5f\n", R2_sse))
  cat(sprintf("Adjusted R² (≈)  : %.5f  [k=%d]\n", Adj_R2, k_eff))
  cat(sprintf("Log Likelihood   : %.2f\n", ll_sdm))
  cat(sprintf("AIC / BIC        : %.0f / %.0f\n", aic_sdm, bic_sdm))
  
  # ----【新增】保存模型（RDS）并在全局命名一份（可选）----
  rds_path <- file.path(out_dir_models, paste0("sdm_knn", k, ".rds"))
  saveRDS(sdm_model, rds_path)
  assign(paste0("sdm_model_knn", k), sdm_model, envir = .GlobalEnv)
  cat("💾 模型已保存：", rds_path, "\n")
  
  # ---- 汇总指标（返回到比较表）----
  result <- data.frame(
    KNN     = k,
    logLik  = ll_sdm,
    AIC     = aic_sdm,
    BIC     = bic_sdm,
    rho     = sdm_model$rho,
    n       = n,
    R2_corr = R2_corr,
    R2_sse  = R2_sse,
    Adj_R2  = Adj_R2,
    stringsAsFactors = FALSE
  )
  return(result)
}


# ---------------------------
# 3) 批量运行 KNN = 4, 8
# ---------------------------
knn_list <- c(4, 8)
results_all <- lapply(knn_list, function(k) run_sdm_for_k(k, gdf, target_var, x_vars))
results_df <- do.call(rbind, results_all)

# ---------------------------
# 4) 保存与展示结果
# ---------------------------
write_csv(results_df, file.path(out_dir_tables, "sdm_knn_results_compare.csv"))
cat("\n✅ 所有模型已运行完成。\n")
cat("📁 结果表已保存： tables/sdm_knn_results_compare.csv\n")

print(results_df)


