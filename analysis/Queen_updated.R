# =========================================================
# SDM（Spatial Durbin Model）完整流程（Queen + eigen）
# 读取 → 清洗 → 邻接 → 标准化 → SDM → impacts → 诊断
# =========================================================
# 运行前请确认工作目录与输入路径

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
# 1) 读取 & 合并（不转小写；如需匹配可统一转大写）
# ---------------------------
df  <- read_csv(csv_path, show_col_types = FALSE)
shp <- st_read(shp_path, quiet = TRUE)

# 保留原大小写；仅去空格；若两边大小写不一致，推荐用 toupper() 统一成大写
shp <- shp %>% mutate(buurtcode = trimws(as.character(buurtcode)))
df  <- df  %>% mutate(buurtcode = trimws(as.character(buurtcode)))

# 如果你知道应该统一成大写（避免 BU/bu 不一致），启用下一行：
# shp$buurtcode <- toupper(shp$buurtcode); df$buurtcode <- toupper(df$buurtcode)

if (nrow(shp) == 0) stop("Shapefile 为空，请检查数据源。")
if (any(duplicated(df$buurtcode))) stop("CSV 中存在重复的 'buurtcode'，请先去重或汇总。")

gdf <- shp %>% left_join(df, by = "buurtcode")
cat("合并后总行数：", nrow(gdf), "\n")

# 列存在性检查
missing_cols <- setdiff(c(target_var, x_vars), names(gdf))
if (length(missing_cols) > 0) stop(paste("缺失列：", paste(missing_cols, collapse=", ")))

# ---------------------------
# 2) 清洗（只保留建模所需列 + geometry）
# ---------------------------
classes <- sapply(st_drop_geometry(gdf)[, c(target_var, x_vars)], class)
cat("变量类型预览：\n"); print(classes)

# 将目标列转为数值；过滤掉非有限值
gdf <- gdf %>%
  dplyr::select(buurtcode, all_of(c(target_var, x_vars)), geometry) %>%
  mutate(across(all_of(c(target_var, x_vars)), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(c(target_var, x_vars)), ~ is.finite(.) & !is.na(.))) %>%
  filter(is.finite(.data[[target_var]]) & !is.na(.data[[target_var]]))

cat("✅ 清洗后候选样本数：", nrow(gdf), "行\n")

# ---------------------------
# 3) Queen 邻接 & 孤岛剔除
# ---------------------------
gdf <- st_make_valid(gdf)
nb0 <- poly2nb(gdf, queen = TRUE, snap = 1e-6)

iso <- which(card(nb0) == 0)
if (length(iso) > 0) {
  cat("⚠️ 检测到", length(iso), "个孤立单元，将剔除并重建邻接。\n")
  gdf_no_islands <- gdf[-iso, , drop = FALSE]
  nb <- poly2nb(gdf_no_islands, queen = TRUE, snap = 1e-6)
} else {
  cat("✅ 无孤立单元。\n")
  gdf_no_islands <- gdf
  nb <- nb0
}

# 可选：只保留最大连通子图（若碎片过多）
# library(igraph)
# g <- igraph::graph_from_adj_list(nb, mode = "undirected")
# comp <- igraph::components(g)$membership
# main <- which.max(tabulate(comp))
# keep <- which(comp == main)
# gdf_no_islands <- gdf_no_islands[keep, ]
# nb <- spdep::subset.nb(nb, keep, remap = TRUE)

# 导出 .gal（可选）
write.nb.gal(nb, file = file.path(out_dir_data, "queen_after_clean_no_islands.gal"))
cat("📁 .gal 已保存：", file.path(out_dir_data, "queen_after_clean_no_islands.gal"), "\n")

# ---------------------------
# 3.x) 保存删除孤岛后的 shp
# ---------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 输出路径（带时间戳避免覆盖）
out_shp <- file.path(out_dir_data, sprintf("gdf_no_islands_%s.shp", timestamp))

# 写出 shapefile，保留属性（包括 buurtcode）
st_write(gdf_no_islands, out_shp, delete_layer = TRUE)

cat("✅ 删除孤岛后的 shp 已保存：", out_shp, "\n")
cat("👉 属性列包括：", paste(colnames(st_drop_geometry(gdf_no_islands)), collapse = ", "), "\n")


# ---------------------------
# 4) listw（行标准化；已剔除孤岛 → zero.policy=FALSE）
# ---------------------------
lw <- nb2listw(nb, style = "W", zero.policy = FALSE)
cat("✅ Queen 权重创建完成：", length(nb), " 个区域。\n")

# ---------------------------
# 5) 标准化 X + log(Y) —— 基于 gdf_no_islands
# ---------------------------
log_eps <- 1e-6
df_no_geom <- st_drop_geometry(gdf_no_islands)

# X：数值化 + scale
X_num    <- as.data.frame(lapply(df_no_geom[, x_vars, drop = FALSE], function(v) as.numeric(v)))
X_scaled <- as.data.frame(scale(X_num, center = TRUE, scale = TRUE))

# Y：逐元素处理（<=0 加小常数），不改变大小写等
y_raw <- as.numeric(df_no_geom[[target_var]])
y_log <- ifelse(y_raw <= 0, log(y_raw + log_eps), log(y_raw))

# 组装回 sf（与 nb 同一批样本）
model_gdf_scaled <- gdf_no_islands
model_gdf_scaled[, x_vars]     <- X_scaled
model_gdf_scaled[[target_var]] <- y_log

cat("✅ 模型数据准备完成：样本数 =", nrow(model_gdf_scaled),
    "；X 列数 =", length(x_vars), "\n")

# ---------------------------
# 6) 顺序对齐（按 nb 的 region.id）
# ---------------------------
rid <- attr(nb, "region.id")
if (is.null(rownames(model_gdf_scaled))) {
  rownames(model_gdf_scaled) <- as.character(seq_len(nrow(model_gdf_scaled)))
}
stopifnot(all(rid %in% rownames(model_gdf_scaled)))
model_gdf_scaled <- model_gdf_scaled[rid, , drop = FALSE]

# 终检：维度一致
stopifnot(nrow(model_gdf_scaled) == length(nb))
stopifnot(identical(rownames(model_gdf_scaled), rid))

# ---------------------------
# 7) 诊断：近零方差 & VIF（基于最终样本）
# ---------------------------
sd_vals <- apply(st_drop_geometry(model_gdf_scaled)[, x_vars, drop = FALSE], 2, sd)
nzv <- sd_vals[sd_vals < 1e-4]
if (length(nzv) > 0) { cat("⚠️ 近零方差变量：\n"); print(nzv) }

ols_form  <- as.formula(paste(target_var, "~", paste(x_vars, collapse = " + ")))
ols_model <- lm(ols_form, data = st_drop_geometry(model_gdf_scaled))
if (requireNamespace("car", quietly = TRUE)) {
  vif_vals <- car::vif(ols_model)
  cat("\n===== VIF 检查 =====\n"); print(vif_vals)
  high_vif <- vif_vals[vif_vals > 5]
  if (length(high_vif) > 0) {
    cat("\n⚠️ 高 VIF 变量（>5）：\n"); print(high_vif)
  } else {
    cat("\n✅ 所有变量 VIF ≤ 5。\n")
  }
} else {
  cat("（提示）未安装 {car}，跳过 VIF。\n")
}

# ---------------------------
# 8) 拟合 SDM（Durbin = type="mixed"）
# ---------------------------
form <- as.formula(paste0(target_var, " ~ ", paste(x_vars, collapse = " + ")))
sdm_model <- spatialreg::lagsarlm(
  formula = form,
  data    = model_gdf_scaled,
  listw   = lw,
  type    = "mixed",   # SDM（Durbin）
  method  = "eigen",
  zero.policy = FALSE
)
cat("\n===== ✅ SDM 结果（Queen + eigen） =====\n")
print(summary(sdm_model))

saveRDS(sdm_model, file = file.path(out_dir_models, "sdm_model_full_eigen.rds"))

# ---------------------------
# 9) 协方差检查 & impacts
# ---------------------------
vcov_matrix <- vcov(sdm_model)
if (is.matrix(vcov_matrix) && isSymmetric(vcov_matrix)) {
  ev <- eigen(vcov_matrix, symmetric = TRUE)$values
  cat("✅ 协方差矩阵近似正定？", min(ev) > -1e-8, "\n")
  cat("ℹ️ 最小特征值：", signif(min(ev), 6), "\n")
}

# ---------------------------
# 8.x) 模型拟合指标（稳健版：不调用 model.frame(sdm_model)）
# ---------------------------

if (!exists("timestamp")) timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 样本量：用最终建模数据（最稳）
n_obs <- nrow(model_gdf_scaled)

# LogLik / AIC / BIC：直接从模型对象取，不会触发变量查找
ll  <- as.numeric(logLik(sdm_model))
aic <- AIC(sdm_model)
bic <- BIC(sdm_model)

# ✅ y：直接从最终数据取（避免公式环境问题）
y_vec <- as.numeric(st_drop_geometry(model_gdf_scaled)[[target_var]])

# residuals：从模型对象取
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

cat("\n===== 📌 Model Fit Statistics =====\n")
print(model_fit)

fit_out <- file.path(out_dir_tables, sprintf("model_fit_stats_%s.csv", timestamp))
readr::write_csv(model_fit, fit_out)
cat("✅ 模型拟合指标已保存：", fit_out, "\n")


# ---------------------------
# 9) impacts 长表导出（更完整：estimate + sd + z）
# ---------------------------
set.seed(1)
R_draws <- 2000
imp <- impacts(sdm_model, listw = lw, R = R_draws, tr = TRUE)
imp_sum <- summary(imp, zstats = TRUE)

cat("\n===== Impacts（Direct / Indirect / Total）=====\n")
print(imp_sum, zstats = TRUE, short = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 保存 RDS（带时间戳）
saveRDS(imp,     file.path(out_dir_models, sprintf("impacts_sdm_eigen_R%d_%s.rds", R_draws, timestamp)))
saveRDS(imp_sum, file.path(out_dir_models, sprintf("impacts_sdm_eigen_R%d_summary_%s.rds", R_draws, timestamp)))

# ---- helper：把 impacts summary 的 res / se / z 合成一张长表 ----
mat_to_long <- function(mat, effect_name, stat_name) {
  if (is.null(mat)) return(NULL)
  df <- as.data.frame(mat)
  df$variable <- rownames(mat)
  df$effect <- effect_name
  names(df)[1] <- stat_name
  df[, c("variable", "effect", stat_name)]
}

# imp_sum$res : estimate（Direct/Indirect/Total）
# imp_sum$se  : Monte Carlo std.error（Direct/Indirect/Total）
# imp_sum$zmat: z stats（Direct/Indirect/Total）
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

# 合并成一张表：variable + effect + estimate + std_error + z
impacts_long <- res_long
if (!is.null(se_long)) impacts_long <- merge(impacts_long, se_long, by = c("variable", "effect"), all.x = TRUE)
if (!is.null(z_long))  impacts_long <- merge(impacts_long, z_long,  by = c("variable", "effect"), all.x = TRUE)

# 加一个两侧 p 值（基于 z；如果没有 z 就为空）
if ("z" %in% names(impacts_long)) {
  impacts_long$p_value <- 2 * pnorm(-abs(impacts_long$z))
}

# 排序：先 variable 再 effect（Direct/Indirect/Total）
effect_order <- c("direct", "indirect", "total")
impacts_long$effect <- tolower(impacts_long$effect)
impacts_long$effect <- factor(impacts_long$effect, levels = effect_order)
impacts_long <- impacts_long[order(impacts_long$variable, impacts_long$effect), ]

# 写出一张“最终长表”
out_csv <- file.path(out_dir_tables, sprintf("impacts_sdm_eigen_long_R%d_%s.csv", R_draws, timestamp))
readr::write_csv(impacts_long, out_csv)

cat("✅ impacts 长表已保存：", out_csv, "\n")
