# ============================================
# SDM (Matrix) — 仅保留最大连通子图 + 按组 Drop-off 比较
# KNN = 6（最大子图上重建），W 行标准化，X 标准化，Y=log(y)
# ============================================

setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(stringr); library(purrr)
})

# ---------------------------
# 0) 路径与变量、分组
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

# 分组（与主流程变量一致）
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
# 1) 读取 & 合并 & 清洗
# ---------------------------
df  <- read_csv(csv_path, show_col_types = FALSE)
shp <- st_read(shp_path, quiet = TRUE)

shp <- shp %>% mutate(buurtcode = trimws(as.character(buurtcode)))
df  <- df  %>% mutate(buurtcode = trimws(as.character(buurtcode)))

stopifnot(nrow(shp) > 0)
if (any(duplicated(df$buurtcode))) stop("CSV 中 'buurtcode' 有重复，请先去重。")

gdf <- shp %>% left_join(df, by = "buurtcode")
missing_cols <- setdiff(c(target_var, x_vars_full), names(gdf))
if (length(missing_cols) > 0) stop(paste("缺失列：", paste(missing_cols, collapse=", ")))

gdf <- gdf %>%
  dplyr::select(buurtcode, all_of(c(target_var, x_vars_full)), geometry) %>%
  mutate(across(all_of(c(target_var, x_vars_full)), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_all(all_of(c(target_var, x_vars_full)), ~ is.finite(.) & !is.na(.))) %>%
  st_make_valid()

cat("✅ 清洗后候选样本：", nrow(gdf), "\n")

# ---------------------------
# KNN = 6（直接在 gdf 上构建，没有最大子图过滤）
# ---------------------------
gdf_proj <- tryCatch(st_transform(gdf, 28992), error = function(e) gdf)
coords6  <- st_coordinates(st_centroid(st_geometry(gdf_proj)))

nbK6 <- spdep::knn2nb(spdep::knearneigh(coords6, k = 6), sym = TRUE)
stopifnot(all(spdep::card(nbK6) >= 1))

lwK6 <- nb2listw(nbK6, style = "W", zero.policy = FALSE)

cat("📌 KNN=6 连通子图数：", spdep::n.comp.nb(nbK6)$nc, 
    "（应为 1）\n")

# 导出 .gal（可复用）
gal_path <- file.path("data", sprintf("knn6_bigcomp_%s.gal", timestamp))
write.nb.gal(nbK6, file = gal_path)
cat("📁 .gal 已保存：", gal_path, "\n")

# 验证新图的连通性（应为 1）
comp_big <- spdep::n.comp.nb(nbK6)
cat("📌 重新构建后连通子图个数：", comp_big$nc, "（应为 1）\n")

# ---------------------------
# 4) 标准化 X + log(Y)（在最大子图上进行）
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

cat("✅ 模型数据准备完成：样本=", nrow(model_gdf), "；X数=", length(x_vars_full), "\n")

# ---------------------------
# 5) 通用拟合函数（Matrix）
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

# 用于阶段性保存：若文件不存在就新建，否则追加
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
# 6) 四个模型：full + drop-off 三组（阶段性保存）
# ---------------------------
cat("\n===== 拟合 Full =====\n")
m_full <- fit_sdm_matrix(x_vars_full)
write_rds(m_full, "sdm_full_matrix_bigcomp")
mr_full <- metrics_row("Full", x_vars_full, m_full); save_metrics(mr_full)

x_drop_demographic <- setdiff(x_vars_full, grp_demographic)
x_drop_building    <- setdiff(x_vars_full, grp_building)
x_drop_landscape   <- setdiff(x_vars_full, grp_landscape)

cat("\n===== 拟合 Drop Demographic =====\n")
m_drop_demographic <- fit_sdm_matrix(x_drop_demographic)
write_rds(m_drop_demographic, "sdm_drop_demographic_matrix_bigcomp")
mr_dem <- metrics_row("Drop Demographic", x_drop_demographic, m_drop_demographic); save_metrics(mr_dem)

cat("\n===== 拟合 Drop Building =====\n")
m_drop_building <- fit_sdm_matrix(x_drop_building)
write_rds(m_drop_building, "sdm_drop_building_matrix_bigcomp")
mr_bld <- metrics_row("Drop Building", x_drop_building, m_drop_building); save_metrics(mr_bld)

cat("\n===== 拟合 Drop Landscape =====\n")
m_drop_landscape <- fit_sdm_matrix(x_drop_landscape)
write_rds(m_drop_landscape, "sdm_drop_landscape_matrix_bigcomp")
mr_lm <- metrics_row("Drop Landscape", x_drop_landscape, m_drop_landscape); save_metrics(mr_lm)

# ---------------------------
# 7) 汇总比较 & 计算 ΔAIC/ΔBIC（相对 Full）
# ---------------------------
summ_tbl <- bind_rows(mr_full, mr_dem, mr_bld, mr_lm) %>%
  mutate(model = factor(model, levels = c("Full","Drop Demographic","Drop Building","Drop Landscape")))

full_AIC <- summ_tbl$AIC[summ_tbl$model == "Full"][1]
full_BIC <- summ_tbl$BIC[summ_tbl$model == "Full"][1]

summ_tbl <- summ_tbl %>%
  mutate(dAIC = AIC - full_AIC,
         dBIC = BIC - full_BIC)

print(summ_tbl)
write_csv2(summ_tbl, "sdm_matrix_bigcomp_dropoff_compare")  # 总表另存一份

# ---------------------------
# 8) 绘图：ΔAIC / ΔBIC 柱状图
# ---------------------------
# -------------------------------------------------
# 横向并排柱状图（两色分别代表 ΔAIC / ΔBIC）
# 依赖对象：summ_tbl  (含列 model, dAIC, dBIC)
# 输出：tables/<run_id>/sdm_dropoff_delta_barplot_styled.png
# -------------------------------------------------

library(ggplot2); library(dplyr); library(tidyr); library(forcats)

# ---------------------------
# 7) 绘图：ΔAIC / ΔBIC 双色横向柱状图（以500为分隔）
# ---------------------------

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

# 自动计算分隔区间（向上取整至最接近500的倍数）
max_val <- ceiling(max(plot_df$delta) / 500) * 500

p <- ggplot(plot_df, aes(x = group, y = delta, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = NULL, y = NULL,
       title = "ΔAIC / ΔBIC 相对 Full（eigen，全样本）") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.06)),
    breaks = seq(0, max_val, by = 500),     # 每500一格
    limits = c(0, max_val)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.4),
    
    # Y轴：字体放大 & 颜色黑色
    axis.text.y = element_text(size = 14, color = "black", margin = margin(r = 6)),
    
    # X轴：字体变黑
    axis.text.x = element_text(size = 12, color = "black"),
    
    # 标题字体更黑
    plot.title = element_text(hjust = 0, face = "bold", color = "black"),
    
    # 轴标题字体也变黑（虽然你这里 x、y 为空，但为了完整性）
    axis.title = element_text(color = "black")
  )



styled_plot_path <- file.path(dir_tables, "sdm_dropoff_delta_barplot_styled_eigen_all_500grid.png")
ggsave(styled_plot_path, p, width = 10, height = 4, dpi = 300)
cat("🖼 已保存：", styled_plot_path, "\n")



# 同时把用于复现的对象打包保存（阶段性成果）
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

# 🔧 与第一段代码保持一致：保存到 data/<run_id>/
dir_data_run <- file.path(out_dir_data, run_id)
dir.create(dir_data_run, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  bundle,
  file.path(dir_data_run, sprintf("bundle_matrix_bigcomp_%s.rds", timestamp))
)

cat("\n🎉 运行完成！主要输出目录：\n")
cat("   models: ", normalizePath(dir_models, winslash = "/"), "\n")
cat("   tables: ", normalizePath(dir_tables, winslash = "/"), "\n")
cat("   logs  : ", normalizePath(dir_logs,   winslash = "/"), "\n")
cat("   data  : ", normalizePath(dir_data_run, winslash = "/"), "\n")

