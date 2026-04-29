setwd("D:/1st/code/project")

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg)
  library(dplyr); library(readr)
})

# --------------------------------------------------------
# 自动找到 Step1 最新 run_id（按修改时间）
# --------------------------------------------------------
run_dirs <- list.files("models", pattern="^sdm_knn6_run_", full.names = TRUE)
if (length(run_dirs) == 0) stop("❌ models/ 下找不到 sdm_knn6_run_ 开头的文件夹")

run_dirs <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
run_dir_latest <- run_dirs[1]
run_id <- basename(run_dir_latest)

cat("📁 读取的 SDM run：", run_id, "\n")

# --------------------------------------------------------
# 时间戳（用于结果文件命名）
# --------------------------------------------------------
time_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat("🕒 当前时间戳：", time_tag, "\n")

# --------------------------------------------------------
# Step1 输入数据
# --------------------------------------------------------
gdf_model <- readRDS("data/03_model_input_knn6.rds")
df_used   <- st_drop_geometry(gdf_model)

lwK6      <- readRDS("data/02_listw_knn6.rds")
sdm_model <- readRDS(file.path(run_dir_latest, "sdm_model.rds"))

# ===========================================================
# 自动识别因变量
# ===========================================================
all_cols   <- colnames(df_used)
coef_names <- names(sdm_model$coefficients)

x_candidates <- coef_names[coef_names != "(Intercept)"]
x_vars <- intersect(all_cols, x_candidates)

target_var <- setdiff(all_cols, x_vars)
target_var <- target_var[sapply(df_used[target_var], is.numeric)][1]

cat("✔ 因变量识别为：", target_var, "\n")

# ===========================================================
# 自动构建公式
# ===========================================================
form_sdm <- as.formula(
  paste0(target_var, " ~ ", paste(x_vars, collapse = " + "))
)

cat("\n📌 使用的回归公式：\n")
print(form_sdm)

# --------------------------------------------------------
# 1) OLS + Moran / LM tests
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
# 2) SAR / SEM（ML）
# --------------------------------------------------------
sar_ml <- lagsarlm(form_sdm, data = df_used, listw = lwK6, method = "eigen")
sem_ml <- errorsarlm(form_sdm, data = df_used, listw = lwK6, method = "eigen")

# --------------------------------------------------------
# 3) SDM 化简检验（SAR / SEM）
# --------------------------------------------------------

# SDM → SAR（theta = 0）
test_theta <- LR.Sarlm(sdm_model, sar_ml)

# SDM → SEM（logLik LR）
LL_SDM <- logLik(sdm_model)
LL_SEM <- logLik(sem_ml)

LR_SEM <- 2 * (LL_SDM - LL_SEM)
df_SEM <- attr(LL_SDM, "df") - attr(LL_SEM, "df")
p_SEM  <- 1 - pchisq(LR_SEM, df_SEM)

lr_out <- data.frame(
  Test = c("SDM → SAR (LR test)",
           "SDM → SEM (LR test)"),
  Statistic = c(test_theta$LR, as.numeric(LR_SEM)),
  df        = c(test_theta$df, as.numeric(df_SEM)),
  p_value   = c(test_theta$p.value, as.numeric(p_SEM))
)

write_csv(
  lr_out,
  paste0("models/LR_tests_SDM_", run_id, "_", time_tag, ".csv")
)

# --------------------------------------------------------
# 4) log-likelihood 汇总
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
# 5) 模型比较表（logLik / AIC / BIC）
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

cat("\n✨ Step2 完成！所有结果已带 run_id + 时间戳保存 ✨\n")
