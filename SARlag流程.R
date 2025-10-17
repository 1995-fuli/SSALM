library(terra)

TAC <- rast('F:/02_SIF/Zhang_SIF/03_zhang_0.05/inter_0.25/SIF_01_20_lag1_60m_slope.nc')


# ==== 热浪：选择变量并用“年份名”筛选 ====
HW <- rast("F:/05_heatwave/MHW_features.nc")

# 方式1：用变量名挑选（若 names(HW) 里带有 'event_count'）
HW_cnt <- subset(HW, grep("event_count", names(HW)))
crs(HW_cnt) <- "EPSG:4326"
#plot(HW_cnt)

HW_cum <- subset(HW, grep("Total_cumulative_intensity", names(HW)))
crs(HW_cum) <- "EPSG:4326"
#plot(HW_cum)

HW_dur = rast("F:/05_heatwave/04Duration_nan.nc")
#plot(HW_dur)
crs(HW_dur) <- "EPSG:4326"
HW_dur <- flip(HW_dur, "vertical")    # 若原始上下颠倒

HW_max <-rast("F:/05_heatwave/05Maximum_intensity_nan.nc")
#plot(HW_max)
crs(HW_max) <- "EPSG:4326"
HW_max <- flip(HW_max, "vertical")

years <- 2001:2020
names(HW_cnt) <- paste0("cnt_", years)
names(HW_dur) <- paste0("dur_", years)
names(HW_cum) <- paste0("cum_", years)
names(HW_max) <- paste0("max_", years)

idx <- which(years >= 2006 & years <= 2020)  # 2006–2020
# 1) 次数（年均）
HWcnt_mean <- mean(HW_cnt[[idx]], na.rm=TRUE)

# 2) 累计强度（年均；若想表征总暴露，可用 sum(...)）
HWcum_mean <- mean(HW_cum[[idx]], na.rm=TRUE)

# 3) 年最大强度（多年平均的“年最大”）
HWmax_mean <- mean(HW_max[[idx]], na.rm=TRUE)

# 4) 持续时间（事件数加权的多年平均）
HWdur_mean <- mean(HW_dur[[idx]], na.rm=TRUE)


# ========= L1：最小增量 SAR-lag（SSALM） =========
# 依赖
suppressPackageStartupMessages({
  library(spdep)
  library(spatialreg)
  library(car)
})

# 将热浪四个聚合指标重采样到 TAC 网格（双线性即可；count 可用邻近，但差异很小）
HWcnt_mean_rs <- resample(HWcnt_mean, TAC, method = "near")
HWcum_mean_rs <- resample(HWcum_mean, TAC, method = "bilinear")
HWmax_mean_rs <- resample(HWmax_mean, TAC, method = "bilinear")
HWdur_mean_rs <- resample(HWdur_mean, TAC, method = "bilinear")

# 2) 统一掩膜：只保留 TAC 与热浪都非 NA 的像元
Xstack <- c(HWcnt_mean_rs, HWcum_mean_rs, HWmax_mean_rs, HWdur_mean_rs)
names(Xstack) <- c("cnt_mean","cum_mean","max_mean","dur_mean")
Xstack <- mask(Xstack, TAC)             # 先用 TAC 掩膜
okmask  <- app(c(TAC, Xstack), fun = function(x) as.integer(all(is.finite(x))))
Xstack  <- mask(Xstack, okmask, maskvalues = 0)
TAC_ok  <- mask(TAC,   okmask, maskvalues = 0)


# 3) 抽到点表（经纬度 + 响应 + 自变量），并可选抽样（大数据建议先抽 2 万像元调通）
dat <- as.data.frame(c(TAC_ok, Xstack), xy = TRUE, na.rm = TRUE)
names(dat)[3] <- "tac_slope"
head(dat)

# 可选：抽样以加速（先调通再全量/分块）
set.seed(1)
n_max <- 20000L
if (nrow(dat) > n_max) dat <- dat[sample(nrow(dat), n_max), ]

# 4) 标准化自变量
zcols <- c("cnt_mean","cum_mean","max_mean","dur_mean")
dat[ zcols ] <- scale(dat[ zcols ])
dat <- na.omit(dat)

# 5) 空间权重：kNN=8（全局半干旱区更稳健；行标准化）
coords <- as.matrix(dat[, c("x","y")])     # 经/纬度（WGS84）
knn    <- knearneigh(coords, k = 8, longlat = TRUE)
nb     <- knn2nb(knn, sym = TRUE)  # 关键：对称！
lw     <- nb2listw(nb, style = "W", zero.policy = TRUE)

# 6) OLS + 空间诊断
f   <- tac_slope ~ cnt_mean + cum_mean + max_mean + dur_mean
ols <- lm(f, data = dat)
cat("\n== OLS 概要 ==\n"); print(summary(ols))
cat("\nVIF:\n")
print(car::vif(ols))

cat("\n== Moran's I of OLS residuals ==\n")
print(lm.morantest(ols, lw, zero.policy = TRUE))
cat("\n== Lagrange Multiplier tests ==\n")
print(lm.LMtests(ols, lw, test = "all", zero.policy = TRUE))

# 7) SAR-lag 与 SEM（对照）并比较 AIC
sar <- lagsarlm(f, data = dat, listw = lw, method = "Matrix", zero.policy = TRUE)
sem <- errorsarlm(f, data = dat, listw = lw, method = "Matrix", zero.policy = TRUE)

cat("\n== AIC 对比（越小越好）==\n")
print(AIC(ols, sar, sem))

cat("\n== SAR-lag 概要（ρ 与热浪系数是关键）==\n")
print(summary(sar))

# 8) 可选：导出系数表，方便放文中
coef_tab <- data.frame(
  term     = c(names(coef(sar)), "rho"),
  estimate = c(unname(coef(sar)), sar$rho))
coef_tab

#write.csv(coef_tab, "SARlag_coefficients_L1.csv", row.names = FALSE)


##cum_mean与max_mean 高度共线
f_drop <- tac_slope ~ cnt_mean + max_mean + dur_mean
sar_drop <- lagsarlm(f_drop, data=dat, listw=lw, method="Matrix", zero.policy=TRUE)
summary(sar_drop) 
AIC(sar_drop)

dat$dur_resid <- resid(lm(dur_mean ~ cnt_mean + max_mean, data=dat))
sar_orth <- lagsarlm(tac_slope ~ cnt_mean + max_mean + dur_resid, data=dat, listw=lw, method="Matrix", zero.policy=TRUE)

summary(sar_orth)

# 和主模型对比 AIC（越小越好）
AIC(sar_drop, sar_orth)

# 基线（已跑）：去掉 cum_mean
sar_drop <- lagsarlm(tac_slope ~ cnt_mean + max_mean + dur_mean, data=dat, listw=lw, method="Matrix", zero.policy=TRUE)
# 进一步简化 1：去掉 dur_mean
sar_drop_nodur <- lagsarlm(tac_slope ~ cnt_mean + max_mean, data=dat, listw=lw, method="Matrix", zero.policy=TRUE)

# 进一步简化 2：只留 max_mean
sar_onlymax <- lagsarlm(tac_slope ~ max_mean, data=dat, listw=lw,
                        method="Matrix", zero.policy=TRUE)

AIC(sar_drop, sar_drop_nodur, sar_onlymax)

suppressPackageStartupMessages({
  library(spatialreg)
  library(spdep)
  library(car)
})
# dat: 包含列 tac_slope, cnt_mean, max_mean, dur_mean（自变量已标准化）
# lw : kNN=8 的行标准化 listw（W），zero.policy=TRUE

####1) 拟合最终模型 + 对照 OLS/SEM，并导出 AIC
# OLS（对照）
ols_keep <- lm(tac_slope ~ cnt_mean + max_mean + dur_mean, data = dat)

# SAR-lag（最终模型）
sar_keep <- lagsarlm(tac_slope ~ cnt_mean + max_mean + dur_mean,
                     data = dat, listw = lw, method = "Matrix", zero.policy = TRUE)

# SEM（对照）
sem_keep <- errorsarlm(tac_slope ~ cnt_mean + max_mean + dur_mean,
                       data = dat, listw = lw, method = "Matrix", zero.policy = TRUE)

cat("\n== AIC 对比（越小越好）==\n")
print(AIC(ols_keep, sar_keep, sem_keep))

# 保存 AIC 表
aic_tab <- AIC(ols_keep, sar_keep, sem_keep)
#write.csv(aic_tab, "AIC_comparison.csv", row.names = TRUE)
aic_tab

####2) 提取ρ（rho）点估与标准误 + 系数表（含 rho）
s_sum <- summary(sar_keep)  # 打印也很重要：print(s_sum)

# ρ 及其 SE（不同版本对象名略有差异，做个兜底）
rho_est <- if (!is.null(s_sum$rho)) s_sum$rho else sar_keep$rho
rho_se  <- if (!is.null(s_sum$rho.se)) s_sum$rho.se else {
  NA_real_
}

# 回归系数（不含 rho）的估计与 SE
coef_est <- s_sum$Coef[, "Estimate"]
coef_se  <- s_sum$Coef[, "Std. Error"]

# 组合成一张干净的表并保存
coef_tab <- data.frame(
  term     = c("rho", names(coef_est)),
  estimate = c(rho_est, unname(coef_est)),
  std_error = c(rho_se, unname(coef_se)),
  row.names = NULL)
#write.csv(coef_tab, "SARlag_coefficients.csv", row.names = FALSE)
coef_tab


## 0) 前提检查
stopifnot(all(c("x","y","tac_slope","cnt_mean","max_mean","dur_mean") %in% names(dat)))

## 小工具：兼容不同包版本的类名/字段名
is_sarlm <- function(x) inherits(x, c("sarlm","Sarlm"))

get_coef_table <- function(s_sum) {
  if (!is.null(s_sum$Coef)) return(s_sum$Coef)
  if (!is.null(s_sum$coefficients)) return(s_sum$coefficients)
  stop("未找到系数字段（既无 $Coef 也无 $coefficients）")
}

## 1) 重建 listw（kNN=8，行标准化，sym=TRUE）
coords <- as.matrix(dat[, c("x","y")])
nb     <- knn2nb(knearneigh(coords, k = 8, longlat = TRUE), sym = TRUE)
lw     <- nb2listw(nb, style = "W", zero.policy = TRUE)
stopifnot(length(lw$neighbours) == nrow(dat))

## 2) 拟合 SAR-lag（显式命名空间，避免同名覆盖）
sar_model <- spatialreg::lagsarlm(tac_slope ~ cnt_mean + max_mean + dur_mean,
                                  data = dat, listw = lw,
                                  method = "Matrix", zero.policy = TRUE)

cat("sar_model class: ", paste(class(sar_model), collapse=", "), "\n")
stopifnot(is_sarlm(sar_model))  # 大小写都接受

## 3) 用迹近似计算 Direct / Indirect / Total（点估；不做MC）
W_sparse <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
MOMENTS  <- 12L                                     # 10–15 足够；越小越快
TR_cache <- trW(W_sparse, type = "MC", moments = MOMENTS)

fast_impacts <- function(model, tr, W) {
  stopifnot(is_sarlm(model))
  n   <- nrow(W)
  ss  <- summary(model)
  rho <- if (!is.null(ss$rho)) ss$rho else model$rho
  
  Coef <- get_coef_table(ss)              # 兼容不同字段名
  beta <- Coef[, "Estimate"]
  beta <- beta[ setdiff(names(beta), "(Intercept)") ]  # 去截距
  vars <- names(beta)
  
  k_seq <- seq_along(tr)
  dirM  <- 1 + sum((rho^k_seq) * (tr / n))  # mean(diag((I - rho W)^-1))
  totM  <- 1 / (1 - rho)                    # 对 style="W" 成立
  indM  <- totM - dirM
  
  data.frame(
    Variable = vars,
    Direct   = as.numeric(beta) * dirM,
    Indirect = as.numeric(beta) * indM,
    Total    = as.numeric(beta) * totM,
    row.names = NULL)
}

imp_out <- fast_impacts(sar_model, tr = TR_cache, W = W_sparse)
print(imp_out)
# write.csv(imp_out, "SARlag_impacts_point.csv", row.names = FALSE)



#####方便投稿的小表（补充“占比”和 multiplier）
# 已有 imp_out: Variable, Direct, Indirect, Total
imp_out$Direct_share   <- imp_out$Direct   / imp_out$Total
imp_out$Indirect_share <- imp_out$Indirect / imp_out$Total

# 从已知关系反解 multiplier（所有变量相同）
dirM <- unique(round(imp_out$Direct / c(1,1,1) / 
                       c(NA, NA, NA), 6)) # 这行仅占位，下面用ρ更稳
# 更稳：直接用 rho 推导
ss  <- summary(sar_keep)  # 或 sar_model
rho <- if (!is.null(ss$rho)) ss$rho else sar_keep$rho
totM <- 1/(1-rho)
# 用 Total = beta*totM、Direct = beta*dirM ⇒ dirM = Direct/Total * totM
# 对每个变量都一样，所以取一行即可
dirM_est <- (imp_out$Direct[1] / imp_out$Total[1]) * totM
indM_est <- totM - dirM_est

mult_row <- data.frame(Multiplier = c("Direct_M","Indirect_M","Total_M"),
                       Value = c(dirM_est, indM_est, totM))
print(imp_out)
print(mult_row)
# write.csv(imp_out, "SARlag_impacts_point_with_share.csv", row.names = FALSE)

#快速补一份“可选 SE”的小代码（两选一）
#A. 子样本 MC（推荐，快）
set.seed(1)
idx <- sample(nrow(dat), 8000)
dat_sub <- dat[idx, ]
coords_sub <- as.matrix(dat_sub[, c("x","y")])
nb_sub  <- knn2nb(knearneigh(coords_sub, k=8, longlat=TRUE), sym=TRUE)
lw_sub  <- nb2listw(nb_sub, style="W", zero.policy=TRUE)

fit_sub <- spatialreg::lagsarlm(tac_slope ~ cnt_mean + max_mean + dur_mean,
        data=dat_sub, listw=lw_sub, method="Matrix", zero.policy=TRUE)

W_sub <- as(as_dgRMatrix_listw(lw_sub), "CsparseMatrix")
tr_sub <- trW(W_sub, type="MC", moments=12)

# 小R的MC
imp_sub <- impacts(fit_sub, listw=lw_sub, tr=tr_sub, R=40)
sum_sub <- summary(imp_sub, zstats=TRUE, short=TRUE)
as.data.frame(sum_sub$res)  # 给出均值/SE（Direct/Indirect/Total）

##全量小R
W_full <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
tr_full <- trW(W_full, type="MC", moments=12)
set.seed(1)
imp_quick <- impacts(sar_keep, listw=lw, tr=tr_full, R=20)  # R=20 或 40
summary(imp_quick, zstats=TRUE, short=TRUE)


####“闭环验证”
e_sar <- residuals(sar_keep, type = "response")
moran.test(e_sar, lw, alternative = "greater", zero.policy = TRUE)

suppressPackageStartupMessages({
  library(spdep); library(spatialreg); library(Matrix)
})


########“热浪指标 → TAC slope”的 SAR-lag（SSALM） 结果汇总
## ============ 前提：lw、sar_keep、ols_keep、sem_keep 已按 k=8 拟合 ============
## ---------- 表A：AIC + ρ±SE + 系数±SE ----------
# AIC 对比
tabA_aic <- AIC(ols_keep, sar_keep, sem_keep)

# ρ 与系数表
s_sum   <- summary(sar_keep)
rho_est <- if (!is.null(s_sum$rho)) s_sum$rho else sar_keep$rho
rho_se  <- if (!is.null(s_sum$rho.se)) s_sum$rho.se else NA_real_

coef_mat <- s_sum$Coef
tabA_coef <- data.frame(
  term      = c("rho", rownames(coef_mat)),
  estimate  = c(rho_est,  coef_mat[, "Estimate"]),
  std_error = c(rho_se,   coef_mat[, "Std. Error"]),
  row.names = NULL)

# 保存（可选）
# write.csv(tabA_aic,  "TableA_AIC.csv")
# write.csv(tabA_coef, "TableA_SARlag_coefficients.csv")

tabA_aic
tabA_coef

## ---------- 表B：影响分解（全量小R） ----------
# 用稀疏矩阵 + 迹近似准备 tr，moments=12（与正文一致）
W_full  <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
tr_full <- trW(W_full, type="MC", moments=12)

set.seed(1)
imp_quick <- impacts(sar_keep, listw=lw, tr=tr_full, R=20)  # “全量小R”
sum_quick <- summary(imp_quick, zstats=TRUE, short=TRUE)

# 点估（Direct/Indirect/Total）
means <- as.data.frame(sum_quick$res)    # 列名: direct/indirect/total
# 统一列名（有些版本用大写首字母）
names(means) <- tolower(names(means))  # 期望得到 direct / indirect / total

#标准误：优先用sres；没有则用z反推（SE = |mean / z|）；都没有则NA
if (!is.null(sum_quick$sres)) {
  se_df <- as.data.frame(sum_quick$sres)
  names(se_df) <- tolower(names(se_df))
} else if (!is.null(sum_quick$zmat)) {
  z_df <- as.data.frame(sum_quick$zmat)
  names(z_df) <- tolower(names(z_df))
  # 防 0：若有 0，用 NA 避免除零
  z_mat <- as.matrix(z_df); z_mat[abs(z_mat) < .Machine$double.eps] <- NA_real_
  se_df <- as.data.frame(abs(as.matrix(means) / z_mat))
  names(se_df) <- paste0("se_", names(means))
} else {
  se_df <- data.frame(direct=NA_real_, indirect=NA_real_, total=NA_real_)
}

# 3) z / p：若 summary 里给了就用；否则用 mean/SE 自行计算
if (!is.null(sum_quick$zmat)) {
  z_df <- as.data.frame(sum_quick$zmat); names(z_df) <- paste0("z_", tolower(names(z_df)))
} else {
  # 用上一步得到的 se_df 反算
  tmp_se <- se_df
  names(tmp_se) <- gsub("^se_", "", names(tmp_se))  # 统一成 direct/indirect/total
  z_df <- as.data.frame(as.matrix(means) / as.matrix(tmp_se))
  names(z_df) <- paste0("z_", names(z_df))
}

if (!is.null(sum_quick$pzmat)) {
  p_df <- as.data.frame(sum_quick$pzmat); names(p_df) <- paste0("p_", tolower(names(p_df)))
} else {
  # 正态近似反算 p
  p_df <- as.data.frame(2 * pnorm(abs(as.matrix(z_df)), lower.tail = FALSE))
  names(p_df) <- gsub("^z_", "p_", names(p_df))
}

# 4) 组装表B
tabB <- data.frame(
  Variable = rownames(means),
  Direct   = means$direct,
  Indirect = means$indirect,
  Total    = means$total,
  row.names = NULL)

# 合并 SE / z / p
# 若 se_df 列名为 direct/indirect/total 或 se_direct... 都兼容处理
names(se_df) <- sub("^se_", "", names(se_df))
tabB$SE_Direct   <- se_df$direct
tabB$SE_Indirect <- se_df$indirect
tabB$SE_Total    <- se_df$total

tabB$z_Direct    <- z_df$z_direct
tabB$z_Indirect  <- z_df$z_indirect
tabB$z_Total     <- z_df$z_total

tabB$p_Direct    <- p_df$p_direct
tabB$p_Indirect  <- p_df$p_indirect
tabB$p_Total     <- p_df$p_total

# 占比
tabB$Direct_share   <- tabB$Direct   / tabB$Total
tabB$Indirect_share <- tabB$Indirect / tabB$Total

# 如果行名不是变量名（比如 "1""2""3"），做个兜底按顺序映射
if (!all(tabB$Variable %in% c("cnt_mean","max_mean","dur_mean"))) {
  tabB$Variable <- c("cnt_mean","max_mean","dur_mean")
}

# 查看/保存
print(tabB)
# write.csv(tabB, "TableB_Impacts_full_R20.csv", row.names = FALSE)

