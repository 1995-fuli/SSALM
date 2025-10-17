suppressPackageStartupMessages({
  library(spdep); library(spatialreg); library(Matrix)
})

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
