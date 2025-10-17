library(terra)

pre = rast("F:/ERA5_land/total precipitation/ERA5 tp-LMAD2001-2020_arid_month_float0.25.nc")
pre

# 1) 单位换算：m -> mm
pre_mm <- pre * 1000
names(pre_mm) <- paste0("mm_", names(pre))

# 2) 取出年月索引
tt   <- time(pre_mm)
yyyy <- as.integer(format(tt, "%Y"))
mmmm <- as.integer(format(tt, "%m"))


# A) 年内季节性（基于“月气候”12层的CV）
# 步骤：对每个历月求20年均值 -> 得到12层气候 -> 在这12层上算CV
# -----------------------------
# 2.1 月气候：tapp 按月份聚合
clim12 <- tapp(pre_mm, index = mmmm, fun = mean, na.rm = TRUE)  # 12层，Jan..Dec
names(clim12) <- month.abb

# 2.2 在12层上计算 CV = sd/mean（对极干像元设置微小epsilon）
eps <- 1e-6
prec_seasonality_cv <- app(clim12, function(v){
  mu <- mean(v, na.rm = TRUE)
  sdv <- sd(v, na.rm = TRUE)
  if (is.na(mu)) return(NA)
  sdv / max(mu, eps)
})
names(prec_seasonality_cv) <- "prec_seasonality_cv"

# -----------------------------
# B) 年际变率（基于“年总降水”的CV）
# 步骤：按年求和 -> 20层年值 -> 在这20层上算CV
# -----------------------------
# 2.3 年总降水
# 生成一个 1..20 的年索引
yr_index <- match(yyyy, sort(unique(yyyy)))   # 2001..2020 -> 1..20
annual_sum <- tapp(pre_mm, index = yr_index, fun = sum, na.rm = TRUE)
names(annual_sum) <- paste0("P", sort(unique(yyyy)))  # P2001..P2020

# 2.4 年际 CV
prec_interannual_cv <- app(annual_sum, function(v){
  mu <- mean(v, na.rm = TRUE)
  sdv <- sd(v, na.rm = TRUE)
  if (is.na(mu)) return(NA)
  sdv / max(mu, eps)
})
names(prec_interannual_cv) <- "prec_interannual_cv"

# 3) 可选：把两个指标合并保存
prec_variability <- c(prec_seasonality_cv, prec_interannual_cv)
#writeRaster(prec_variability, "prec_variability_ERA5_2001_2020.tif", overwrite = TRUE)

# 4) 简单检查（推荐）
plot(prec_seasonality_cv, main="Precip seasonality CV (2001–2020)")
plot(prec_interannual_cv, main="Precip interannual CV (2001–2020)")
global(prec_variability, mean, na.rm=TRUE)


writeCDF(prec_seasonality_cv, filename = "F:/revised/SARlag/prec_seasonality_cv.nc", varname = "seasonality_cv")
writeCDF(prec_interannual_cv, filename = "F:/revised/SARlag/prec_interannual_cv.nc", varname = "interannual_cv")



