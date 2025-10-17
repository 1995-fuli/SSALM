library(terra)

##SPEI-03 → 干/湿事件特征
spei03 <- rast("F:/revised/SARlag/spei03.nc")
idx <- which(format(time(spei03), "%Y-%m") >= "2001-01" &
               format(time(spei03), "%Y-%m") <= "2020-12")
spei03_01_20 <- spei03[[idx]]   # 取你需要的2001–2020
nlyr(spei03_01_20)
plot(spei03_01_20)


slope <- rast('F:/02_SIF/Zhang_SIF/03_zhang_0.05/inter_0.25/SIF_01_20_lag1_60m_slope.nc')

tt2 <- time(spei03_01_20)
idx_05_20 <- which(format(tt2, "%Y-%m") >= "2005-12" & format(tt2, "%Y-%m") <= "2020-12")
spei03_05_20 <- spei03_01_20[[idx_05_20]]   # 16*12 = 192 层

spei_crop   <- crop(spei03_05_20, ext(slope))
spei_match  <- resample(spei_crop, slope, method = "bilinear")  # 0.25°
spei_match 

## 3) 干/湿月阈值与掩膜
thr_dry <- -1   # 干旱月阈值：SPEI < -1
thr_wet <-  1   # 湿润月阈值：SPEI > +1

# 干旱月掩膜：满足条件=1，否则=NA
dry_mask <- classify(spei_match, rcl = matrix(c(-Inf, thr_dry, 1,
                                                thr_dry, Inf,  NA), ncol=3, byrow=TRUE))
# 湿润月掩膜
wet_mask <- classify(spei_match, rcl = matrix(c(-Inf, thr_wet,  NA,
                                                thr_wet, Inf, 1), ncol=3, byrow=TRUE))

## 4) 单像元事件特征函数（向量化）
# 传入：c(SPEI序列, 掩膜序列)；返回：
# 干旱：事件次数、平均持续（月）、平均强度（各段最小SPEI的绝对值均值）
get_dry_stats_vec <- function(v){
  n <- length(v) / 2
  x <- v[1:n]
  m <- v[(n+1):(2*n)]
  if (all(is.na(m))) return(c(0, NA, NA))
  idx <- which(!is.na(m))
  runs <- split(idx, cumsum(c(1, diff(idx) != 1)))   # 连续月份分段
  n_evt <- length(runs)
  dur   <- sapply(runs, length)
  x_evt <- x[idx]
  seg   <- split(x_evt, rep(seq_along(runs), dur))
  inten <- sapply(seg, function(z) abs(min(z, na.rm=TRUE)))
  c(n_evt,
    mean(dur, na.rm=TRUE),
    mean(inten, na.rm=TRUE))
}

# 湿润：总严重度 = 各湿润段内 SPEI（正值）求和的绝对值，再对所有段求和
get_wet_totsev_vec <- function(v){
  n <- length(v) / 2
  x <- v[1:n]
  m <- v[(n+1):(2*n)]
  if (all(is.na(m))) return(NA)
  idx <- which(!is.na(m))
  runs <- split(idx, cumsum(c(1, diff(idx) != 1)))
  dur   <- sapply(runs, length)
  x_evt <- x[idx]
  seg   <- split(x_evt, rep(seq_along(runs), dur))
  sev   <- sapply(seg, function(z) abs(sum(z, na.rm=TRUE)))
  sum(sev, na.rm=TRUE)
}

## 5) 逐像元计算四个指标
dry_input <- c(spei_match, dry_mask)
wet_input <- c(spei_match, wet_mask)

dry_stats <- app(dry_input, get_dry_stats_vec)     # 三层
names(dry_stats) <- c("D_freq","D_avdur","D_avint")

wet_totsev <- app(wet_input, get_wet_totsev_vec)   # 一层
names(wet_totsev) <- "Wet_totsev"

spei_features <- c(dry_stats, wet_totsev)          # 合并四层
spei_features <- mask(spei_features, slope)        # 用 slope 的有效范围做最终掩膜（-60~60）

## 3) 干/湿月阈值与掩膜
thr_dry <- -1   # 干旱月阈值：SPEI < -1
thr_wet <-  1   # 湿润月阈值：SPEI > +1

# 干旱月掩膜：满足条件=1，否则=NA
dry_mask <- classify(spei_match, rcl = matrix(c(-Inf, thr_dry, 1,
                                                 thr_dry, Inf,  NA), ncol=3, byrow=TRUE))
# 湿润月掩膜
wet_mask <- classify(spei_match, rcl = matrix(c(-Inf, thr_wet,  NA,
                                                  thr_wet, Inf, 1), ncol=3, byrow=TRUE))

## 4) 单像元事件特征函数（向量化）
# 传入：c(SPEI序列, 掩膜序列)；返回：
# 干旱：事件次数、平均持续（月）、平均强度（各段最小SPEI的绝对值均值）
get_dry_stats_vec <- function(v){
  n <- length(v) / 2
  x <- v[1:n]
  m <- v[(n+1):(2*n)]
  if (all(is.na(m))) return(c(0, NA, NA))
  idx <- which(!is.na(m))
  runs <- split(idx, cumsum(c(1, diff(idx) != 1)))   # 连续月份分段
  n_evt <- length(runs)
  dur   <- sapply(runs, length)
  x_evt <- x[idx]
  seg   <- split(x_evt, rep(seq_along(runs), dur))
  inten <- sapply(seg, function(z) abs(min(z, na.rm=TRUE)))
  c(n_evt,
    mean(dur, na.rm=TRUE),
    mean(inten, na.rm=TRUE))
}

# 湿润：总严重度 = 各湿润段内 SPEI（正值）求和的绝对值，再对所有段求和
get_wet_totsev_vec <- function(v){
  n <- length(v) / 2
  x <- v[1:n]
  m <- v[(n+1):(2*n)]
  if (all(is.na(m))) return(NA)
  idx <- which(!is.na(m))
  runs <- split(idx, cumsum(c(1, diff(idx) != 1)))
  dur   <- sapply(runs, length)
  x_evt <- x[idx]
  seg   <- split(x_evt, rep(seq_along(runs), dur))
  sev   <- sapply(seg, function(z) abs(sum(z, na.rm=TRUE)))
  sum(sev, na.rm=TRUE)
}

## 5) 逐像元计算四个指标
dry_input <- c(spei_match, dry_mask)
wet_input <- c(spei_match, wet_mask)

dry_stats <- app(dry_input, get_dry_stats_vec)     # 三层
names(dry_stats) <- c("D_freq","D_avdur","D_avint")

wet_totsev <- app(wet_input, get_wet_totsev_vec)   # 一层
names(wet_totsev) <- "Wet_totsev"

spei_features <- c(dry_stats, wet_totsev)          # 合并四层
spei_features <- mask(spei_features, slope)        # 用 slope 的有效范围做最终掩膜（-60~60）

plot(spei_features, nc=2)
global(spei_features, mean, na.rm=TRUE)

writeCDF(spei_features[[1]], filename = "F:/revised/SARlag/D_freq.nc", varname = "D_freq_count")
writeCDF(spei_features[[2]], filename = "F:/revised/SARlag/D_avdur.nc", varname = "D_avdur_month")
writeCDF(spei_features[[3]], filename = "F:/revised/SARlag/D_avint.nc", varname = "D_avint")
writeCDF(spei_features[[4]], filename = "F:/revised/SARlag/Wet_totsev.nc", varname = "Wet_totsev")



