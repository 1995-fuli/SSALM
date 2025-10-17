library(terra)
library(sf)
library(exactextractr)

tmpl <- rast('F:/02_SIF/Zhang_SIF/03_zhang_0.05/inter_0.25/SIF_01_20_lag1_60m_slope.nc')

# 1) 用“唯一ID”生成逐格多边形（关键：dissolve=FALSE）
id_r <- tmpl
values(id_r) <- 1:ncell(id_r)              # 给每个0.25°格子一个独立ID
grid025_ll <- as.polygons(id_r, dissolve = FALSE) |>
  st_as_sf() |>
  st_make_valid()
names(grid025_ll)[1] <- "cell_id"          # 第一列就是刚才的唯一ID

# 只保留tmpl的有效格（若tmpl有大量NA）
keep <- which(!is.na(values(tmpl)))
grid025_ll <- grid025_ll[grid025_ll$cell_id %in% keep, ]

dir_hfi <- "F:/revised/SARlag/hfp/16571064"  # hfpYYYY.tif 所在文件夹
years <- 2006:2020 

# --- 2) 获取 HFI 的 Mollweide 投影，并把网格面重投影过去 ---
hfi_sample <- rast(file.path(dir_hfi, sprintf("hfp%d.tif", years[1])))
crs_mw <- crs(hfi_sample)
grid025_mw <- st_transform(grid025_ll, crs_mw)

# ---3) 单年份：在 Mollweide 中做“覆盖面积加权均值”的聚合 ---
# 说明：exactextractr::exact_extract(..., fun='mean') 会按单元格与多边形的
# 覆盖比例做权重（在等积投影中即为面积权重），得到真实的面积平均值。
agg_one_year <- function(y){
  f <- file.path(dir_hfi, sprintf("hfp%d.tif", y))
  stopifnot(file.exists(f))
  r <- rast(f)   # Mollweide, ~1 km
  # 覆盖比例加权均值（在等积投影下 = 面积加权）
  vals <- exact_extract(r, grid025_mw, 'mean')
  grid025_mw$val <- vals
  # 投回 WGS84，并栅格化到 0.25° 模板（按 cell_id 放回去）
  grid_ll <- st_transform(grid025_mw, crs(tmpl))
  rasterize(vect(grid_ll), tmpl, field = "val")
}

# --- 4) 批量聚合为年度 0.25° 栈 ---
cat("Aggregating HFI to 0.25° (area-weighted)...\n")

hfi_list <- lapply(years, agg_one_year)
hfi_stack_025 <- rast(hfi_list)
names(hfi_stack_025) <- paste0("HFI", years)
hfi_stack_025 <- mask(hfi_stack_025, tmpl)        # 只保留 tmpl 有效域
plot(hfi_stack_025)
hfi_stack_025
# --- 5) 主变量：2006–2020 多年均值 ---
hfi_mean_0620 <- app(hfi_stack_025, mean, na.rm=TRUE)

# --- 7) 快速检查 ---
print(hfi_mean_0620)
plot(hfi_mean_0620)
print(global(hfi_mean_0620, "mean", na.rm=TRUE))

writeCDF(hfi_stack_025, filename = "F:/revised/SARlag/hfp/hfi_stack_025_2006-2020.nc", varname = "hfp")
writeCDF(hfi_mean_0620, filename = "F:/revised/SARlag/hfp/hfi_mean_0620.nc", varname = "hfp")

