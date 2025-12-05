#Rscript ".\subset_nc_time_first10.r" ".\hunan.nc" ".\hunan10.nc" "tas" "time"

suppressPackageStartupMessages({
  library(ncdf4)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("用法: Rscript subset_nc_time_first10.r <hunan.nc> [out_nc] [var_name] [time_name]\n")
  cat("  说明: 读取NetCDF，保留time前10个切片，写出新文件。\n")
  quit(status = 1)
}

in_nc <- args[1]
out_nc <- ifelse(length(args) >= 2 && nzchar(args[2]), args[2], file.path(dirname(in_nc), "hunan10.nc"))
var_name_arg <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], "")
time_name_arg <- ifelse(length(args) >= 4 && nzchar(args[4]), args[4], "")

if (!file.exists(in_nc)) stop(paste("找不到NC文件:", in_nc))

nc <- ncdf4::nc_open(in_nc, readunlim = TRUE)
on.exit(ncdf4::nc_close(nc), add = TRUE)

# Pick variable with 3 dims if not provided
pick_var <- function(nc) {
  if (nzchar(var_name_arg) && var_name_arg %in% names(nc$var)) return(var_name_arg)
  # Prefer 'tas' if exists
  if ("tas" %in% names(nc$var)) return("tas")
  # Else first var with a time dimension
  for (vn in names(nc$var)) {
    vd <- nc$var[[vn]]$dim
    if (length(vd) >= 3) return(vn)
  }
  # Fallback to first var
  names(nc$var)[1]
}
var_name <- pick_var(nc)
var <- nc$var[[var_name]]

# Identify dims and time dim
get_time_name <- function(nc, var, time_name_arg) {
  if (nzchar(time_name_arg)) return(time_name_arg)
  # Common names
  cand <- c("time", "Time", "TIME")
  for (nm in cand) {
    if (!is.null(nc$dim[[nm]])) return(nm)
  }
  # Look in var dims
  for (d in var$dim) {
    if (grepl("time", d$name, ignore.case = TRUE)) return(d$name)
  }
  stop("未找到time维度名称，请提供 time_name")
}

time_name <- get_time_name(nc, var, time_name_arg)

# Read dims
dims <- var$dim
dim_names <- sapply(dims, function(d) d$name)
# Ensure order lon/lat/time or any order; we'll slice last dim that matches time_name
time_idx <- match(time_name, dim_names)
if (is.na(time_idx)) stop("变量中未包含time维度")

# Read full array then subset
arr <- ncdf4::ncvar_get(nc, var_name)
dim_arr <- dim(arr)
nt <- dim_arr[time_idx]
keep <- min(10, nt)

# Create start/count for slicing
start <- rep(1, length(dim_arr))
count <- dim_arr
count[time_idx] <- keep
arr10 <- ncdf4::ncvar_get(nc, var_name, start = start, count = count)

# Prepare output dims (copy original dims, shrink time)
make_dimdef <- function(d) {
  ncdf4::ncdim_def(name = d$name, units = d$units %||% "", vals = d$vals)
}
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

out_dims <- vector("list", length(dims))
for (i in seq_along(dims)) {
  di <- dims[[i]]
  vals <- di$vals
  if (i == time_idx) {
    vals <- vals[seq_len(keep)]
  }
  out_dims[[i]] <- ncdf4::ncdim_def(name = di$name, units = di$units %||% "", vals = vals)
}

# Define var (copy units and longname)
get_chr_attr <- function(nc, vname, aname, default = "") {
  at <- tryCatch(ncdf4::ncatt_get(nc, vname, aname), error = function(e) list(value = default))
  val <- at$value
  if (is.null(val)) return(default)
  # Ensure a single string
  if (is.raw(val)) return(default)
  if (is.list(val)) val <- tryCatch(val$value, error = function(e) default)
  if (length(val) == 0) return(default)
  as.character(val)[1]
}
var_units <- get_chr_attr(nc, var_name, "units", default = "")
var_long  <- get_chr_attr(nc, var_name, "long_name", default = var_name)
miss <- tryCatch(nc$var[[var_name]]$missval, error = function(e) NA_real_)
prec <- tryCatch(nc$var[[var_name]]$prec, error = function(e) "double")
if (is.null(prec) || !nzchar(prec)) prec <- "double"
var_def <- ncdf4::ncvar_def(name = var_name, units = var_units, dim = out_dims, missval = miss, longname = var_long, prec = prec)

# Create and write
nc_out <- ncdf4::nc_create(out_nc, vars = list(var_def))
# Copy global attrs
for (an in names(nc$gatts)) {
  ncdf4::ncatt_put(nc_out, 0, an, nc$gatts[[an]]$value)
}
# Write var data
ncdf4::ncvar_put(nc_out, var_def, arr10)
# Copy var attrs beyond units/long_name when present
for (an in names(nc$var[[var_name]]$att)) {
  if (an %in% c("units", "long_name")) next
  ncdf4::ncatt_put(nc_out, var_name, an, nc$var[[var_name]]$att[[an]]$value)
}

ncdf4::nc_close(nc_out)
message(paste("已写出:", out_nc, "变量:", var_name, "保留time前", keep, "个切片"))
