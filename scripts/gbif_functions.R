
# Countries test ----------------------------------------------------------


cc_coun_2 <- function (x, lon = "decimalLongitude", lat = "decimalLatitude", 
          iso3 = "countrycode", value = "clean", ref = NULL, ref_col = "iso_a3", 
          verbose = TRUE, buffer = NULL) 
{
  match.arg(value, choices = c("clean", "flagged"))
  if (!iso3 %in% names(x)) {
    stop("iso3 argument missing, please specify")
  }
  if (verbose) {
    message("Testing country identity")
  }
  if (is.null(ref)) {
    if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
      stop("Install the 'rnaturalearth' package or provide a custom reference", 
           call. = FALSE)
    }
    ref <- terra::vect(rnaturalearth::ne_countries(scale = "medium", 
                                                   returnclass = "sf"))
  }
  else {
    if (any(is(ref) == "Spatial") | inherits(ref, "sf")) {
      ref <- terra::vect(ref)
    }
    if (!(inherits(ref, "SpatVector") & terra::geomtype(ref) == 
          "polygons")) {
      stop("ref must be a SpatVector with geomtype 'polygons'")
    }
    ref <- reproj(ref)
  }
  dat <- terra::vect(x[, c(lon, lat)], geom = c(lon, lat), 
                     crs = "+proj=longlat +datum=WGS84")
  if (is.numeric(buffer)) {
    buffer <- ifelse(buffer == 0, 1e-11, buffer)
    ref_buff <- terra::buffer(ref, buffer)
    ref <- terra::vect(stats::na.omit(terra::geom(ref_buff)), 
                       type = "polygon", crs = "+proj=longlat +datum=WGS84")
    terra::values(ref) <- terra::values(ref_buff)
  }
  country <- terra::extract(ref, dat)
  count_dat <- as.character(unlist(x[, iso3]))
  if (is.numeric(buffer)) {
    out <- logical(length(dat))
    for (i in seq_along(dat)) {
      out[i] <- count_dat[i] %in% country[country[, 1] == 
                                            i, ref_col]
    }
  }
  else {
    country <- country[, ref_col]
    out <- as.character(country) == count_dat
    out[is.na(out)] <- FALSE
  }
  if (verbose) {
    if (value == "clean") {
      message(sprintf("Removed %s records.", sum(!out)))
    }
    else {
      message(sprintf("Flagged %s records.", sum(!out)))
    }
  }
  switch(value, clean = return(x[out, ]), flagged = return(out))
}


# Zeros test --------------------------------------------------------------


cc_zero_2 <- function (x, lon = "decimalLongitude", lat = "decimalLatitude", 
          buffer = 0.5, value = "clean", verbose = TRUE) 
{
  match.arg(value, choices = c("clean", "flagged"))
  if (verbose) {
    message("Testing zero coordinates")
  }
  t1 <- !(x[[lon]] == 0 | x[[lat]] == 0)
  dat <- terra::vect(x[, c(lon, lat)], geom = c(lon, lat))
  if (buffer == 0) {
    buffer <- 1e-14
  }
  buff <- terra::buffer(terra::vect(data.frame(lat = 0, lon = 0), crs = "+proj=longlat +datum=WGS84"), 
                        width = buffer)
  ext_dat <- terra::extract(buff, dat)
  t2 <- is.na(ext_dat[!duplicated(ext_dat[, 1]), 2])
  out <- Reduce("&", list(t1, t2))
  if (verbose) {
    if (value == "clean") {
      message(sprintf("Removed %s records.", sum(!out)))
    }
    else {
      message(sprintf("Flagged %s records.", sum(!out)))
    }
  }
  switch(value, clean = return(x[out, ]), flagged = return(out))
}
