library(rgeomorphon)
library(soilDB)
library(terra)
library(rgeedim)

b <- terra::vect("POLYGON((
  -120.9613 38.1600,
  -120.9259 38.1600,
  -120.9259 38.1743,
  -120.9613 38.1743,
  -120.9613 38.1600))",
  crs = "OGC:CRS84"
) |>
  terra::project("EPSG:5070") |>
  terra::as.polygons(ext = TRUE)

x <- soilDB::SDA_spatialQuery(
  b,
  what = "mupolygon",
  geomIntersection = TRUE
) |>
  terra::project("EPSG:5070")


# search and download from USGS 1m lidar data collection
project_id <- Sys.getenv("GOOGLE_CLOUD_QUOTA_PROJECT", "rgeedim-demo")
gd_initialize(project = project_id)

r <- gd_region(x)

# note that resampling is done on the images as part of compositing/before download
y <- "USGS/3DEP/1m" |>
  gd_collection_from_name() |>
  gd_search(region = r) |>
  gd_composite(resampling = "bilinear") |>
  gd_download(
    region = r,
    bands = list("elevation"),
    crs = "EPSG:5070",
    scale = 1
  ) |>
  rast()

# inspect
plot(
  terra::shade(
    terra::terrain(y$elevation, unit = "radians"),
    terra::terrain(y$elevation, "aspect", unit = "radians")
  ),
  col = grey((0:100) / 100), legend = FALSE, mar = c(2, 2, 1, 4)
)

plot(project(x, y), add = TRUE)

y10 <- terra::aggregate(y, 10)

plot(y10)

g <- geomorphons(y10, search = 7, skip = 0, flat_angle_deg = 2)

plot(g)

r <- rasterize(x, y10, "mukey")

geomorphon_concordance(g, r)

geomorphon_ptrm(g, mode = "class")

p <- geomorphon_ptrm_focal(g, radius = 10, mode = "class")

plot(p)

# query mapunit and muaggatt attributes for categorical raster levels
mukeys <- sort(unique(x$mukey))
mu_tab <- soilDB::SDA_query(sprintf("
  SELECT
    mu.mukey,
    mu.musym,
    mu.muname,
    mu.mukind,
    mu.farmlndcl,
    ma.slopegraddcp,
    ma.slopegradwta,
    ma.drclassdcd,
    ma.hydgrpdcd,
    ma.aws0100wta
  FROM mapunit AS mu
  LEFT JOIN muaggatt AS ma ON mu.mukey = ma.mukey
  WHERE mu.mukey IN (%s)
  ORDER BY mu.mukey
", paste(mukeys, collapse = ",")))
attr(mu_tab, "SDA_id") <- NULL

# package 10m DEM and rasterized SSURGO mukey grid as matrices with CRS and extent attributes
dem_mat <- matrix(round(as.matrix(y10), 2), nrow = nrow(y10), ncol = ncol(y10), byrow = TRUE)
mukey_mat <- matrix(as.integer(as.matrix(r)), nrow = nrow(r), ncol = ncol(r), byrow = TRUE)

scorpion <- list(
  dem = dem_mat,
  mukey = mukey_mat,
  mapunit = mu_tab
)
attr(scorpion, "crs") <- terra::crs(y10)
attr(scorpion, "extent") <- terra::ext(y10)[1:4]

usethis::use_data(scorpion, compress = "bzip2", overwrite = TRUE)
