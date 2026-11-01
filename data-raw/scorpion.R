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

geomorphon_ptrm(g)

p <- geomorphon_ptrm_focal(g, 10)
