library(rgeomorphon)
library(terra)
library(fields, warn.conflicts = FALSE)
# source: Bathymetric Contours (1 foot) - Salton Sea [ds426]
# https://map.dfg.ca.gov/metadata/ds0426.html

v <- vect("~/Geodata/BATHYMETRY/SaltonSea/ds426.shp")

plot(v, "Contour")
plot(as.points(v), "Contour", breaks = 20)
o0 <- as.polygons(v[v$Contour==-221,][1,])
o <- hull(o0, "concave_ratio")

r <- rast(round(ext(as.polygons(v, ext = TRUE)), -2), resolution = 100, crs = crs(v))
r

ri <- rasterize(as.points(v), r, "Contour")
plot(ri)

set.seed(123)
sv <- spatSample(ri, size=10000, method="regular", na.rm=TRUE, as.points=TRUE)
tps_model <- fields::Tps(crds(sv), sv$last)

r2 <- interpolate(r, tps_model) |>
    mask(o0)

plot(r2)
plet(r2)
contour(r2, filled=TRUE)

r3 <- aggregate(r2, 3)

plot(r3)
plet(r3)
contour(r3, filled=TRUE)

# inspect geomorphons output
plot(geomorphons(r3, flat_angle_deg = 0.1, search = 10, skip = 3))
plot(geomorphons(r3, flat_angle_deg = 0.1, search = 10, skip = 3, forms = "forms6"))
plot(geomorphons(r3, flat_angle_deg = 0.1, search = 10, skip = 3, forms = "forms5"))
plot(geomorphons(r3, flat_angle_deg = 0.1, search = 10, skip = 3, forms = "forms4"))
plot(geomorphons(-r3, flat_angle_deg = 0.1, search = 10, skip = 3))

# convert to matrix with elevation in feet
salton <- matrix(as.matrix(r3), nrow = nrow(r3), ncol = ncol(r3), byrow = TRUE)
attr(salton, "crs") <- terra::crs(r3)
attr(salton, "extent") <- terra::ext(r3)[1:4]

# inspect
str(salton)

usethis::use_data(salton, compress = "bzip2", overwrite = TRUE)

