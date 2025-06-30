library(terra)
library(rgeomorphon)

SEARCH = 3       # outer search radius (cells)
SKIP = 0         # inner skip radius (cells)
DIST = 0         # flatness distance (cells)
FLAT = 1         # flat angle threshold
MODE = "anglev1" # comparison mode
INTERIOR = FALSE # check only interior of raster (avoid difference at edges)

## classic volcano
dem <- terra::rast(volcano)
terra::crs(dem) <- terra::crs("EPSG:2193")
terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)

## new volcano
# data(volcano2, package = "tidyterra")
# dem <- terra::rast(volcano2)
# terra::crs(dem) <- terra::crs("EPSG:2193")
# terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
# names(dem) <- "volcano"

# tahoe lidar
# dem <- rast(system.file("extdata", "tahoe_lidar_bareearth.tif", package="gdalUtilities"))
# dem <- project(dem, "EPSG:5070")
# dem_t <- rast(dem)
# res(dem_t) <- 1
# ext(dem_t) <- round(ext(dem_t))
# dem <- project(dem, dem_t)

system.time({
    rg <- geomorphons(
        dem,
        search = SEARCH,
        skip = SKIP,
        dist = DIST,
        flat = FLAT,
        comparison_mode = MODE
    )
})

system.time({
    tg <- geomorphons(
        dem,
        search = SEARCH,
        skip = SKIP,
        dist = DIST,
        flat = FLAT,
        comparison_mode = MODE,
        forms = FALSE,
        ternary = TRUE
    )
})
library(rgrass)
initGRASS(SG = dem, override = TRUE)
write_RAST(dem, "volcano2", flags=c("o", "overwrite"))
system.time({
    execGRASS(
        "r.geomorphon",
        elevation = "volcano2",
        search = SEARCH,
        skip = SKIP,
        dist = DIST,
        flat = FLAT,
        comparison = MODE,
        forms = "volcano2_forms",
        flags = "overwrite"
    )
})
grassg <- read_RAST("volcano2_forms")


system.time({
    execGRASS(
        "r.geomorphon",
        elevation = "volcano2",
        search = SEARCH,
        skip = SKIP,
        dist = DIST,
        flat = FLAT,
        comparison = MODE,
        ternary = "volcano2_ternary",
        flags = "overwrite"
    )
})
terng <- read_RAST("volcano2_ternary")

if (INTERIOR) {
    # extract cells that are at least SEARCH from border
    rgc <- geomorphon_theme(crop(rg, ext(rg) - SEARCH * res(rg)[1]))
    grassgc <- geomorphon_theme(crop(grassg, ext(grassg) - SEARCH * res(grassg)[1]))
} else {
    rgc <- geomorphon_theme(rg)
    grassgc <- geomorphon_theme(grassg)
}

c(rgc, grassgc) |> plot()
plot((rgc - grassgc) == 0)

vals_r <- values(rgc)[,1] |>
    factor(levels = 1:10, labels = geomorphon_categories()$form)
vals_grass <- values(grassgc)[,1] |>
    factor(levels = 1:10, labels = geomorphon_categories()$form)

comparison_df <- data.frame(
    r_val = vals_r,
    grass_val = vals_grass,
    match = ifelse(is.na(vals_r) & is.na(vals_grass), TRUE,
                   ifelse(is.na(vals_r) | is.na(vals_grass), FALSE,
                          vals_r == vals_grass))
)
message(paste("Number of matching cells:", sum(comparison_df$match, na.rm = TRUE)))
message(paste("Number of non-matching cells:", sum(!comparison_df$match, na.rm = TRUE)))

# Confusion matrix / crosstab
message("Cross-tabulation of forms:")
ctf <- table(R_Forms = vals_r, GRASS_Forms = vals_grass, useNA = "ifany")
print(ctf)
message("Accuracy: ", round(sum(diag(ctf)) / sum(ctf) * 100, 2), "%")

ctf <- table(R_Forms = vals_r, GRASS_Forms = vals_grass, useNA = "no")
print(ctf)
message("Accuracy (no-NA): ", round(sum(diag(ctf)) / sum(ctf) * 100, 2), "%")

# Identify specific differing cells (example: first 10)
differing_indices <- which(!comparison_df$match)
if (length(differing_indices) > 0) {
    message("Examples of differing cells (R vs GRASS):")
    head(comparison_df[differing_indices, ], 10)

    # You can convert differing_indices back to row/col if needed for debugging
    # coords_diff <- terra::xyFromCell(r_forms, differing_indices)
}

# ## test whiteboxtools geomorphon algorithm
library(whitebox)
writeRaster(dem, "test.tif", overwrite=TRUE)
wbt_geomorphons(
    "test.tif",
    "geomorphons.tif",
    search = SEARCH,
    skip = SKIP,
    threshold = FLAT,
    fdist = 0,
    forms = TRUE
)

wbtg <- rast("geomorphons.tif")
if (INTERIOR) {
    wbtgc <- crop(wbtg, ext(wbtg) - SEARCH * res(wbtg)[1])
} else {
    wbtgc <- wbtg
}
plot(geomorphon_theme(wbtgc))
plet((wbtgc - grassgc) == 0)
plet((wbtgc - rgc) == 0)

