# Create a `forms_matrix` object

This constructor function wraps a 9x9 integer matrix and associates it
with a set of levels, creating a 'forms_matrix' object.

## Usage

``` r
forms_matrix(x, levels = get_forms_grass_enum())
```

## Arguments

- x:

  Integer. A 9x9 matrix.

- levels:

  Named integer vector. Map of integer values to their string names.
  Default: `get_forms_grass_enum()`

## Value

An object of class `c("forms_matrix", "matrix", "array")`.

## Details

This function is intended for custom classification matrix based on
positive and negative overlooks. See
[`forms_matrix_get()`](http://humus.rocks/rgeomorphon/reference/forms_matrix_get.md)
for a convenient accessor for the standard classification systems with
4, 5, 6 or 10 forms.

## Examples

``` r
library(terra)
#> terra 1.8.93
library(rgeomorphon)

# default values
x <- forms_matrix_get(num_forms = 10, levels = get_forms_grass_enum())

# inspect
x
#> <forms_matrix> object
#> 
#>       pos=0  pos=1  pos=2  pos=3  pos=4  pos=5  pos=6  pos=7  pos=8 
#> neg=0 G_FL   G_FL   G_FL   G_FS   G_FS   G_VL   G_VL   G_VL   G_PT  
#> neg=1 G_FL   G_FL   G_FS   G_FS   G_FS   G_VL   G_VL   G_VL   G_NONE
#> neg=2 G_FL   G_SH   G_SL   G_SL   G_HL   G_HL   G_VL   G_NONE G_NONE
#> neg=3 G_SH   G_SH   G_SL   G_SL   G_SL   G_HL   G_NONE G_NONE G_NONE
#> neg=4 G_SH   G_SH   G_SP   G_SL   G_SL   G_NONE G_NONE G_NONE G_NONE
#> neg=5 G_RI   G_RI   G_SP   G_SP   G_NONE G_NONE G_NONE G_NONE G_NONE
#> neg=6 G_RI   G_RI   G_RI   G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE
#> neg=7 G_RI   G_RI   G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE
#> neg=8 G_PK   G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE

# create a 9-class system where PEAK is combined with RIDGE
x[x == 2] <- 3
a <- get_forms_grass_enum()
a <- a[!names(a) == "G_PK"]

# create a forms matrix with custom levels
fm <- forms_matrix(x, a)

# run geomorphon algorithm
SEARCH = 7       # outer search radius (cells)
SKIP = 1         # inner skip radius (cells)
DIST = 0         # flatness distance (cells)
FLAT = 1         # flat angle threshold
MODE = "anglev1" # comparison mode

## classic volcano
data("volcano", package = "datasets")
dem <- terra::rast(volcano)
terra::crs(dem) <- terra::crs("EPSG:2193")
terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
names(dem) <- "elevation"

# include original forms, positive, and negative output
res <- geomorphons(
    dem,
    search = SEARCH,
    skip = SKIP,
    dist = DIST,
    flat = FLAT,
    comparison_mode = MODE,
    forms = TRUE,
    positive = TRUE,
    negative = TRUE
)

 # apply custom classification to positive and negative
 res2 <- geomorphon_theme(
   forms_matrix_apply(
       x = res[[c("positive", "negative")]],
       rcl = fm
   )
 )

 # compare with default
 terra::plot(terra::rast(c(`10 form`=res$forms, `9 form`=res2)))
```
