# Apply Geomorphon Theme to Result Object

Applies standard class names and colors to a SpatRaster, or creates a
factor matrix. Input values should be integers between 1 and 10.

## Usage

``` r
geomorphon_categories()

geomorphon_colors()

geomorphon_theme(x, forms = "forms10")
```

## Arguments

- x:

  A SpatRaster or matrix object.

- forms:

  character. One of: `"forms10"` (default), `"forms6"`, `"forms5"`,
  `"forms4"`. These are themes corresponding to the built-in 10-form,
  6-form, 5-form, and 4-form `"forms`" outputs from
  [`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md).

## Value

A SpatRaster or matrix object with geomorphon class names (and colors
for SpatRaster) applied.

## Details

When `x` is a matrix the result is a factor using
`geomorphon_categories()`. Values are integers 1 to 10 and labels are
the geomorphon form names.

## Examples

``` r
geomorphon_theme(1:10)
#>  [1] "Flat"      "Peak"      "Ridge"     "Shoulder"  "Spur"      "Slope"    
#>  [7] "Hollow"    "Footslope" "Valley"    "Pit"      
```
