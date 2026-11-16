# Pairwise Dissimilarity of Concordance Profiles

Computes the pairwise dissimilarity matrix between conditional
probability profiles from a contingency table or
[`geomorphon_concordance()`](http://humus.rocks/rgeomorphon/reference/geomorphon_concordance.md)
output using exact base-2 Jensen-Shannon divergence (default) or any
method supported by `philentropy`.

## Usage

``` r
geomorphon_profile_dist(
  x,
  margin = 1,
  method = "jensen-shannon",
  unit = "log2"
)
```

## Arguments

- x:

  A `geomorphon_concordance` result list (containing
  `contingency_table`) or a 2D contingency `table` or numeric `matrix`.

- margin:

  Integer indicating the profile direction: `1` for row profiles (e.g.
  \\P(\text{Category} \mid \text{Geomorphon})\\), or `2` for column
  profiles (e.g. \\P(\text{Geomorphon} \mid \text{Category})\\). Default
  is `1`.

- method:

  Character string specifying distance metric. Default is
  `"jensen-shannon"`.

- unit:

  Logarithm unit (`"log2"`, `"log"`, `"log10"`). Default is `"log2"`.

## Value

A symmetric numeric distance matrix with row and column names matching
the profile labels.

## Examples

``` r
ct <- matrix(c(10, 0, 5, 5, 0, 10), nrow = 3, byrow = TRUE)
rownames(ct) <- c("Flat", "Slope", "Ridge")
colnames(ct) <- c("CatA", "CatB")
geomorphon_profile_dist(ct, margin = 1)
#>            Flat     Slope     Ridge
#> Flat  0.0000000 0.3112781 1.0000000
#> Slope 0.3112781 0.0000000 0.3112781
#> Ridge 1.0000000 0.3112781 0.0000000
```
