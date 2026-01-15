# Get a `forms_matrix` for Geomorphon Classification

Gets one of the internally defined forms matrices. A form matrix is
defined for the classic 10-form output (default; Jasiewicz & Stepinski,
2013) as well as three simplified classes: 4-form, 5-form, and 6-form
(Masetti et al., 2018)

## Usage

``` r
forms_matrix_get(num_forms = 10, levels = get_forms_grass_enum())
```

## Arguments

- num_forms:

  Integer. The number of forms to classify, one of `4`, `5`, `6`, or
  `10` (default).

- levels:

  Named integer with values between 0 and 10 corresponding to form class
  labels. Default:
  [`get_forms_grass_enum()`](http://humus.rocks/rgeomorphon/reference/forms_matrix.md)

## Value

An object of class `forms_matrix`

## Details

For creating custom classification systems see the
[`forms_matrix()`](http://humus.rocks/rgeomorphon/reference/forms_matrix.md)
constructor.

## References

Stepinski, T., Jasiewicz, J., 2011, Geomorphons - a new approach to
classification of landform, in : Eds: Hengl, T., Evans, I.S., Wilson,
J.P., and Gould, M., Proceedings of Geomorphometry 2011, Redlands,
109-112. Available online:
<https://www.geomorphometry.org/uploads/pdf/pdf2011/StepinskiJasiewicz2011geomorphometry.pdf>

Jasiewicz, J., Stepinski, T., 2013, Geomorphons - a pattern recognition
approach to classification and mapping of landforms, Geomorphology, vol.
182, 147-156.
([doi:10.1016/j.geomorph.2012.11.005](https://doi.org/10.1016/j.geomorph.2012.11.005)
)

Masetti, G., Mayer, L. A., & Ward, L. G. 2018, A Bathymetry- and
Reflectivity-Based Approach for Seafloor Segmentation. Geosciences,
8(1), 14.
([doi:10.3390/geosciences8010014](https://doi.org/10.3390/geosciences8010014)
)

## See also

[`forms_matrix()`](http://humus.rocks/rgeomorphon/reference/forms_matrix.md)

## Examples

``` r
forms_matrix_get()
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
```
