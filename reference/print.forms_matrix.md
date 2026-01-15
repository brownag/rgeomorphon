# Print method for a forms_matrix object

Controls how the 'forms_matrix' object is displayed in the console.

## Usage

``` r
# S3 method for class 'forms_matrix'
print(x, show_values = FALSE, ...)
```

## Arguments

- x:

  The `forms_matrix` object to print.

- show_values:

  A logical value. If `FALSE` (default), prints enum names. If `TRUE`,
  prints the underlying integer values.

- ...:

  Additional arguments passed to print (not used here).

## Value

Invisibly returns the original object `x`.

## Examples

``` r
print(forms_matrix_get(num_forms = 4))
#> <forms_matrix> object
#> 
#>       pos=0  pos=1  pos=2  pos=3  pos=4  pos=5  pos=6  pos=7  pos=8 
#> neg=0 G_FL   G_FL   G_FL   G_SL   G_VL   G_VL   G_VL   G_VL   G_VL  
#> neg=1 G_FL   G_FL   G_SL   G_SL   G_VL   G_VL   G_VL   G_VL   G_NONE
#> neg=2 G_FL   G_SL   G_SL   G_SL   G_SL   G_VL   G_VL   G_NONE G_NONE
#> neg=3 G_SL   G_SL   G_SL   G_SL   G_SL   G_SL   G_NONE G_NONE G_NONE
#> neg=4 G_RI   G_RI   G_SL   G_SL   G_SL   G_NONE G_NONE G_NONE G_NONE
#> neg=5 G_RI   G_RI   G_RI   G_SL   G_NONE G_NONE G_NONE G_NONE G_NONE
#> neg=6 G_RI   G_RI   G_RI   G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE
#> neg=7 G_RI   G_RI   G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE
#> neg=8 G_RI   G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE G_NONE
```
