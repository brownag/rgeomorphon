# Changelog

## rgeomorphon 0.3.0

CRAN release: 2025-09-16

- Added support for tiled processing for rasters too large to fit in
  memory

- Added vignette: Parallel and Distributed Processing with ‘rgeomorphon’

- Build with `-DARMA_USE_CURRENT` for Armadillo 15+

## rgeomorphon 0.2.0

- Add `positive`, `negative`, and `ternary` arguments and output types
  (for [\#3](https://github.com/brownag/rgeomorphon/issues/3))

- Add support for 4-, 5-, and 6-form classifications via `forms`
  argument

- Add `forms_matrix_*` methods for classifying positive and negative
  output using arbitrary forms matrices (lookup tables)

  - [`forms_matrix_get()`](http://humus.rocks/rgeomorphon/reference/forms_matrix_get.md)
    to obtain forms matrices hard coded in the algorithm

  - [`forms_matrix()`](http://humus.rocks/rgeomorphon/reference/forms_matrix.md)
    to create a new object of class `forms_matrix` (basically a factor
    matrix)

  - [`forms_matrix_apply()`](http://humus.rocks/rgeomorphon/reference/forms_matrix_apply.md)
    to classify positive and negative outputs using a `forms_matrix`

- Add new Salton Sea bathymetry dataset; see
  [`?salton`](http://humus.rocks/rgeomorphon/reference/salton.md)

## rgeomorphon 0.1.0

- Initial GitHub release
