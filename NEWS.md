# rgeomorphon 0.4.0

 - Added `geomorphon_ptrm()` and `geomorphon_ptrm_focal()` to measure perceived terrain realism from a DEM or geomorphon map

   - Supports global and circular focal-window calculations, `forms4`/`forms5`/`forms6`/`forms10` classifications, and tunable model parameters (`coefficients`, `intercept`, `normalizer`).

 - Added `geomorphon_concordance()` to compare map structure between geomorphon and categorical rasters

   - Includes Normalized Mutual Information (`nmi`) with directional `homogeneity` and `completeness`, directional boundary-alignment summaries (`edge_offset`, `edge_concordance`), symmetric shared-domain boundary metrics (`edge_f1`, `edge_mhd`, `edge_chamfer`), and an optional deterministic 128-shift toroidal `spatial_null` ensemble.

 - Added the `scorpion` sample dataset (10 m USGS 3DEP DEM and co-registered SSURGO `mukey` raster and map unit attributes) and the vignette *Terrain Realism (PTRM) and Structural Concordance with 'rgeomorphon'*.

 - Fixed simplified `geomorphons(..., simplify = TRUE)` output when given a `SpatRaster` input and added internal information-theory helpers and tests; no new package dependencies.

# rgeomorphon 0.3.0

 - Added support for tiled processing for rasters too large to fit in memory
 
 - Added vignette: Parallel and Distributed Processing with 'rgeomorphon'
 
 - Build with `-DARMA_USE_CURRENT` for Armadillo 15+
 
# rgeomorphon 0.2.0
 
 - Add `positive`, `negative`, and `ternary` arguments and output types (for #3)
 
 - Add support for 4-, 5-, and 6-form classifications via `forms` argument
 
 - Add `forms_matrix_*` methods for classifying positive and negative output using arbitrary forms matrices (lookup tables)
 
   - `forms_matrix_get()` to obtain forms matrices hard coded in the algorithm
 
   - `forms_matrix()` to create a new object of class `forms_matrix` (basically a factor matrix)
 
   - `forms_matrix_apply()` to classify positive and negative outputs using a `forms_matrix`
 
 - Add new Salton Sea bathymetry dataset; see `?salton`

# rgeomorphon 0.1.0

* Initial GitHub release
