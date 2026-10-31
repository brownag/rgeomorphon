# rgeomorphon 0.4.0

 - Added `geomorphon_ptrm()` for calculating the Perceived Terrain Realism Metric (PTRM) from a DEM or an existing geomorphon classification

 - Added `geomorphon_ptrm_focal()` for spatially varying PTRM calculated over a circular moving window

 - Added `geomorphon_concordance()` for measuring structural agreement between geomorphon and categorical maps

 - Support for `forms4`, `forms5`, and `forms6` form systems in the PTRM functions

 - Tunable `intercept` and `normalizer` PTRM model parameters

 - Corrected focal PTRM window geometry, border handling, and radius interpretation

 - Corrected concordance distance units and boundary metrics; removed the stale `...` argument

# rgeomorphon 0.3.2

 - Added `geomorphon_ptrm()` function to calculate Perceived Terrain Realism Metric (Rajasekaran et al. 2022)
   - Quantifies human-perceived terrain realism on a 0-1 scale
   - Supports calculation from DEMs or pre-computed geomorphon classifications
   - Based on proportional distribution of 10 standard geomorphon types

 - Added `geomorphon_concordance()` function to measure structural agreement between geomorphon and categorical maps
   - Normalized Mutual Information (NMI) for statistical association
   - Boundary alignment metrics (edge offset and edge concordance)
   - No assumptions about 1-to-1 class equivalency between maps

 - Internal helper functions for entropy and mutual information calculations (no new package dependencies)

 - Comprehensive unit tests for new analysis functions

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
