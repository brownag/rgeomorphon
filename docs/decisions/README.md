# Architecture Decision Records

This directory is the living catalog of enduring architectural decisions in
`rgeomorphon`. Records use a MADR-compatible format with YAML front matter.
Statuses are `Proposed`, `Accepted`, `Deprecated`, or `Superseded`.

| ID | Title | Status | Date | Supersedes | Superseded By | Context |
|---|---|---|---|---|---|---|
| [0001](0001-cpp-rcppparallel-engine.md) | Use a C++ RcppParallel worker for geomorphon classification | Accepted | 2025-06-03 | - | - | Initial public commit; C++ engine |
| [0002](0002-terra-optional-matrix-core.md) | Keep runtime dependencies minimal and make terra optional | Accepted | 2025-06-03 | - | - | Initial public commit; dependency policy |
| [0003](0003-forms-matrix-classification.md) | Represent reduced form systems as composable classification matrices | Accepted | 2025-06-25 | - | - | Reduced form-system work |
| [0004](0004-tiled-raster-processing.md) | Process large rasters as buffered tiles with two-level parallelism | Accepted | 2025-06-30 | - | - | Tiled processing milestone |
| [0005](0005-grass-class-codes-theme.md) | Use GRASS-compatible integer class codes with a shared display theme | Accepted | 2025-06-03 | - | - | Initial public commit; class taxonomy |
