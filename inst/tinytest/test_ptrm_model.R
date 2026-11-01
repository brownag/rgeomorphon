# Tests for the shared PTRM model helper

library(tinytest)
library(rgeomorphon)

paper <- c(
    flat = 25.12,
    peak = 1.75,
    ridge = 7.59,
    shoulder = 7.31,
    spur = 9.02,
    slope = 28.95,
    hollow = 6.71,
    footslope = 7.63,
    valley = 9.61,
    depression = 3.55
)

# forms10 defaults are the paper weights unchanged
model10 <- .get_ptrm_model("forms10")
expect_equal(model10$coefficients, paper)
expect_identical(model10$classes, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
expect_identical(
    model10$labels,
    c(
        "flat", "peak", "ridge", "shoulder", "spur",
        "slope", "hollow", "footslope", "valley", "depression"
    )
)
expect_equal(model10$intercept, -38.02)
expect_equal(model10$normalizer, 69.96)

# Reduced systems use the reachable GRASS codes in class-code order
model4 <- .get_ptrm_model("forms4")
model5 <- .get_ptrm_model("forms5")
model6 <- .get_ptrm_model("forms6")

expect_identical(model4$classes, c(1, 3, 6, 9))
expect_identical(model4$labels, c("flat", "ridge", "slope", "valley"))
expect_identical(length(model4$coefficients), 4L)
expect_identical(length(model4$labels), 4L)

expect_identical(model5$classes, c(1, 2, 3, 6, 9))
expect_identical(
    model5$labels,
    c("flat", "peak", "ridge", "slope", "valley")
)
expect_identical(length(model5$coefficients), 5L)
expect_identical(length(model5$labels), 5L)

expect_identical(model6$classes, c(1, 3, 4, 6, 8, 9))
expect_identical(
    model6$labels,
    c("flat", "ridge", "shoulder", "slope", "footslope", "valley")
)
expect_identical(length(model6$coefficients), 6L)
expect_identical(length(model6$labels), 6L)

# Each reduced default is the mean of the ten-form paper weights over the
# multiset of reachable (negative, positive) pairs mapping to that class.
# forms4 "slope" is the mean of codes 4,4,4,8,6,6,5,8,8,6,6,6,5,7,6,6,7:
# 278.93 / 17.
expect_equal(model4$coefficients[["slope"]], 16.40764705882353)
expect_equal(model4$coefficients[["ridge"]], 7.13818181818182)
expect_equal(model5$coefficients[["ridge"]], 7.677)
expect_equal(model6$coefficients[["slope"]], 20.11833333333333)
expect_equal(model6$coefficients[["valley"]], 8.61444444444444)

# A full ten-form custom vector is reduced by the same multiset rule
for (system in c("forms4", "forms5", "forms6")) {
    custom <- .get_ptrm_model(system, coefficients = paper)
    derived <- .get_ptrm_model(system)
    expect_equal(custom$coefficients, derived$coefficients)
}

# A reduced-system named vector is accepted case-insensitively and returned in
# class-code order
custom4 <- c(FLAT = 1, ridge = 2, slope = 3, valley = 4)
expect_equal(
    .get_ptrm_model("forms4", coefficients = custom4)$coefficients,
    c(flat = 1, ridge = 2, slope = 3, valley = 4)
)

# Validation errors
expect_error(
    .get_ptrm_model("forms9"),
    "must be one of"
)
expect_error(
    .get_ptrm_model(
        "forms4",
        coefficients = c(flat = 1, ridge = 2, slope = 3, peak = 4)
    ),
    "must be named for"
)
expect_error(
    .get_ptrm_model(
        "forms4",
        coefficients = c(flat = 1, ridge = 2, slope = 3)
    ),
    "must be named for"
)
expect_error(
    .get_ptrm_model(
        "forms4",
        coefficients = c(flat = 1, ridge = 2, slope = NA, valley = 4)
    ),
    "finite"
)
expect_error(
    .get_ptrm_model(
        "forms4",
        coefficients = c(1, 2, 3, 4)
    ),
    "named numeric vector"
)
expect_error(
    .get_ptrm_model("forms10", normalizer = 0),
    "non-zero"
)
expect_error(
    .get_ptrm_model("forms10", intercept = Inf),
    "finite"
)
expect_error(
    .get_ptrm_model(
        "forms10",
        coefficients = paper[names(paper) != "depression"]
    ),
    "must be named for"
)
