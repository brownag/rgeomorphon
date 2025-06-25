#' Get a `forms_matrix` for Geomorphon Classification
#'
#' Gets one of the internally defined forms matrices. A form matrix is defined
#' for the classic 10-form output (default) as well as three simplified classes:
#' 4-form, 5-form, and 6-form.
#'
#' @param num_forms Integer. The number of forms to classify, one of `4`, `5`,
#'   `6`, or `10`.
#' @param levels Named integer with values between 0 and 10 corresponding to
#'   form class labels. Default: `get_forms_grass_enum()`
#'
#' @returns An object of class `forms_matrix`
#' @export
forms_matrix_get <- function(num_forms = 10, levels = get_forms_grass_enum()) {
    forms_matrix(get_forms_matrix_cpp(num_forms), levels)
}

#' Apply a `forms_matrix` to Positive and Negative Overlooks
#'
#' @param x SpatRaster containing two layers with names specified in `positive`
#'   and `negative`.
#' @param rcl forms_matrix. Matrix to use for classification of x. Rows are
#'   "negative" and columns are "positive".
#' @param positive Character. Layer name of positive count. Default:
#'   `"positive"`.
#' @param negative Character. Layer name of negative count. Default:
#'   `"negative"`.
#' @param ... Additional arguments passed to [terra::classify()].
#'
#' @returns A SpatRaster containing the classification result.
#' @export
#' @examplesIf requireNamespace("terra")
#' library(terra)
#' library(rgeomorphon)
#'
#' SEARCH = 7       # outer search radius (cells)
#' SKIP = 1         # inner skip radius (cells)
#' DIST = 0         # flatness distance (cells)
#' FLAT = 1         # flat angle threshold
#' MODE = "anglev1" # comparison mode
#'
#' ## classic volcano
#' data("volcano", package = "datasets")
#' dem <- terra::rast(volcano)
#' terra::crs(dem) <- terra::crs("EPSG:2193")
#' terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
#' names(dem) <- "elevation"
#'
#' res <- geomorphons(
#'     dem,
#'     search = SEARCH,
#'     skip = SKIP,
#'     dist = DIST,
#'     flat = FLAT,
#'     comparison_mode = MODE,
#'     forms = TRUE,
#'     ternary = TRUE,
#'     positive = TRUE,
#'     negative = TRUE
#' )
#'
#' res2 <- terra::rast(lapply(c(4, 5, 6), function(n) {
#'   geomorphon_theme(
#'     forms_matrix_apply(
#'         x = res[[c("positive", "negative")]],
#'         rcl = forms_matrix_get(n)
#'     )
#'   )
#' }))
#' names(res2) <- c("forms4", "forms5", "forms6")
#'
#' terra::plot(c(res, res2))
forms_matrix_apply <- function(x,
                               rcl = forms_matrix_get(),
                               positive = "positive",
                               negative = "negative",
                               ...) {
    p <- x[[positive]]
    n <- x[[negative]]
    d <- nrow(rcl)
    r <- n * d + p

    lut <- expand.grid(negative = seq(0, d - 1), positive = seq(0, d - 1))
    lut$form <- unlist(apply(rcl, MARGIN = 2, function(a) a, simplify = FALSE))

    lut <- cbind(ID = lut$negative * d + lut$positive, lut)
    lut <- lut[order(lut$ID), ]
    cls <- terra::subst(r, lut$ID, lut$form)
    names(cls) <- "form"
    cls
}

#' Create a `forms_matrix` object
#'
#' This constructor function wraps a 9x9 integer matrix and associates it
#' with a set of levels, creating a 'forms_matrix' object.
#'
#' @param x Integer. A 9x9  matrix.
#' @param levels Named integer vector. Map of integer values to their string
#'   names. This is typically the output of `get_forms_grass_enum()`.
#'
#' @return An object of class `c("forms_matrix", "matrix", "array")`.
#' @export
forms_matrix <- function(x, levels) {

    if (!is.matrix(x) || !is.numeric(x)) {
        stop("Input 'x' must be a numeric matrix.", call. = FALSE)
    }
    if (!all(dim(x) == c(9, 9))) {
        stop("Input matrix 'x' must be 9x9.", call. = FALSE)
    }
    if (is.null(names(levels)) || !is.numeric(levels)) {
        stop("'levels' must be a named numeric vector.", call. = FALSE)
    }

    if (!all(x %in% levels)) {
        stop("Matrix contains values not found in the provided levels.", call. = FALSE)
    }

    attr(x, "levels") <- levels

    class(x) <- c("forms_matrix", "matrix", "array")

    return(x)
}

#' Print method for a forms_matrix object
#'
#' Controls how the 'forms_matrix' object is displayed in the console.
#'
#' @param x The `forms_matrix` object to print.
#' @param show_values A logical value. If `FALSE` (default), prints enum names.
#'   If `TRUE`, prints the underlying integer values.
#' @param ... Additional arguments passed to print (not used here).
#' @return Invisibly returns the original object `x`.
#' @export
print.forms_matrix <- function(x, show_values = FALSE, ...) {
    levels <- attr(x, "levels")
    display_matrix <- x

    if (!show_values) {
        level_names <- names(levels)
        names(level_names) <- levels
        char_matrix <- matrix(level_names[as.character(x)], nrow = 9, ncol = 9)
        display_matrix <- format(char_matrix, justify = "left")
    } else {
        display_matrix <- format(x, justify = "left")
    }

    rownames(display_matrix) <- paste0("neg=", 0:8)
    colnames(display_matrix) <- paste0("pos=", 0:8)

    cat("<forms_matrix> object\n\n")
    print(display_matrix, quote = FALSE)

    invisible(x)
}
