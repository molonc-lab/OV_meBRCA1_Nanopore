#illumina manual processing EPIC v2 method from minfi
.highestType <- function(...) {
  dots <- list(...)
  types <- vapply(dots, DelayedArray::type, character(1L))
  vector <- do.call(c, lapply(types, vector))  # guaranteed to be atomic
  typeof(vector)
}

.NA_type <- function(type) {
  c(vector(type), NA)
}
setGeneric(
  ".preprocessRawMeth",
  function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...)
    standardGeneric(".preprocessRawMeth"),
  signature = c("Red", "Green"))

# `...` are additional arguments passed to methods.
setGeneric(
  ".preprocessRawUnmeth",
  function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...)
    standardGeneric(".preprocessRawUnmeth"),
  signature = c("Red", "Green"))

# `...` are additional arguments passed to methods.
setGeneric(
  ".preprocessRaw",
  function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...)
    standardGeneric(".preprocessRaw"),
  signature = c("Red", "Green"))

setMethod(
  ".preprocessRawMeth",
  c("matrix", "matrix"),
  function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII) {
    # Set up output matrix with appropriate dimension and type
    type <- .highestType(Red, Green)
    M <- matrix(.NA_type(type),
                nrow = length(locusNames),
                ncol = ncol(Red),
                dimnames = list(locusNames, colnames(Red)))
    
    # Fill output matrix
    M[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressB, ]
    M[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressB, ]
    M[TypeII$Name, ] <- Green[TypeII$AddressA, ]
    
    # Return output matrix
    M
  }
)

setMethod(
  ".preprocessRawUnmeth",
  c("matrix", "matrix"),
  function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...) {
    # Set up output matrix with appropriate dimension and type
    type <- .highestType(Red, Green)
    U <- matrix(.NA_type(type),
                nrow = length(locusNames),
                ncol = ncol(Red),
                dimnames = list(locusNames, colnames(Red)))
    
    # Fill output matrix
    U[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressA, ]
    U[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressA, ]
    U[TypeII$Name, ] <- Red[TypeII$AddressA, ]
    
    # Return output matrix
    U
  }
)

setMethod(
  ".preprocessRaw",
  c("matrix", "matrix"),
  function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...) {
    list(
      M = .preprocessRawMeth(
        Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII),
      U = .preprocessRawUnmeth(
        Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII))
  }
)

blockMapplyWithRealization <- function(FUN, ..., MoreArgs = NULL, sinks = NULL,
                                       dots_grids = NULL, sinks_grids = NULL,
                                       BPREDO = list(), BPPARAM = bpparam()) {
  FUN <- match.fun(FUN)
  dots <- unname(list(...))
  # Check valid `sinks`
  stopifnot(is.null(sinks) || is.list(sinks))
  # Check conformable dots_grids and sinks_grids
  if (is.null(dots_grids)) {
    dots_grids <- replicate(length(dots), NULL)
  } else {
    stopifnot(is.list(dots_grids))
  }
  dots_grids <- mapply(
    FUN = DelayedArray:::normarg_grid,
    dots_grids,
    dots,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE)
  if (is.null(sinks_grids)) {
    sinks_grids <- replicate(length(sinks), NULL)
  } else {
    stopifnot(is.list(sinks_grids))
  }
  sinks_grids <- mapply(
    FUN = DelayedArray:::normarg_grid,
    sinks_grids,
    sinks,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE)
  grids_dims <- lapply(c(dots_grids, sinks_grids), dim)
  all_same_grids_dims <- all(
    vapply(X = grids_dims,
           FUN = function(dim) all(dim == grids_dims[[1L]]),
           FUN.VALUE = logical(1L)))
  if (!all_same_grids_dims) {
    stop("non-conformable 'dots_grids' and 'sinks_grids'")
  }
  stopifnot(length(dots) == length(dots_grids),
            length(sinks) == length(sinks_grids))
  
  # Loop over blocks of `dots` and write to `sinks`
  nblock <- length(dots_grids[[1]])
  bplapply(seq_len(nblock), function(b) {
    if (DelayedArray:::get_verbose_block_processing()) {
      message("Processing block ", b, "/", nblock, " ... ",
              appendLF = FALSE)
    }
    input_viewports <- lapply(dots_grids, function(grid) grid[[b]])
    output_viewports <- lapply(sinks_grids, function(grid) grid[[b]])
    blocks <- mapply(
      FUN = function(x, grid, viewport) {
        block <- read_block(x, viewport)
        attr(block, "from_grid") <- grid
        attr(block, "block_id") <- b
        block
      },
      x = dots,
      grid = dots_grids,
      viewport = input_viewports,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE)
    block_ans <- do.call(FUN, c(blocks, MoreArgs))
    if (!is.list(block_ans)) {
      block_ans <- list(block_ans)
    }
    # NOTE: This is the only part different from blockMapply()
    if (!is.null(sinks)) {
      mapply(function(ans, sink, viewport) {
        write_block(sink, viewport, ans)
      }, ans = block_ans, sink = sinks, viewport = output_viewports,
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
      block_ans <- NULL
    }
    if (DelayedArray:::get_verbose_block_processing()) {
      message("OK")
    }
    block_ans
  },
  BPREDO = ,
  BPPARAM = bpparam())
}

