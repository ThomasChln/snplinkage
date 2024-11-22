#' @import gdsfmt
#' @import ggplot2
#' @import GWASTools
#' @importFrom biomaRt getBM useMart
#' @importFrom cowplot theme_cowplot
#' @importFrom data.table fread
#' @importFrom ggrepel geom_label_repel
#' @importFrom grDevices grey rgb
#' @importFrom grid grid.draw unit
#' @importFrom gtable gtable_add_grob gtable_add_rows gtable_filter
#' @importFrom knitr kable
#' @importFrom methods as callNextMethod is new
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
#' @importFrom SNPRelate snpgdsIBS snpgdsLDMat snpgdsLDpruning snpgdsOpen
#'                       snpgdsSampMissRate snpgdsSelectSNP snpgdsSNPRateFreq
#' @importFrom stats chisq.test na.omit p.adjust setNames
#' @importFrom utils globalVariables head tail unzip
NULL

#' Pipe
#'
#' Pipe an object forward into a function or call expression.
#' Magrittr imported function, see details and examples in the magrittr package.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return Result of rhs applied to lhs, see details in magrittr package.
#' @export
NULL

#' Assignment pipe
#'
#' Pipe an object forward into a function or call expression and update the
#' `lhs` object with the resulting value.
#' Magrittr imported function, see details and examples in the magrittr package.
#'
#' @importFrom magrittr %<>%
#' @name %<>%
#' @rdname compound
#' @param lhs An object which serves both as the initial value and as target.
#' @param rhs a function call using the magrittr semantics.
#' @return None, used to update the value of lhs.
#' @export
NULL

#' Exposition pipe
#'
#' Expose the names in `lhs` to the `rhs` expression.
#' Magrittr imported function, see details and examples in the magrittr package.
#'
#' @importFrom magrittr %$%
#' @name %$%
#' @rdname exposition
#' @param lhs A list, environment, or a data.frame.
#' @param rhs An expression where the names in lhs is available.
#' @return Result of rhs applied to one or several names of lhs.
#' @export
NULL

