#' Create RFPCA input from a data frame
#'
#' @param data A data frame
#' @param IDvar The column name that stands for the ID
#' @param tvar The column name that stands for the time
#' @param yvar The column names that stands for the values for Lt. If not specified, all columns except for IDvar and tvar are used.
#' @export

MakeInput <- function(data, IDvar, tvar, yvar) {
  ID <- data[[IDvar]]
  tt <- data[[tvar]]
  if (missing(yvar)) {
    yvar <- setdiff(names(data), c(IDvar, tvar))
  }
  yy <- data[, yvar]
  data <- data[order(data[, IDvar, drop=TRUE],
                     data[, tvar, drop=TRUE] ), ]

  Ly <- plyr::dlply(data, IDvar, function(d1) {
    dd <- t(as.matrix(d1[, yvar, drop=FALSE]))
    dd
  })
  Lt <- plyr::dlply(data, IDvar, function(d1) {
    d1[, tvar, drop=TRUE]
  })
  attr(Ly, 'split_type') <- NULL
  attr(Ly, 'split_labels') <- NULL
  attr(Lt, 'split_type') <- NULL
  attr(Lt, 'split_labels') <- NULL
  list(Lid=as.list(unique(data[, IDvar, drop=TRUE])), 
       Ly=Ly, Lt=Lt)
}
