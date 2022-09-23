library(testthat)

mfdList <- listAvailMfd()
mfdFinite <- mfdList[vapply(mfdList, is.finiteDim, TRUE)]
