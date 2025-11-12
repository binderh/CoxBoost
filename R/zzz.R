#' @useDynLib CoxBoost, .registration = TRUE
NULL

#' @import Matrix
#' @importFrom survival Surv survfit
#' @importFrom graphics abline axis legend lines points text
#' @importFrom stats as.dist cor hclust heatmap model.frame model.matrix model.response pchisq predict sd terms
#' @importFrom grDevices gray
"_PACKAGE"
