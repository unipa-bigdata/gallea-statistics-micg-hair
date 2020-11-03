library(huge)
library(corrplot)
library(gRim)
library(corrplot)

estimateModel <- function(learn_data, varRedFcn, nvar = 15, type = 'decomposable') {

  distance_metric <- varRedFcn(learn_data, 'ET', 'SIC')

  learn_data <- dropColumns(learn_data, names(distance_metric[(nvar + 1):ncol(learn_data)]))

  glist <- ~.^1
  unconditioned_model <- mmod(glist, data = learn_data)
  unconditioned_model
  plot(unconditioned_model)
  fitted_model <- stepwise(unconditioned_model, details = 0, direction = "forward", type = type)
  fitted_model
  plot(fitted_model)

  ## The canonical parameters are obtained using coef(, type = "ghk")
  ## Get the estimates
  coef_mod <- coef(fitted_model, type = "pms")

  ## p(i)
  pi <- coef_mod$p

  ## m(i)
  mi <- coef_mod$mu
  colnames(mi) <- names(pi)
  rownames(mi) <- fitted_model$coef_mod$cont.names

  i = 2
  plot(1:length(pi), mi[i,], ylab = rownames(mi)[i], axes = FALSE, xlab = "Code")
  Axis(side = 2, labels = TRUE, ylab = rownames(mi)[i])
  axis(1, at = 1:ncol(mi), labels = colnames(mi), las = 1)

  ## K
  K <- round(solve(coef_mod$Sigma), 5)
  colnames(K) <- rownames(K) <- coef_mod$cont.names
  PartialCorrelation <- conc2pcor(K)

  ord <- corrMatOrder(PartialCorrelation, order = "hclust")
  cormat <- PartialCorrelation[ord, ord]
  corrplot(cormat, type = "lower", tl.col = "black")

  returnValue(list(
    "fitted_model" = fitted_model,
    "learn_data" = learn_data,
    "variables" = names(distance_metric[1:nvar]),
    "coef_mod" = coef_mod,
    "K" = K,
    "varRedFcn" = varRedFcn
  ))
}