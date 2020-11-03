gatherModelMetrics <- function(maxvars = 15) {
  deviances = matrix(ncol = 4, nrow = maxvars)
  logLs = matrix(ncol = 4, nrow = maxvars)
  dfs = matrix(ncol = 4, nrow = maxvars)
  Bics = matrix(ncol = 4, nrow = maxvars)
  colnames(deviances) <- colnames(logLs) <- colnames(dfs) <- colnames(Bics) <- methods <- c('conditionalMutInfo', 'kldivergence', 'kstest', 'meandistance')

  total <- (maxvars-1) * length(methods)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  for (i in 2:maxvars) {
    for (j in 1:length(methods)) {
      setTxtProgressBar(pb, (i-2) * length(methods) + j)
      model <- estimateModel(filtered_data, varRedFcn = get(methods[j]), nvar = i)
      deviances[i, j] <- model$fitted_model$fitinfo$ideviance
      logLs[i, j] <- model$fitted_model$fitinfo$logL
      dfs[i, j] <- model$fitted_model$fitinfo$dimension['df']
      Bics[i, j] <- model$fitted_model$fitinfo$bic
    }
  }



  returnValue(list(
    "deviances" = deviances,
    "logLs" = logLs,
    "dfs" = dfs,
    "Bics" = Bics
  ))
}

makeMetricPlot2 <- function(values, ylabel = 'Value') {
  nn <- ncol(values)
  layout(matrix(c(1, 2), nrow = 1), width = c(4, 1))
  par(mar = c(5, 4, 4, 0)) #No margin on the right side
  matplot(values, type = "l", lwd = 2, ylab = ylabel)
  par(mar = c(5, 0, 4, 2)) #No margin on the left side
  plot(c(0, 1), type = "n", axes = F, xlab = "", ylab = "")
  legend("center", colnames(values), col = seq_len(nn), cex = 0.8, fill = seq_len(nn), lty = 1:nn)
}

makeMetricPlot <- function(values, ylabel = 'Value', xlabel = '# variables') {
  plot(values[,1], type = "l", lwd = 2, ylab = ylabel, xlab = xlabel)
  lines(values[,4], type = "l", lwd = 2, col = 'red', lty = 2)
  lines(values[,2], type = "l", lwd = 2, col = 'green', lty = 3)
  lines(values[,3], type = "l", lwd = 2, col = 'blue', lty = 4)

  legend("topleft", legend=c("Mutual Information", "Means distance", "Kullback-Leibler divergence", "Kolmogorov-Smirnov test"),
       col=c("black", "red", "green", "blue"), lty=1:4, cex=0.8)
}