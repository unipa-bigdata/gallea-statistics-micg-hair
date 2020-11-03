library(infotheo)
library(ggplot2)
library(plyr)
library(gridExtra)
library(gRim)
library(rlang)
library(msos)
library(e1071)
library(LaplacesDemon)
library(dgof)

dropColumns <- function(data, cols) data[, !(names(data) %in% cols)]

showVarPlot <- function(data, y, x_var, type) {
  a <- ggplot(data, aes_string(x = x_var))
  #mu <- ddply(data, "Code", summarise, grp.mean = mean(Co))
  return(a +
           (if (type == 'density')
             geom_density(aes(color = Code))
           else
             geom_boxplot(aes(color = Code))) +
           #geom_vline(data = mu, aes(xintercept = grp.mean, color = Code), linetype = "dashed") +
           theme(legend.position = "none"))
}

whichpart <- function(x, n = 30) {
  nx <- length(x)
  p <- nx - n
  xp <- sort(x, partial = p)[p]
  which(x > xp)
}

crossMutInfo <- function(data) {

  nvar <- ncol(data) - 1

  mi = matrix(nrow = nvar, ncol = nvar)

  for (i in 1:nvar) {
    for (j in 1:i) {
      if (j == i) mi[i, j] <- 0
      else {
        var1 <- discretize(data[, i + 1])
        var2 <- discretize(data[, j + 1])
        mi[i, j] <- mi[j, i] <-
          mutinformation(var1, var2) /
            (sqrt(entropy(var1)) * sqrt(entropy(var2)))
      }
    }
  }
  colnames(mi) = colnames(data[, -1])
  rownames(mi) = colnames(data[, -1])

  heatmap(exp(mi), Colv = NA, Rowv = NA, scale = "column", xlab = "variable 1", ylab = "variable 2", main = "heatmap")

  returnValue(mi)
}

whichpart <- function(x, n = 20) {
  nx <- length(x)
  p <- nx - n
  xp <- sort(x, partial = p)[p]
  returnValue(which(x > xp))
}

selectVariables <- function(mi_matrix, n = 20) {
  max_vector <- apply(mi_matrix, 1, sum)
  print(max_vector)
  #max_idxs <- whichpart(max_vector, n)
  max_vector = sort(max_vector, decreasing = TRUE)

  returnValue(list(
    'max_mi' = max_vector,
    'less_informative_vars' = names(max_vector[1:n])
  ))
}

conditionalPlots <- function(data, type) {

  plots = list()
  var_names <- colnames(data[, -1])
  i <- 1
  for (var in var_names) {
    plot <- showVarPlot(data, y = "Code", x_var = var, type)
    plots[[i]] <- plot

    if ((i %% 16) == 0) {
      grid.arrange(grobs = duplicate(plots), ncol = 4, nrow = 4)
      plot.new()
      i <- 0
      plots <- list()
    }

    i <- i + 1
  }
}

conditionalMutInfo <- function(data, var1, var2) returnValue(sort(distanceFunction(data, var1, var2, nmi), decreasing = FALSE))

kldivergence <- function(data, var1, var2) returnValue(distanceFunction(data, var1, var2, kld))

kstest <- function(data, var1, var2) returnValue(distanceFunction(data, var1, var2, ksd))

meandistance <- function(data, var1, var2) returnValue(distanceFunction(data, var1, var2, means))

means <- function(x1, x2) {
  returnValue(abs(mean(x1) - mean(x2)))
}

ksd <- function(x1, x2) {
  options(warn=-1)
  res <- ks.test(x1, x2)$statistic
  options(warn=0)
  returnValue(res)
}

nmi <- function(x1, x2) {
  d1 <- discretize(density(x1)$y)
  d2 <- discretize(density(x2)$y)

  returnValue(mutinformation(d1, d2) /
                (sqrt(entropy(d1)) * sqrt(entropy(d2))))
}

kld <- function(x1, x2) {
  returnValue(KLD(
    density(x1)$y,
    density(x2)$y
  )$mean.sum.KLD)
}

distanceFunction <- function(data, var1, var2, fun) {
  mis <- c()
  for (i in colnames(data[, -1])) {
    d1 <- discretize(density(data[, i][data$Code == var1])$y)
    d2 <- discretize(density(data[, i][data$Code == var2])$y)
    #mis[i] <- mutinformation(d1, d2) /
    #  (sqrt(entropy(d1)) * sqrt(entropy(d2)))
    mis[i] <- fun(data[, i][data$Code == var1], data[, i][data$Code == var2])
  }
  names(mis) <- colnames(data[, -1])
  mis <- sort(mis, decreasing = TRUE)
  plot(mis)
  returnValue(mis)
}

negentropyTest <- function(data, classes, thresh = 0.1) returnValue(metricTest(negent, data, classes, thresh, 'Negentropy'))

kurtosisTest <- function(data, classes, thresh = 0.1) returnValue(metricTest(kurtosis, data, classes, thresh, 'Kurtosis'))

metricTest <- function(fun, data, classes, thresh = 0.1, label) {
  res <- matrix(nrow = length(classes), ncol = ncol(data) - 1)
  rownames(res) = classes
  colnames(res) = colnames(data[, -1])
  for (i in names(data[, -1])) {
    for (j in classes) {
      res[j, i] <- fun(data[data$Code == j, i])
    }
  }

  par(las = 2)
  if (thresh == 'auto') {
    thresh = mean(res) + 3 * sd(res)
  }

  colors = c("black", "gray", "lightgray")
  pchs = c(15, 18, 17)

  matplot(t(res), xlab = 'Elements', ylab = label, pch = pchs, axes = FALSE, col = colors)
  #show <- which(res > thresh)
  #print(show)
  show = c()
  for (i in classes) {
    show <- c(show, which(res[i,] > thresh))
    #plot(res[i,], xlab = 'Elements', ylab = 'Negentropy', pch = 20, axes = FALSE, col = colors[i])
  }
  axis(1, at = show, labels = colnames(data[i, show + 1]))
  axis(2)

  returnValue(list(
    "excluded_vars" = unique(colnames(data[i, show + 1])),
    "included_vars" = colnames(data[i, -c(1, show + 1)]),
    "values" = res
  ))
}

communalityTest <- function(data, threshold = 0.1) {
  data = data[, -1]

  data.fa <- factanal(data[, -1], factors = 1)
  sq_loadings = sort(apply(data.fa$loadings, 1, function(x) sum(x^2)), decreasing = TRUE)

  returnValue(list(
    'correlated_variables' = names(sq_loadings[sq_loadings >= threshold]),
    'uncorrelated_variables' = names(sq_loadings[sq_loadings < threshold]),
    'loadings' = sq_loadings
  ))


}