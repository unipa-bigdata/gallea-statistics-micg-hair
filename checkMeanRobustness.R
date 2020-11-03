library(MASS)
library(ggplot2)

checkMeanRobustness <- function(learn_data, nvar, nsynth, ntrials) {

  robustnesses = matrix(nrow = ntrials, ncol = 4)
  colnames(robustnesses) <- methods <- c('conditionalMutInfo', 'kldivergence', 'kstest', 'meandistance')

  pb <- txtProgressBar(min = 0, max = ntrials * 4, style = 3)
  for (i in 1:length(methods)) {
    model <- estimateModel(learn_data, varRedFcn = get(methods[i]), nvar = nvar)
    for (j in 1:ntrials) {
      setTxtProgressBar(pb, (i - 1) * ntrials + j)
      robustnesses[j, i] <- checkRobustness(model, var_red_fcn = get(methods[i]), n_syn_vars = nsynth)
    }
  }

  returnValue(list(
    'samples' = robustnesses,
    'means' = apply(robustnesses, FUN = mean, MARGIN = 2),
    'vars' = apply(robustnesses, FUN = var, MARGIN = 2)
  ))
}

checkRobustness <- function(model, var_red_fcn, n_syn_vars = 100) {
  learn_data <- model$learn_data
  p <- n_syn_vars
  mu <- runif(p, -10, 10)
  names(mu) <- paste("X", sep = "_", 1:p)
  sigma <- apply(learn_data[, c(2:ncol(learn_data))], 2, var)
  sigma1 <- runif(p, min(sigma), max(sigma))
  Sigma <- diag(p) * sigma1
  dati_gen <- mvrnorm(nrow(learn_data), mu = mu, Sigma = Sigma)

  ## Can we still conduct the same analysis?
  learn_data_new <- data.frame(Code = learn_data[, 1], scale(learn_data[, -1]), dati_gen)
  dim(learn_data_new)

  mis <- var_red_fcn(learn_data_new, 'ET', 'SIC')
  cols <- names(mis)
  returnValue(which(cols == tail(model$variables, 1)))
}

max.na <- function(x) ifelse(!all(is.na(x)), max(x, na.rm = T), NA)

makeRobustnessPlot <- function(data) {
  data <- as.data.frame(data)
  colnames(data) <- c('value')

  p <- ggplot(data, aes(x = value)) +
    geom_histogram()
  print(p)
}
