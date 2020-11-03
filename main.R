rm(list = ls())


source('functions.R')
source('estimateModel.R')
source('gatherModelMetrics.R')
source('checkMeanRobustness.R')

load('data/dataset.Rdata')

# sort variables by elements name
data <- cbind(dati[, 1:2], dati[, sort(colnames(dati[, -c(1, 2)]))])

# drop Gender variable
filtered_data <- dropColumns(data, 'Sesso')

# scale data (not required since values are already scaled)
filtered_data <- data.frame(Code = filtered_data[, 1], scale(filtered_data[, -1]))

summary(filtered_data)

## conditional plot, for assessing quasi-normality
conditionalPlots(data = filtered_data, type = 'density')    # conditional pdf plots
#conditionalPlots(data = filtered_data, type = 'boxplot')   # boxplots, not used in report

# Gaussianity tests
negentropy <- negentropyTest(filtered_data, c('ET', 'SIC'), thresh = 'auto')
kurtosity <- kurtosisTest(filtered_data, c('ET', 'SIC'), thresh = 'auto')

# Communality test, to determine which variables are not correlated, it relies on 1-factor factor analysis, and sorts
# variables by loading value, which is a measure or communality, i.e. the amount of explained variance in common with
# other variables
factor_communalities <- communalityTest(filtered_data, threshold = 0.04)
factor_communalities$uncorrelated_variables
factor_communalities$correlated_variables

# Drop elements which are known to be synthetic, thus cannot be related to volcanic activity
synthetic_elements <- c('Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'Ls', 'No', 'Rf', 'Db', 'Sg', 'Hs', 'Mt', 'Ds', 'Rg',
                        'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'Bh')
filtered_data <- dropColumns(filtered_data, synthetic_elements)

# Some other elements are very rare in nature but they are not filtered out. However, we expect the estimated model not
# to include them
#very_uncommon_elements <- c('Tc', 'Ra', 'Sm', 'Md')
#filtered_data <- dropColumns(filtered_data, very_uncommon_elements)

# Model estimation, using different variable reduction methods

# Conditional Mutual Information
m1 <- estimateModel(filtered_data, conditionalMutInfo)
summary(m1$fitted_model)
plot(m1$fitted_model)

# Means distance
m2 <- estimateModel(filtered_data, meandistance)
summary(m2$fitted_model)
plot(m2$fitted_model)

# Kullback-Leibler Divergence
m3 <- estimateModel(filtered_data, kldivergence)
summary(m3$fitted_model)
plot(m3$fitted_model)

# Kolmogorov-Smirnov Test
m4 <- estimateModel(filtered_data, kstest)
summary(m4$fitted_model)
plot(m4$fitted_model)

####################################################################################################################
#
# From this section computation are very time consuming, uncomment these, and come back after a couple of hours... :D
# There are progress bars, but they don't progress linearly, due to increasing model complexity
#
# As a (suggested) alternative, use the provided pre-computed workspace image
#
####################################################################################################################

load('ws.Rdata')

#########################
# Performance measurement
#########################
#
# Generate various performance plots for the different variable reduction methods
# Metrics shown are Deviance, LogL, DF size and BIC

#metrics <- gatherModelMetrics(maxvars = 20)
makeMetricPlot(metrics$deviances, ylabel = 'Deviance')
makeMetricPlot(metrics$logLs, ylabel = 'LogL')
makeMetricPlot(metrics$dfs, ylabel = 'DF size')
makeMetricPlot(metrics$Bics, ylabel = 'BIC')


#########################
# Robustness measurement
#########################
#
# After injecting nsynth synthetic variables X_i, robustness is measured as the amount of variables required before
# including the whole real variables subset. For example, if nvar = 3 and and nsynth = 2, variables are 'Bo', 'Sr' and
# 'V' and the procedure injects `X_1', 'X_2', if variable reduction sorts variables as `Bo', 'X_1', 'V', 'Sr', 'X_2',
# robustness measure is 4, because the initial subset is complete after taking the first 4 variables. Hence, lower
# values correspond to better robustness.
#
# The range of robustness value is [nvar, nvar+nsynth]
#
# A simulation is ran 1000 times using nsytnh = 100 synthetic variables and models with nvar=15. The resulting
# robustness is binned in an histogram.

#robustnesses <- checkMeanRobustness(learn_data = filtered_data, nvar = 15, nsynth = 100, ntrials = 1000)
makeRobustnessPlot(robustnesses$samples[,1]) # Mutual Information
makeRobustnessPlot(robustnesses$samples[,2]) # Kullback-Leibler Divergence
makeRobustnessPlot(robustnesses$samples[,3]) # Kolmogorov-Smirnov Test
makeRobustnessPlot(robustnesses$samples[,4]) # Means distance

#save.image('ws.Rdata')

