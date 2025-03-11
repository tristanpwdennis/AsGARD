library(tidyverse)
library(OptM)

install.packages('OptM')
folder <- file.path(path="~/Projects/AsGARD/data/TreeMix_20240926/infer_migration/")                     #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points'

test.optM = optM(folder, tsv ="Evanno.variance.txt")              #another option is the Evanno method - see optM package description for detailed information on output
#if data is robust and all runs have the same likelihoods, SD will be 0 and this function will give an error as it can't produce the ad hoc statistic. 
#in this case you might want to increase variance by varying -k (SNP block size), change permutation methods etc.
plot_optM(test.optM, method = "Evanno")                           #plot the proportion of variation explained by each migration event. Calculates deltaM, which is a second-order rate of change in likelihood weighted by the standard deviation
