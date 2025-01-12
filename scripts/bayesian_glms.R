# Load libs
library(ggplot2)
library(rstanarm)
library(loo)
library(projpred)
library(posterior)
library(bayesplot)
library(corrplot)
#Running Bayesian GLMs, performing model fits, and selecting important varaibles, using RSTANARM
# This script mostly follows the vignette of: https://mc-stan.org/projpred/articles/projpred.html#modtypes


# Load fst, the response variable

# Sample metadata
df_samples <- read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_qcpass.20240914.csv')
points <- df_samples %>%
  dplyr::filter(country %in% c('Djibouti', 'Sudan', 'Ethiopia')) %>%
  dplyr::select(longitude, latitude, location, country) %>%
  unique()

# Load data
dd_analysis <- read.csv('~/Projects/AsGARD/data/environmental_modelling/environmental_data.csv')#Drop columns we don't want (cattle has too much missing data)

#Corr mat
tiff("~/Projects/AsGARD/figures/corrplot.tiff", units="in", width=5, height=5, res=300)
corrplot(cor(dd_analysis))
dev.off()

#Specify full model
# Set this manually if desired:
ncores <- parallel::detectCores(logical = FALSE)

# Now try using horseshow peior, and different model formulations
n=dim(dd_analysis)[1]
p=dim(dd_analysis)[2]
p0 <- 2 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
hs_prior <- hs(df=1, global_df=1, global_scale=tau0)
t_prior <- student_t(df = 7, location = 0, scale = 2.5)

post1 <- stan_glm(formula=fst_lin ~ pdist +mean.wind.traveltime+mdr.circuit + evi.circuit+pop.circuit+ele.circuit+hurs.circuit+friction.circuit, 
                  data = dd_analysis,
                  family = gaussian(link = "identity"),
                  prior = hs_prior, prior_intercept = t_prior,
                  seed = 121347, adapt_delta = 0.999,
                  chains = 4, iter = 4000,
                  QR = TRUE, refresh = 0)
plot(post1)
(loo1 <- loo(post1))

post2<-stan_glm(fst_lin ~ pdist+pop.circuit*pdist+mdr.circuit*pdist+slope.circuit*pdist+mean.wind.traveltime*pdist+ evi.circuit*pdist+hurs.circuit*pdist+ele.circuit*pdist,
    data = dd_analysis,
    family = gaussian(link = "identity"),
    prior = hs_prior, prior_intercept = t_prior,
    seed = 121347, adapt_delta = 0.999,
    chains = 4, iter = 4000,
    QR = TRUE, refresh = 0)

plot(post2)
(loo2 <- loo(post2))

loo_compare(loo1, loo2) #controlling for distance does not improve model fit

#using the prior improves loo?

# VAREIABLE SELECTION (takes ages)
varsel2 <- cv_varsel(post1, method='forward', cv_method='loo', nloo = n, nclusters_pred=100, ndraws_pred=2000, nterms_max=15, parallel=FALSE) #run variable selection (this takes quite a long time)
(nsel<-suggest_size(varsel2))
plot(varsel2, stats = "mlpd", deltas = TRUE)+theme_classic()

n_sel <- suggest_size(varsel2, stat = "mlpd") #get suggested no vars

rk <- ranking(varsel2) #rank variable importance and plot
( pr_rk <- cv_proportions(rk) )  
plot(pr_rk)
plot(cv_proportions(rk, cumulate = TRUE))
(predictors_final <- head(rk[["fulldata"]], n_sel+1) ) #get final predictors (eg 3 plus md because we want ti)
prj <- project(
  post2,
  predictor_terms = predictors_final,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = FALSE
  ###
)

# Examine the down-projected model
prj_mat <- as.matrix(prj)
prj_drws <- as_draws_matrix(prj_mat)
prj_smmry <- summarize_draws(
  prj_drws,
  "median", "mad", function(x) quantile(x, probs = c(0.025, 0.975))
)
# Coerce to a `data.frame` because pkgdown versions > 1.6.1 don't print the
# tibble correctly:
prj_smmry <- as.data.frame(prj_smmry)
print(prj_smmry, digits = 1)

# Using the 'bayesplot' package, plot predicted values of fst vs actual
prj_predict <- proj_predict(prj)
ppc_dens_overlay(y = dd_analysis$fst_lin, yrep = prj_predict)
bayesplot_theme_set(ggplot2::theme_bw())


#so a final sensible model looks like it may be
postfinal <- stan_glm(formula=fst_lin ~  mean.wind.traveltime+mdr.circuit + evi.circuit+friction.circuit,
                  data = dd_analysis,
                  family = gaussian(link = "identity"),
                  prior = hs_prior, prior_intercept = t_prior,
                  seed = 121347, adapt_delta = 0.999,
                  chains = 4, iter = 4000,
                  QR = TRUE, refresh = 0)
posterior_plot <- mcmc_areas(as.matrix(postfinal),
                             pars = names(fixef(postfinal)),
                             prob = 0.95) + 
ggtitle("Posterior Distributions with 95% Credible Intervals")
posterior_plot


##### plotting model predictions

#individual submodels + predplots

m.wind <- stan_glm(formula=fst_lin ~  mean.wind.traveltime,
                      data = dd_analysis,
                      family = gaussian(link = "identity"),
                      prior = hs_prior, prior_intercept = t_prior,
                      seed = 121347, adapt_delta = 0.999,
                      chains = 4, iter = 4000)
draws <- as.data.frame(as.matrix(m.wind))
colnames(draws)[1:2] <- c("a", "b")
ggplot(dd_analysis, aes(x = mean.wind.traveltime, y = fst_lin)) +
  geom_point(size = 1) +
  geom_abline(data = draws, aes(intercept = a, slope = b),
              color = "skyblue", size = 0.2, alpha = 0.25) +
  geom_abline(intercept = coef(m.wind)[1], slope = coef(m.wind)[2],
              color = "skyblue4", size = 1)+
  theme_classic()



m.mdr <- stan_glm(formula=fst_lin ~  mdr.circuit,
                   data = dd_analysis,
                   family = gaussian(link = "identity"),
                   prior = hs_prior, prior_intercept = t_prior,
                   seed = 121347, adapt_delta = 0.999,
                   chains = 4, iter = 4000)
draws <- as.data.frame(as.matrix(m.mdr))
colnames(draws)[1:2] <- c("a", "b")
ggplot(dd_analysis, aes(x = mdr.circuit, y = fst_lin)) +
  geom_point(size = 1) +
  geom_abline(data = draws, aes(intercept = a, slope = b),
              color = "skyblue", size = 0.2, alpha = 0.25) +
  geom_abline(intercept = coef(m.mdr)[1], slope = coef(m.mdr)[2],
              color = "skyblue4", size = 1)+
  theme_classic()

m.evi<- stan_glm(formula=fst_lin ~  evi.circuit,
                  data = dd_analysis,
                  family = gaussian(link = "identity"),
                  prior = hs_prior, prior_intercept = t_prior,
                  seed = 121347, adapt_delta = 0.999,
                  chains = 4, iter = 4000)
draws <- as.data.frame(as.matrix(m.evi))
colnames(draws)[1:2] <- c("a", "b")
ggplot(dd_analysis, aes(x = evi.circuit, y = fst_lin)) +
  geom_point(size = 1) +
  geom_abline(data = draws, aes(intercept = a, slope = b),
              color = "skyblue", size = 0.2, alpha = 0.25) +
  geom_abline(intercept = coef(m.evi)[1], slope = coef(m.evi)[2],
              color = "skyblue4", size = 1)+
  theme_classic()


m.fri<- stan_glm(formula=fst_lin ~  friction.circuit,
                 data = dd_analysis,
                 family = gaussian(link = "identity"),
                 prior = hs_prior, prior_intercept = t_prior,
                 seed = 121347, adapt_delta = 0.999,
                 chains = 4, iter = 4000)
draws <- as.data.frame(as.matrix(m.fri))
colnames(draws)[1:2] <- c("a", "b")
ggplot(dd_analysis, aes(x = friction.circuit, y = fst_lin)) +
  geom_point(size = 1) +
  geom_abline(data = draws, aes(intercept = a, slope = b),
              color = "skyblue", size = 0.2, alpha = 0.25) +
  geom_abline(intercept = coef(m.fri)[1], slope = coef(m.fri)[2],
              color = "skyblue4", size = 1)+
  theme_classic()




