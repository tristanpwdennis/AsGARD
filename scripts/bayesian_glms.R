library(tidyverse)
library(caret)
library(GGally)
library(ggplot2)
library(marginaleffects)
library(corrplot)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(rstanarm)
options(mc.cores = 1)
library(loo)
library(projpred)
SEED=14124869



# Load data and scale it
# I modelled scaled and unscaled and it shouldn't make a difference
# But it will make it easier to compare effect sizes
# And when we are reporting final models we will backtranslate coefficients etc to be in og units.
env_data_unscaled <- read.csv("~/Projects/AsGARD/data/environmental_modelling/environmental_data_unscaled.csv")
all_vars <- c('fst_lin','pdist','mean.wind.traveltime','slope.circuit','friction.circuit','evi.circuit','hurs.circuit','ele.circuit','mdr.circuit','pop.circuit')
data_vars <- env_data_unscaled[all_vars]

# Calculate z scores for each variable
z_transformation <- function(col){
  r_mean <-mean(col, na.rm=TRUE)
  r_sd <- sd(col, na.rm=TRUE)
  r_z <- (col - r_mean) / r_sd
  return(r_z)
}

rescaled_data <- lapply(data_vars, z_transformation) # transform numeric vars
rescaled_data <- lapply(rescaled_data, as.numeric) # coerce all to numeric
rescaled_data$fst_lin <- env_data_unscaled$fst_lin
rescaled_data$invasion <- as.factor(env_data_unscaled$invasion)
rescaled_data <- data.frame(rescaled_data)

#load metadata
df_samples = read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_noqc.20241228.csv')

#Corr mat
tiff("~/Projects/AsGARD/figures/corrplot.tiff", units="in", width=5, height=5, res=300)
corrplot(cor(data_vars))
dev.off()

# convenience vars
p <- length(all_vars) -1 # Number of predictors in og model
t_prior.1.5 <- student_t(df = 7, location = 0, scale = 1.5) # the prior that gives the smallest 95% CI of posterior distributions
n <- nrow(rescaled_data) # Here we will use a regularised horseshoe prior with the assumption of multicolinearity and useless variables
p0 <- 3 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
rhs_prior <- hs(global_scale=tau0)



# Fit full bayesian model with rhs prior
post1_hs <- stan_glm(formula=fst_lin ~ pdist +mean.wind.traveltime+mdr.circuit + evi.circuit+pop.circuit+ele.circuit+hurs.circuit+friction.circuit, 
                  data = rescaled_data,
                  family = gaussian(),
                  prior = rhs_prior, QR=TRUE,
                  adapt_delta = 0.95,
                  seed = SEED, refresh=0, cores=6)

# Fit full bayesian model with t prior
post1_t <- stan_glm(formula=fst_lin ~ pdist +mean.wind.traveltime+mdr.circuit + evi.circuit+pop.circuit+ele.circuit+hurs.circuit+friction.circuit, 
                     data = rescaled_data,
                     family = gaussian(),
                     prior = t_prior.1.5, QR=TRUE,
                     adapt_delta = 0.95,
                     seed = SEED, refresh=0, cores=6)


(loohs <- loo(post1_hs, save_psis = TRUE,k_threshold = 0.7))
(loot<- loo(post1_t, save_psis = TRUE, k_threshold = 0.7))

loo_compare(loohs, loot) # the HS prior provides better model fit and predictive performance than the t prior 
summary(post1_hs, digits=6, probs = c(0.025, 0.975))
summary(post1_t, digits=6, probs = c(0.025, 0.975))

# take a look at posterior predictive checks
yrep <- posterior_predict(post1_hs, draws = 50)
ppc_dens_overlay(rescaled_data$fst_lin, yrep) # looks OK., model is smoothing to negative despite no negative values in data - model is not constrained to be positive but the error looks very small and we can continue
# take a look at posterior predictive checks
yrep <- posterior_predict(post1_t, draws = 50)
ppc_dens_overlay(rescaled_data$fst_lin, yrep) 

mcmc_areas(post1_t) # plot marginal posterior: we can see some likely colinearity (prob btw mdr and ele as per corrplot). we can see we definitely want to include friction and evi, but it's not clear which others
mcmc_areas(post1_hs) # plot marginal posterior: we can see some likely colinearity (prob btw mdr and ele as per corrplot). we can see we definitely want to include friction and evi, but it's not clear which others


#still wide CI for effect of sudan, let's try variable selection now but if this produces nonsense I will produce a more streamlined model first and then attempt
varsel_hs <- cv_varsel(post1_hs, method='forward', cv_method='loo', nclusters_pred=100, ndraws_pred=2000, nterms_max=15, parallel=TRUE) #run variable selection (this takes quite a long time)
saveRDS(varsel_hs,file='~/Projects/AsGARD/data/environmental_modelling/varsel_post1_hs.Rds') #save result
plot(varsel_hs, stats = c('rmse','elpd'))
#varsel_hs <- load('~/Projects/AsGARD/data/environmental_modelling/varsel_post1_hs.Rds')
#still wide CI for effect of sudan, let's try variable selection now but if this produces nonsense I will produce a more streamlined model first and then attempt
#varsel_t<- cv_varsel(post1_t, method='forward', cv_method='loo', nclusters_pred=100, ndraws_pred=2000, nterms_max=15, parallel=TRUE) #run variable selection (this takes quite a long time)
#saveRDS(varsel_hs,file='~/Projects/AsGARD/data/environmental_modelling/varsel_t.rds') #save result
#varselplot <- plot(varsel_t, stats = c('rmse','elpd'))
ggsave(plot=last_plot(), '~/Projects/AsGARD/figures/varselplot.svg',width = 10, height=6)

#let's select our veriables based on the RHS prior model
n_sel <- suggest_size(varsel_hs, stat = "rmse") #get suggested no vars

# Looks like we retain wind, evi, friction
#rank variable importance and plot
rk <- ranking(varsel_hs) 
( pr_rk <- cv_proportions(rk) )  
plot(pr_rk)
plot(cv_proportions(rk, cumulate = TRUE))
(predictors_final <- head(rk[["fulldata"]], n_sel+1) ) #get final predictors (eg 3 plus md because we want ti)
prj <- project(
  post1_hs,
  predictor_terms = predictors_final,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = FALSE
  ###
)

# Examine the down-projected model
prj_mat <- as.matrix(prj)
prj_drws <- as_draws_matrix(prj_mat)
prj_smmry <- posterior::summarize_draws(
  prj_drws,
  "median", "mad", function(x) quantile(x, probs = c(0.025, 0.975))
)

# Coerce to a `data.frame` because pkgdown versions > 1.6.1 don't print the
# tibble correctly:
prj_smmry <- as.data.frame(prj_smmry)
print(prj_smmry, digits = 1)

# Using the 'bayesplot' package, plot predicted values of fst vs actual
prj_predict <- proj_predict(prj)
ppc_dens_overlay(y = rescaled_data$fst_lin, yrep = prj_predict)
bayesplot_theme_set(ggplot2::theme_classic())
ggsave(plot = last_plot(), '~/Projects/AsGARD/figures/bayesplot.svg')



# Create final selected model object
post_proj <- stan_glm(formula=fst_lin ~ evi.circuit+friction.circuit+mdr.circuit, data = rescaled_data,
                  family = gaussian(), 
                  prior = rhs_prior, QR=FALSE,
                  cores = 6,
                  chains=4,
                  seed = SEED, refresh=0)
posterior_plot <- mcmc_areas(as.matrix(post_proj),
                             pars = names(fixef(post_proj)),
                             prob = 0.95) + 
ggtitle("Posterior Distributions with 95% Credible Intervals")
ggsave(plot=posterior_plot, filename='~/Projects/AsGARD/figures/postplot_finalmodel.svg', width = 6, height=6)

summary(post_proj, digits=3, probs = c(0.025, 0.975))
mcmc_areas(post_proj) # plot marginal posterior: we can see some likely colinearity (prob btw mdr and ele as per corrplot). we can see we definitely want to include friction and evi, but it's not clear which others
(loo_final <- loo(post_proj, save_psis = TRUE))

# Fit model of isolation by distance alone
post_ibd <- stan_glm(formula=fst_lin ~ pdist, 
                  data = rescaled_data,
                  family = gaussian(),
                  prior = t_prior.1.5,
                  seed = SEED, refresh=0, cores=6)
summary(post_ibd, digits=3, probs=c(0.025, 0.975))
(looibd <- loo(post_ibd, save_psis = TRUE))
mcmc_areas(post_ibd) # plot marginal posterior
yrep <- posterior_predict(post_ibd, draws = 50)
ppc_dens_overlay(rescaled_data$fst_lin, yrep) # looks ok too


# Compare models to see whether model with other predictors is a better fit than isolation by distance alone plus invasion

# Fit model to see if adding invasion as an interaction improves model fit
post_inv <- stan_glm(formula=fst_lin ~ invasion*(evi.circuit+friction.circuit+mdr.circuit), 
                      data = rescaled_data,
                     family = gaussian(),
                     prior = t_prior.1.5,
                     seed = SEED, refresh=0, cores=6)
summary(post_inv, digits=6, probs=c(0.025, 0.975))
(looinv <- loo(post_inv, save_psis = TRUE))
mcmc_areas(post_inv) # plot marginal posterior
yrep <- posterior_predict(post_inv, draws = 50)
ppc_dens_overlay(rescaled_data$fst_lin, yrep) # looks ok too


# Fit model to see if adding invasion as an interaction improves model fit
post_spat <- stan_glm(formula=fst_lin ~ pdist*(evi.circuit+friction.circuit+mdr.circuit), 
                     data = rescaled_data,
                     family = gaussian(),
                     prior = t_prior.1.5,
                     seed = SEED, refresh=0, cores=6)
summary(post_spat, digits=6, probs=c(0.025, 0.975))
(loospat <- loo(post_spat, save_psis = TRUE, k_threshold=0.7))
mcmc_areas(post_spat) # plot marginal posterior
yrep <- posterior_predict(post_spat, draws = 50)
ppc_dens_overlay(rescaled_data$fst_lin, yrep) # looks ok too


loo_compare(looibd, loo_final, looinv, loospat) #our submodel is a better fit than isolation by distance by itself, and when taking invasion into account




####
# Now we want to take a look at our model predictions while w ehold the other variables constant (10th and 90th percentile, median)
# plotting them against thge original, unscaled data
####



plot_predictions(post_proj)

plot_predictions(post_proj, condition = list("fst_lin", "friction.circuit" = "threenum"))

?plot_predictions




#####
# Finally, let's plot the rasters for the three retained predictors
#####

library(windscape)
library(terra)
points <- df_samples %>% filter(country %in% c('Sudan','Ethiopia','Djibouti')) %>%  dplyr::select(country, location,pop_code, latitude, longitude) %>% unique()

# Collect three plots
raster_plot_list = list()
vars <- c('mdr','evi','friction')
for (i in seq_along(vars)){
  var <- vars[[i]]
  raster <- rast(paste0('~/Projects/AsGARD/data/environmental_modelling/circuitscape_data/',var,'/',var,'.asc'))
  r_df <- as.data.frame(raster, xy = TRUE)
  colnames(r_df) <- c('x','y','layer')
  plot <- ggplot(r_df, aes(x=x,y=y, fill=layer))+
    geom_raster()+
    geom_point(data = points, aes(x = longitude, y = latitude), inherit.aes = FALSE, colour='darkgrey', alpha=0.6, size=3) +
    scale_fill_distiller(palette = "Blues") +
    theme_classic()+
    coord_sf(xlim=c(27, 47),ylim=c(6,21), expand=FALSE)+
    labs(fill=paste0('Scaled ',var))+
    borders("world", colour = "black", fill = NA)+  # Add world borders
    theme(#legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12))
  ggsave(filename=paste0('~/Projects/AsGARD/figures/',var,'.svg'), plot, width=5, height=5, units = 'in')
  
  raster_plot_list[[i]] <- plot
}


# Load necessary libraries
library(ggplot2)
library(dplyr)
install.packages('jtools')
  