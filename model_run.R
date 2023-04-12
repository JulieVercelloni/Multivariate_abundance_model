# Model associated to the paper Kennedy et al. 2020: Coral Reef Community Changes in Karimunjawa National Park, Indonesia: Assessing the Efficacy 
#of Management in the Face of Local and Global Stressors published in Journal of Marine Science and Engineering 
# Code written by Julie Vercelloni
# Contact: j.vercelloni@aims.gov.au 

rm(list=ls())

source("R/packages.R")
source("R/functions.R")

nreef = 18 # 18 reefs
nyear = 2 # 2 years of observations 
ngroup = 5 # 5 benthic groups

# Create synthetic data with proportions of "ngroup" benthic groups across "nyear" year of observations and "nreef" reefs. 
make_data(nreef,nyear)

# Add random management zones associated to the reef (6 zones) 

abundance_table <- abundance_table %>%
  mutate(Zone = rep(rep(1:6,each = 3),2)) %>%
  mutate(DHW = c(rtruncnorm(n = 18,a = 1,b = 8,mean = 2,sd = 1),
                 rtruncnorm(n = 18,a = 1,b = 8,mean = 6,sd = 4)))
# Model prep. 

# response variables 
y<-ceiling(abundance_table[,1:5]*50)%>%as.matrix() # transform proportions into counts as per the study methodology to run the negative binomial model 

# covariates
X <- abundance_table %>% 
  dplyr::select(year,Zone, DHW)
X$year <- as.numeric(X$year)

# MCMC control

# for testing 
niter = 900
nburn = 300
nthin = 10

# for running final model 
#niter = 9000
#nburn = 3000
#nthin = 10

mcmc_control <- list(n.burnin = nburn, n.iteration = niter,
                    n.thin = nthin, n.chains = 3)

# Model 
mod_nb <- boral(y,X=X, family = "negative.binomial", row.eff = "random",
                mcmc.control = mcmc_control,save.model=TRUE,
                lv.control = list(num.lv = 2, type = "independent", distmat = NULL),
                prior.control = list(type = c("normal","normal","normal","uniform"),
                hypparams = c(10, 10, 10, 10)))

plot(mod_nb)

########## Extract posterior distributions for responses 

group<-rownames(mod_nb$X.coefs.mean)

mcmc.raw<-get.mcmcsamples(mod_nb)

mcmc.nb<-coda_df(mcmc.raw)%>%
  tidyr::gather(key=Parameter,value = value) 

nrow <- (niter - nburn)/nthin

resp <- mcmc.nb %>%
  filter(str_detect(Parameter, 'X.coefs'))%>%
  mutate(group = rep(rep(group,each=nrow),3))%>%
  mutate(niter = rep(1:nrow,3*ngroup))%>%
  mutate(cov = rep(colnames(X), each =ngroup*nrow))

# Visualization
ggplot(resp, aes(y=group, x=value)) +
  facet_wrap(~cov) +
  stat_halfeye(.width = c(.90, .5)) +
  geom_vline(xintercept = 0, linetype = "dashed") 

# Summary table 
summary <- resp %>% 
  group_by(group,cov)%>%
  median_hdci(value) %>%
  mutate(Sig = ifelse(.lower<0 & .upper<=0 | .lower>=0 & .upper>0,1,0))

########## Latent variables
lv1 <- grep("^lv.c.+1]$", unique(mcmc.nb$Parameter), value=TRUE)

mcmc.lv1<-mcmc.nb%>%filter(Parameter%in%lv1)%>%
  mutate(group=rep(group,each=dim(mcmc.raw)[[1]]))%>%
  mutate(param = "lv1")

lv2 <- grep("^lv.c.+2]$", unique(mcmc.nb$Parameter), value=TRUE)

mcmc.lv2<-mcmc.nb%>%filter(Parameter%in%lv2)%>%
  mutate(group=rep(group,each=dim(mcmc.raw)[[1]]))%>%
  mutate(param = "lv2")

# Dispersion 
disp_names <- grep("^lv.c.+4]$", unique(mcmc.nb$Parameter), value=TRUE)

mcmc.disp<-mcmc.nb%>%filter(Parameter%in%disp_names)%>%
  mutate(group=rep(group,each=dim(mcmc.raw)[[1]]))%>%
  mutate(param = "dispersion")

# Summarize latent variables and dispersion parameter 

sum_table <- rbind(mcmc.lv1, mcmc.lv2,mcmc.disp)%>%group_by(group,param)%>%
  median_hdci(value) %>% 
  mutate(Sig = ifelse(.lower<0 & .upper<=0 | .lower>=0 & .upper>0,1,0))

# Covariance associated to the latent variables 
res.cor<-get.residual.cor(mod_nb, est="mean")

corrplot(res.cor$sig.cor, type = "lower",diag = FALSE,
         title = "Residual correlations",
         mar = c(3,0.5,2,1), tl.srt = 45)

# Covariance associated to the covariates
env.cor<-get.enviro.cor(mod_nb, est="mean")

corrplot(env.cor$sig.cor, type = "lower",diag = FALSE,
title = "Correlations due to covariates",
mar = c(3,0.5,2,1), tl.srt = 45)

# Residual ordination biplot - grouping by reef (change L82 in "function.R" if other grouping)
ordination_plot(mod_nb)
