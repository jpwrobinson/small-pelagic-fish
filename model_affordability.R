## R script to fit Bayesian model to catch-price IHH dataset and model affordability of species groups in wild-caught fisheries
## James PW Robinson, Lancaster University, 14-October-2022

## package load
library(tidyverse)
library(rethinking)

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
  
## read price/catch dataset
price<-read.csv(file = 'price_minimum_dataset.csv')

with(price, table(fg)) 

## plot response histogram
hist(price$cost_100g_relative_to_cheapest_daily_diet)
hist(log(price$cost_100g_relative_to_cheapest_daily_diet))

## log price to NA
price$value<-price$cost_100g_relative_to_cheapest_daily_diet
price$logvalue<-log(price$cost_100g_relative_to_cheapest_daily_diet)

# log catch and lmax
price$log_catch<-log(price$catch_avg_annual_tonnes)
price$log_lmat<-log(price$Lmat)

## scale cont. variables
price$log_catch_scaled<-scale(price$log_catch)[,1]
price$log_lmat_scaled<-scale(price$log_lmat)[,1]

## convert data to list, characters to factors for rethinking
price<-droplevels(price %>% ungroup() %>% mutate_if(is.character, as.factor))
focal.list<-as.list(price)

# add list of region ID for each country name
countryLookup<-price %>% distinct(country, subregion)
focal.list$countryLookup <- countryLookup %>% pull(subregion) %>% as.factor()

subregionLookup <-price %>% distinct(subregion, region) 
focal.list$subregionLookup <- subregionLookup %>% pull(region) %>% as.factor()

fgLookup<-price %>% distinct(species, species_functional_group) 
focal.list$fgLookup<-fgLookup %>% pull(species_functional_group) %>% as.factor()

print(paste(
  'Model fitted to:', 
  length(price$catch_avg_annual_tonnes), 'catch records', 
  round(sum(price$catch_avg_annual_tonnes),0), 'tonnes', 
  length(unique(price$country)), 'countries',
  length(unique(price$species)), 'species'))

## run mode
m<-ulam(
  alist(
    # response variable
    cost_100g_relative_to_cheapest_daily_diet ~ dlnorm(mu, sigma),
    
    ## observation level model (fisheries)
    mu <- #beta_fishery[country_sector] +  
        beta_country[country] +
        beta_sp[species]*sigma_sp + beta_fg[species_functional_group]*sigma_fg +
        b_lmat*log_lmat_scaled + b_catch*log_catch_scaled,
  
    ## country-level intercepts
    beta_country[country] ~ dnorm(mu_country, sigma_country),
    ## countries nested in subregions
    mu_country <- beta_subregion[countryLookup],
    beta_subregion[countryLookup] ~ dnorm(mu_subregion, sigma_subregion),
    
    ## region drawn from global mean
    mu_subregion <- beta_region[subregionLookup],
    beta_region[subregionLookup] ~ dnorm(mu_region, sigma_region),

    mu_region <- y,
    y ~ dnorm(3, 10),
  
    beta_sp[species] ~ dnorm(0, 1),
    beta_fg[species_functional_group] ~ dnorm(0, 1),

    b_lmat ~ dnorm(0, 1),
    b_catch ~ dnorm(0, 1),
    
    c(sigma, sigma_sp, sigma_fg) ~ dexp(1),
    c(sigma_subregion, sigma_country) ~ dunif(0, 20),
    c(sigma_region) ~ dunif(0, 10)
  ),
  data=focal.list, iter=5000, warmup=1500, 
  chains=3, cores = 6, 
  control=list(adapt_delta=0.99), log_lik = TRUE)


# inspect output
precis(m,3) %>% data.frame() %>% filter(Rhat4 > 1.01) 
precis(m)
dashboard(m)
ll<-link(m)
m_mu<-apply(ll$mu, 2, median)
plot(m_mu, focal.list$logvalue)
abline(0, 1)


# posterior predictions for each sector + fg
## get fg preds
tops<-price 
tops$species_functional_group<-factor(tops$species_functional_group)
tops$species<-factor(tops$species)
tops$country <- factor(tops$country)
tops$subregion <- factor(tops$subregion)
tops$region <- factor(tops$region)
tops$countryLookup <- countryLookup$subregion[match(tops$country, countryLookup$country)]
tops$subregionLookup <- subregionLookup$region[match(tops$subregion, subregionLookup$subregion)]

# tops$log_catch_scaled = 0
# tops$log_lmax_scaled = 0

# get posterior samples of mu
post = extract.samples(m)
meds<-link(m, data = tops, post = post)$mu

## median and CIs
tops$mu<-apply(meds, 2, median)
tops$lower95 <- apply( meds , 2 , HPDI , prob=0.95 )[1,]
tops$upper95 <- apply( meds , 2 , HPDI , prob=0.95 )[2,]
tops$lower50 <- apply( meds , 2 , HPDI , prob=0.50 )[1,]
tops$upper50 <- apply( meds , 2 , HPDI , prob=0.50 )[2,]
      
## estimate metrics
tops <- tops %>% 
      rename_with(~tolower(str_replace_all(.x, '_mu', ''))) %>%
      rename(omega_3 = omega3, vitamin_a = vitamina) %>% 
      ungroup() %>% 
      rowwise() %>% mutate(n_targets = sum(c_across(ca_rda:om_rda) > 10)) %>%  ## estimate nutrition targets (10% RDA) for each species)
      mutate(
            cost_caloric_adequacy = exp(mu), # predicted cost of 100g portion, relative to starchy staples (%)
            nutrient_adequacy = rowMeans(across(c(ca_rda, fe_rda, se_rda, zn_rda, vita_rda, om_rda))), # nutrient adequacy = mean RDA contributions
            
            portion_adq33 = 33 / nut_adq * 100, # = portion size to get 33% adequacy
            usd_adq33 = exp(mu) * (portion_adq33/100) ## the cost to reach 33% nutrient adequacy, relative to starchy staples
            )


write.csv(tops, file = 'predicted_species_nutrient_affordability.csv', row.names=FALSE)