# small-pelagic-fish
Data and code accompanying Robinson et al. (in press). 

* catch_minimum_dataset.csv is average annual catch of species in 39 countries
* price_minimum_dataset.csv is average annual catch and ex-vessel price of species in 39 countries, standardised by relative food prices (caloric adequacy)
* model_affordability.R fits Bayesian model to price_minimum_dataset.csv and samples posterior for predictions, generating:
* predicted_species_nutrient_affordability.csv is predicted affordability of each species group in each country
