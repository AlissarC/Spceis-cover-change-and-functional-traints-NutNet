
library(dplyr)   # Data manipulation: filtering, grouping, summarizing,
# creating variables, and joining datasets.

library(tidyr)   # Data cleaning: replacing missing values (NA)
# and organizing datasets.

library(stringr) # String manipulation: standardizing species names,
# converting text to lowercase and sentence case.

################ load cover data ###################################
Cover_data <- read.csv("./data/full-cover_2025-01-31_v2.csv")
names(Cover_data)
nrow(Cover_data) #351776

######### Select only data of plant species (live =1)################
Cover_data_live = subset(Cover_data, live==1)
nrow(Cover_data_live) #312995
table(Cover_data_live$live)

######### remove data of Bryophytes and Lichens and unknown species (Family = NULL)################
Cover_data_live_clean = subset(Cover_data_live, Family!="NULL")
nrow(Cover_data_live_clean)#307205

Cover_data_live_clean$species <- str_to_sentence(str_to_lower(Cover_data_live_clean$species)) ## to be coherent with standard species format

names(Cover_data_live_clean)

### average for each species per plot the max_cover of subplots (which are similar)
Cover_by_plot <- Cover_data_live_clean %>%
  group_by(Year, site_name, site_code, block, plot,trt, year_trt, species) %>%
  summarise(
    max_cover = mean(max_cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::select(Year, site_name, site_code, block, plot, trt,year_trt, species, max_cover)

names(Cover_by_plot)
nrow(Cover_by_plot) # 306851

########### select only sites with 3 or more years of measurements 
Cover_by_plot_select <- Cover_by_plot %>%
  group_by(site_code) %>%
  filter(n_distinct(Year) >= 3) %>%
  ungroup()
n_distinct(Cover_by_plot_select$site_code) # 95 sites

names(Cover_by_plot_select)

Cover_by_plot_select_nYears<- Cover_by_plot_select %>%
  group_by(site_code) %>%
  mutate(
    n_years = n_distinct(Year[year_trt >= 1])
  ) %>%
  ungroup()
names(Cover_by_plot_select_nYears)

############## select only measurements post treatment 
Cover_post <- Cover_by_plot_select_nYears %>%
  filter(year_trt >= 1)
n_distinct(Cover_post$site_code) # 95 sites
names(Cover_post)

################### calculate lambda  ############################################
# Define the function to compute lambda per plot group: if a species in a plot disappear some years, 
#set its max cover at 0
# calculate also the relative abundance, and lambda relative abundance 
compute_lambda_per_plot <- function(plot_df, keys, c = 0.001) {
  all_years <- sort(unique(plot_df$Year))
  all_species <- unique(plot_df$species)
  
  full_grid <- expand.grid(
    Year = all_years,
    species = all_species
  )
  
  full_plot <- full_grid %>%
    dplyr::left_join(
      plot_df %>% dplyr::select(Year, species, max_cover),
      by = c("Year", "species")
    ) %>%
    dplyr::mutate(
      max_cover = tidyr::replace_na(max_cover, 0)
    )
  
  full_plot <- full_plot %>%
    dplyr::group_by(Year) %>%
    dplyr::mutate(
      total_cover = sum(max_cover, na.rm = TRUE),
      Ci = dplyr::if_else(total_cover > 0, max_cover / total_cover, 0),
      one_minus_Ci = 1 - Ci
    ) %>%
    dplyr::ungroup()
  
  full_plot <- full_plot %>%
    dplyr::arrange(species, Year) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      lambda_cover = (max_cover + c) / (dplyr::lag(max_cover) + c),
      log_lambda_cover = log(lambda_cover),
      
      lambda_Ci = (Ci + c) / (dplyr::lag(Ci) + c),
      log_lambda_Ci = log(lambda_Ci),
      
      lambda_cover = dplyr::if_else(dplyr::row_number() == 1, NA_real_, lambda_cover),
      log_lambda_cover = dplyr::if_else(dplyr::row_number() == 1, NA_real_, log_lambda_cover),
      lambda_Ci = dplyr::if_else(dplyr::row_number() == 1, NA_real_, lambda_Ci),
      log_lambda_Ci = dplyr::if_else(dplyr::row_number() == 1, NA_real_, log_lambda_Ci)
    ) %>%
    dplyr::ungroup()
  
  return(full_plot)
}

lambda_cover_NutNet <- Cover_post %>%
  dplyr::group_by(site_name, site_code, block, plot,trt,n_years) %>%
  dplyr::group_modify(~ compute_lambda_per_plot(.x, .y, c = 0.001)) %>%
  dplyr::ungroup()
names(lambda_cover_NutNet)
nrow(lambda_cover_NutNet) # 534669
n_distinct(lambda_cover_NutNet$species) # 2945

############################ Average across all years #####################################

lambda_avg_allYears <- lambda_cover_NutNet %>%
  dplyr::group_by(
    site_name, site_code, n_years,block, plot, trt, species
  ) %>%
  dplyr::summarise(
    n_all_transitions = sum(!is.na(lambda_cover)),
    
    mean_lambda_cover_all = mean(lambda_cover, na.rm = TRUE),
    mean_log_lambda_cover_all = mean(log_lambda_cover, na.rm = TRUE),
    mean_lambda_Ci_all = mean(lambda_Ci, na.rm = TRUE),
    mean_log_lambda_Ci_all = mean(log_lambda_Ci, na.rm = TRUE),
    
    .groups = "drop"
  )
names(lambda_avg_allYears)
nrow(lambda_avg_allYears) # 63394


n_distinct(lambda_avg_allYears$species) # 2945


##################### Combine with climate ###################
sites <- read.csv("./data/sites_coord.csv")
names(sites)
length(unique(sites$site_code)) # 175

lambda_Avg_coord <- lambda_avg_allYears %>%
  left_join(sites, by = c("site_code", "site_name"))
names(lambda_Avg_coord) 
nrow(lambda_Avg_coord) # 63394

climate_summary <- read.csv("./data/climate.csv")
names(climate_summary)

climate_summary_select <- climate_summary %>% 
  dplyr::select(
    -continent, -country, -habitat,
    -first_year, -final_year, -active,
    -number_experiment_years, -experiment_type)
names(climate_summary_select)
names(lambda_Avg_coord)

lambda_Avg_climate <- lambda_Avg_coord %>%
  left_join(climate_summary_select, by = c("site_name", "site_code", "latitude", "longitude"))
names(lambda_Avg_climate)
nrow(lambda_Avg_climate) # 63394

n_distinct(lambda_Avg_climate$species) # 2930

lambda_Avg_climate <- lambda_Avg_climate %>% filter(!is.na(aridity))
nrow(lambda_Avg_climate) #  63300

names(lambda_Avg_climate)

############# merge with all traits  ############################### 
length(unique(lambda_Avg_climate$species)) # 2930

All_traits_merged <- read.csv("./data/All_traits_merged.csv")
names(All_traits_merged)

lambda_cover_NutNet_traits <- lambda_Avg_climate %>%
  left_join(All_traits_merged, by = c("species"))
names(lambda_cover_NutNet_traits) 
nrow(lambda_cover_NutNet_traits) # 63300


lambda <- lambda_cover_NutNet_traits 

lambda$trt

lambda$Ntrt = 0 ## creation of a column with Ntrt= 0 
lambda$Ntrt[lambda$trt == 'N' | lambda$trt == 'NP' | 
              lambda$trt == 'NK' | lambda$trt == 'NPK' |
              lambda$trt == 'NPK+Fence'] = 1 ## set Ntrt=1 if any point receives N

### repeat the same for the other treatments 
lambda$Ptrt = 0
lambda$Ptrt[lambda$trt == 'P' | lambda$trt == 'NP' |
              lambda$trt == 'PK' | lambda$trt == 'NPK' | 
              lambda$trt == 'NPK+Fence'] = 1

lambda$Ktrt = 0
lambda$Ktrt[lambda$trt == 'K' | lambda$trt == 'NK' | 
              lambda$trt == 'PK' | lambda$trt == 'NPK' | 
              lambda$trt == 'NPK+Fence'] = 1


lambda <- lambda %>%
  mutate(
    site_code = as.factor(site_code),
    block = as.factor(block),
    plot = as.factor(plot),     
    trt = as.factor(trt)         
  )

names(lambda)
nrow(lambda) # 63300

lambdanofence <- lambda %>%
  filter(!trt %in% c("Other", "Fence", "NPK+Fence"))
length(unique(lambdanofence$species[!is.na(lambdanofence$Nmass_NutNet)])) # 1723 
table(lambdanofence$trt)


lambdanofence_gram_forb<- lambdanofence %>%
  filter(growth_form %in% c("forb", "graminoid"))
length(unique(lambdanofence_gram_forb$species[!is.na(lambdanofence_gram_forb$Nmass_NutNet)])) # 1545


table(lambdanofence_gram_forb$growth_form)
table(lambdanofence_gram_forb$photosynthetic_pathway)

#Keep only well-known and determined C3 vs C4
final_df <- lambdanofence_gram_forb %>%
  filter(photosynthetic_pathway %in% c("C3", "C4"))
length(unique(final_df$species[!is.na(final_df$Nmass_NutNet)])) # 1514

table(final_df$photosynthetic_pathway)
table(final_df$growth_form)

write.csv(final_df, "./output/final_df.csv", row.names=FALSE)
nrow(final_df) # 44419
names(final_df)
n_distinct(final_df$site_code) #94

