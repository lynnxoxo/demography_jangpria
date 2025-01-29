## This is the primary file which subsequently runs all further analysis
## decentralised in other r scripts

# Load tidyverse
library(tidyverse) # to pipe (%>%) and map across each file
library(purrr) # for mapping across files
library(renv)
library(tibble)
library(dplyr)
library(readr)
library(ggplot2)
library(flextable)
library(officer)
library(mortAAR)
library(ggthemes)

renv::restore()
renv::snapshot()
dependencies_list <- renv::dependencies()
renv::record(as.list(dependencies_list$Package))

options(readr.show_col_types = FALSE)


# List files and source each
#list.files("./demography_scripts/", full.names = TRUE) %>% map(source)

source("./demography_scripts/functions_MortDem.R")

##data modeling start##


# Load burial and harmonization data
jangdat_indiv <- read_csv("jang_data/csv_cemeteries_data_indiv.csv")
jangdat_aggreg <- read_csv("jang_data/csv_cemeteries_data_aggreg.csv")
jangdat_meta <- read_csv("jang_data/csv_cemeteries_data_meta.csv",
                         col_types = cols_only(
                           site_name = col_character(),
                           coord_y = col_double(),
                           coord_x = col_double(),
                           cemetery_country = col_character(),
                           cemetery_region = col_character(),
                           cemetery_admin = col_character(),
                           relative_chronology_system = col_character(),
                           cemetery_start_ia = col_character(),
                           cemetery_end_ia = col_character(),
                           cemetery_size = col_integer(),
                           cemetery_anthro = col_character(),
                           cemetery_grav_dat = col_character(),
                           dig_date = col_character(),
                           cemetery_quality = col_character()
                         ))

harmonization_data <- read_csv("jang_data/harmonization_jangpria.csv")


# Create harmonization table from the raw harmonization data
harmonization_table <- create_harmonization_table(harmonization_data, 
                                                  system_name_col = "system_name", 
                                                  relative_chron_col = "relative_chron", 
                                                  horizon_start_col = "horizon_start", 
                                                  horizon_end_col = "horizon_end", 
                                                  fade_in_start_col = "fade_in_start", 
                                                  fade_out_end_col =  "fade_out_end")

  # Create a plot of the harmonization table
harmonization_plot <- plot_harmonization_table(harmonization_table,  # Corrected to use harmonization_table
                                                 system_name_col = "system_name", 
                                                 relative_chron_col = "relative_chron", 
                                                 horizon_start_col = "horizon_start", 
                                                 horizon_end_col = "horizon_end", 
                                                 fade_in_start_col = "fade_in_start", 
                                                 fade_out_end_col =  "fade_out_end")


# Prepare individual mortuary data using the harmonization table
prepared_indiv_data <- prepare_mortuary_data_indiv(jangdat_indiv, 
                                                   harmonization_table = harmonization_table,
                                                   site_name = "site_name", 
                                                   system_name = "system_name", 
                                                   id_burial_col = "id_burial",  # Corrected argument names
                                                   period_col = "period", 
                                                   age_col = "age", 
                                                   sex_gender_col = "sex_gender")

prepared_aggreg_data <- prepare_mortuary_data_aggreg(
  data = jangdat_aggreg,
  harmonization_table = harmonization_table,
  site_name = "site_name",
  system_name = "system_name",
  period = "period",
  amount = "amount",
  sex_gender = "sex_gender",
  age = "age"
)

harmonized_data <- harmonize_chronologies(
  prepared_indiv_data, prepared_aggreg_data,
  harmonization_table = harmonization_table,
  site_name_col = "site_name", 
  period_col = "period", 
  horizon_start_col = "horizon_start", 
  horizon_end_col = "horizon_end", 
  harmonized_system_name_col = "system_name", 
  harmonized_relative_chron_col = "relative_chron", 
  horizon_bracket_size = 10)


## plotting start ##

#print(harmonization_plot)

results <- analyze_occupancy(
  harmonized_data, 
  jangdat_meta,
  site_name_col = "site_name", 
  horizon_bucket_col = "horizon_bucket",
  cemetery_size_col = "cemetery_size",
  age_col = "age",
  sex_col = "sex_gender",
  cemetery_region_col = "cemetery_region",
  coord_x_col = "coord_x",
  coord_y_col = "coord_y",
  fuzzy_match = TRUE,        # Enable fuzzy matching if needed
  max_fuzzy_distance = 0.2,  # Threshold for fuzzy matching
  change_threshold = 2,      # Minimum cemeteries for simultaneous change
  distance_threshold = 50000, # Distance in meters for spatial neighbors
  #horizon_bracket_size = 40,
  beast_min_length = 3 # Minimum length of TS for Rbeast
)

##print the scatterplots and density estimation
print(results[[1]])
print(results[[2]])
print(results[[3]])
print(results[[4]])
print(results[[5]])

p_animation <- results$p_animation

library(gapminder)
library(gganimate)
library(gifski)

hoagiffed <- animate(p_animation, nframes = 480, fps = 24, height=750, width=650, renderer=gifski_renderer(), end_pause = 24, start_pause = 24)

anim_save("hoangif.gif", hoagiffed)





## mortAAR start ##

# Replace existing age categories in `harmonized_data` with numeric ranges
jang_mortTable <- harmonized_data %>%
  mutate(age = case_when(
    age == "inf1" ~ "0-6",
    age == "inf2" ~ "7-14",
    age == "juv" ~ "15-20",
    age == "fad" ~ "21-29",
    age == "ad" ~ "30-40",
    age == "mat" ~ "41-60",
    age == "sen" ~ "61-85",
    TRUE ~ age # Keep the age as is if it's not one of these labels
  ))

# Sort specific ages into appropriate age brackets
jang_mortTable <- jang_mortTable %>%
  mutate(age = case_when(
    as.numeric(age) >= 0 & as.numeric(age) <= 6 ~ "0-6",
    as.numeric(age) >= 7 & as.numeric(age) <= 14 ~ "7-14",
    as.numeric(age) >= 15 & as.numeric(age) <= 20 ~ "15-20",
    as.numeric(age) >= 21 & as.numeric(age) <= 29 ~ "21-29",
    as.numeric(age) >= 30 & as.numeric(age) <= 40 ~ "30-40",
    as.numeric(age) >= 41 & as.numeric(age) <= 60 ~ "41-60",
    as.numeric(age) >= 61 & as.numeric(age) <= 85 ~ "61-85",
    TRUE ~ age # Leave age unchanged if it doesn't fit into any of these ranges
  ))

jang_mortTable <- jang_mortTable %>%
  dplyr::filter(!is.na(age))

jang_mortTable <- jang_mortTable %>%
  tidyr::separate(age, c("from", "to"))%>%
  transform(from = as.numeric(from),to = as.numeric(to))

td_prepared <- prep.life.table(
  jang_mortTable, 
  dec = NA, 
  agebeg = "from",
  ageend = "to", 
  group = "site_name", 
  method = "Standard",
  agerange = "included"
)

td_result <- td_prepared %>% 
  life.table()


##plot demographic data
td_result %>% plot(display = c("qx", "dx", "lx", "ex", "rel_popx"), line_vis = "color", prefer.ggplot = TRUE) 


### Site specific data like MI MMR F_M_R 
masc_ind_mutate <- na.omit(jang_mortTable$sex_gender, jang_mortTable$from, jang_mortTable$to) 
  comb_life_prep <- filter(jang_mortTable) 
  camm_life_prep <- filter(jang_mortTable, site_name == "Cammer") 
  kolb_life_prep <- filter(jang_mortTable, site_name == "Kolbow") 
  much_life_prep <- filter(jang_mortTable, site_name == "Muchow") 
  nich_life_prep <- filter(jang_mortTable, site_name == "Nichel") 
  sod_life_prep <- filter(jang_mortTable, site_name == "Soderstorf") 
  soh_life_prep <- filter(jang_mortTable, site_name == "Søhale") 
  soer_life_prep <- filter(jang_mortTable, site_name == "Sörup II")
  veld_life_prep <- filter(jang_mortTable, site_name == "Veldbækvej III")
  
  
masc_ind_put <- prep.life.table(
  put_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_camm <- prep.life.table(
  camm_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_kolb <- prep.life.table(
  kolb_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_much <- prep.life.table(
  much_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_nich <- prep.life.table(
  nich_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_sod <- prep.life.table(
  sod_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_soh <- prep.life.table(
  soh_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_soer <- prep.life.table(
  soer_life_prep,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to",
) 
masc_ind_comb <- prep.life.table(
  masc_ind_mutate,
  group ="sex_gender",
  agebeg = "from",
  ageend = "to"
) 

jang_life <- life.table(masc_ind_comb)

lt.sexrelation(jang_life$f, jang_life$indet)
  

lt.representativity(jang_life)


