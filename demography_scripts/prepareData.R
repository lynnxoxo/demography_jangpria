library(readxl)
library(readr)
library(tidyverse)
library(dplyr)
library(magrittr)
require(reshape2)
library(forcats)

cemdat_wide <- read_excel("jang_data/cemeteries_JANGPRIA_large_data.xlsm")

cemdat_long <- prepare_mortuary_data(cemdat_wide, site_name = "cemetery_name", start_col =  24, end_col =  359, category_order = "period, age, sex_gender", additional_cols = "relative_chronology_system, cemetery_start_local, cemetery_end_local, cemetery_start_lt, cemetery_end_lt, cemetery_size, cemetery_anthro, cemetery_grav_dat, dig_date, cemetery_quality", filter_conditions = "value > 0")
