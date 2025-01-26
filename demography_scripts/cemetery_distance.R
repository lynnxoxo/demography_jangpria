library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)

distanceData <- read_csv("jang_data/west_jutland_cemetery_distances.csv") %>% 
  pivot_longer(
    cols = !starts_with("Site"),  # All columns that represent distances
    names_to = "Feature_Type",
    values_to = "Distance"
  )

ggplot(distanceData, aes(x = Feature_Type, y = Distance)) +
  geom_boxplot() +
  labs(
  #  x = "Feature Type",
    y = "Distance of cemeteries to... in km"
  ) +
  theme_minimal()+
  theme(
    axis.title.x = element_blank()
  ) 


