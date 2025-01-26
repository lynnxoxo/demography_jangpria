library(ggplot2)
library(tidyverse)
require(reshape2)
library(forcats)
library(paletteer)
library(readr)

#colour
scale_colour_paletteer_d("MetBrewer::Greek")
scale_color_paletteer_d("MetBrewer::Greek")
scale_fill_paletteer_d("MetBrewer::Greek")
paletteer_d("MetBrewer::Greek")

#coords
coord_cartesian(ylim=c(0,100))

#data
term_frequencies_HoR <- read_csv("./jang_data/term_frequencies_HoR_BA_en1870_2024_tidy.csv")
publdata <- term_frequencies_HoR

#tidy long
long_publdata <- publdata %>%
  pivot_longer(cols = !starts_with("w"), names_to = "decade", values_to = "occurence")
  

#bar geom_text(aes(label=occurence), position=position_dodge(width=0.9), vjust=-0.25, size = 1, accuracy = 0.01) +
ggplot(long_publdata, aes(x = decade, y = occurence, fill=word)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  facet_wrap(vars(word), strip.position = "bottom", ncol = 4) +
  scale_fill_paletteer_d("MetBrewer::Greek") +
  theme(axis.title.y=element_blank(), axis.title.x = ) +
  labs (fill = "", x= "", y = "", title = "Research trends in archaeological literature", caption = "Dataset built with constellate.org, based on JSTOR indexed publications filtered by keyword 'archaeology' between 1870 and 2024")+
  guides(fill = "none")+
  theme_minimal()

