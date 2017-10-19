# have a look at Alex's data of species coexistence

# clear workspace
mise::mise(vars = TRUE, console = TRUE, figs = TRUE, pkgs = TRUE)

# load packages ####
library(ggplot2)
library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(viridis)

# weekly transfers
# 10^-5 dilution of culture (25 µl)
# sp_1 - smooth, opaque yellow colonies
# sp_2 - smooth, translucent brown-cream colonies
# sp_3 - fuzzy, translucent yellow colonies
# sp_4 - smooth, opaque, white-cream colonies
# sp_5 - fuzzy/wrinkly, translucent pale brown colonies

# load in data
d <- read_excel('data/Community_coexistence_weekly_transfers.xlsx', range = 'A6:G54')
names(d) <- c('week', 'microcosm', 'sp_1', 'sp_2', 'sp_3', 'sp_4', 'sp_5')

d <- gather(d, 'sp', 'density', starts_with('sp_'))

d_mean <- group_by(d, sp, week) %>%
  summarise(., mean = mean(density, na.rm = TRUE)) %>%
  ungroup()

# plots
ggplot(d) +
  geom_line(aes(week, density, col = sp)) +
  theme_bw() +
  ggtitle('Dynamics of community coexistence through time per microcosm') +
  facet_wrap(~ microcosm) +
  scale_color_viridis(discrete = TRUE)

ggplot(d) +
  geom_line(aes(week, density, col = sp, group = interaction(microcosm, sp))) +
  theme_bw() +
  ggtitle('Dynamics of community coexistence through time per species') +
  facet_wrap(~ sp) +
  scale_color_viridis(discrete = TRUE)


ggplot(d) +
  geom_line(aes(week, density, col = sp, group = interaction(microcosm, sp)), alpha = 0.5) +
  geom_line(aes(week, mean, col = sp), size = 1.5, d_mean) +
  theme_bw() +
  ggtitle('Dynamics of community coexistence through time') +
  scale_color_viridis(discrete = TRUE)

