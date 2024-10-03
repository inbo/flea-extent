library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsankey)
library(terra)
git_root <- rprojroot::find_root(rprojroot::is_git_root)
flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root)


temporal_stratification <- rast(file.path(flea_data, "data/2013_2016_2019", "temporal_stratification.tif"))

changes_df <- cats(temporal_stratification)[[1]] |> as_tibble()

(max_possible <- 9^3)
nrow(changes_df)

summary(changes_df$count)

# rastercel = 1 are (10 m x 10 m)
# number of transition types that represent less than 100 are = 1 ha
sum(changes_df$count < 100) # 1ha
sum(changes_df$count < 1000) # 10ha
sum(changes_df$count < 10000) # 100ha

# proportion of all cells that changed?
changes_df %>%
  group_by(stable) %>%
  summarize(
    pixelcount = sum(count)) %>%
  mutate(
    proportion = pixelcount / sum(pixelcount)
  )

# change cats per land use
# incl stable
changes_df %>%
  select(count, contains("changecat")) %>%
  pivot_longer(cols = contains("changecat")) %>%
  group_by(name, value) %>%
  summarize(
    pixelcount = sum(count)) %>%
  mutate(
    proportion = pixelcount / sum(pixelcount),
    name = reorder(name, pixelcount)) %>%
  ggplot() +
  geom_bar(aes(x = name, weight = pixelcount, fill = value)) +
  coord_flip()
# excl stable
changes_df %>%
  select(count, contains("changecat")) %>%
  pivot_longer(cols = contains("changecat")) %>%
  group_by(name, value) %>%
  filter(!grepl("^Stable", value)) %>%
  summarize(
    pixelcount = sum(count)) %>%
  mutate(
    proportion = pixelcount / sum(pixelcount),
    name = reorder(name, pixelcount),
    pixelcount2 = if_else(
      grepl("gain", value, ignore.case = TRUE),
      pixelcount, -pixelcount)) %>%
  ggplot() +
  geom_bar(aes(x = name, weight = pixelcount2, fill = value)) +
  coord_flip()

# all transitions

df <- changes_df %>%
  ggsankey::make_long(
    lg2013_label, lg2016_label, lg2019_label,
    value = count
  )

df2 <-  df %>%
  group_by(x, node) %>%
  summarise(n = sum(value))

df3 <- df %>%
  left_join(df2)

p <- df3 %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(node,": n = ", n),
             value = value)) +
  geom_sankey(alpha = 0.5) +
  geom_sankey_label(alpha = 0.5, colour = "black") +
  theme_sankey() +
  theme(legend.position = "none")

p

# only changes
df <- changes_df %>%
  filter(stable == "changed") %>%
  ggsankey::make_long(
    lg2013_label, lg2016_label, lg2019_label,
    value = count
  )

df2 <-  df %>%
  group_by(x, node) %>%
  summarise(n = sum(value))

df3 <- df %>%
  left_join(df2)

p <- df3 %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(node,": n = ", n),
             value = value)) +
  geom_sankey(alpha = 0.5) +
  geom_sankey_label(alpha = 0.5, colour = "black") +
  theme_sankey() +
  theme(legend.position = "none")

p




