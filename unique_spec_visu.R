library(tidyverse)

species_rich <- obs_m %>% 
  mutate(
    month_date = as.Date(month_date),                 # POSIXct → Date
    richness   = map_int(species_list, length)        # n unique spp / month
  ) %>% 
  ggplot(aes(month_date, richness)) +
  geom_col() +                                        # bars = richness
  scale_x_date(date_breaks = "1 month",
               date_labels = "%b\n%Y") +
  labs(title = "Monthly species richness",
       x = NULL, y = "Number of species") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



unique_by_month <- obs_m %>% 
  ## 1. explode the list-column: one row = one species-month record
  unnest_longer(species_list, values_to = "species") %>% 
  
  ## 2. find species that appear in *exactly one* month
  add_count(species, name = "n_months") %>% 
  filter(n_months == 1) %>%                   # keep the rare ones
  
  ## 3. re-collect by month
  group_by(month_date) %>% 
  summarise(
    unique_species = list(species),           # character vector per month
    n_unique       = n(),                     # how many such species
    .groups = "drop"
  ) %>% 
  arrange(desc(n_unique))                     # months with most uniques first

print(unique_by_month, n = Inf)   # shows month, count, and the species list

unique_by_month %>% 
  ggplot(aes(as.Date(month_date), n_unique)) +
  geom_col() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(title = "Number of species seen in only that month",
       x = NULL, y = "Unique-to-month species") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
