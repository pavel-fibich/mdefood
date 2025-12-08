rm(list=ls())
# Load libraries----
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggh4x)
library(glue)
# Data setup ----
fwlist <- c("Plowmanetal", "Krkonose", "CameroonW", "CameroonD")
fwlist2 <- c("Myrmecophytic\nant-plant dataset",
             "Temperate\nplant-pollinator dataset",
             "Tropical wet-season\nplant-pollinator dataset",
             "Tropical dry-season\nplant-pollinator dataset")

fwlen <- c(4, 4, 4, 4)

elevations_plowman <- c(700, 935, 1170, 1400)
elevations_krkonose <- c(450, 600, 800, 1000)
elevations_cameroon <- c(650, 1100, 1450, 2200)
names(elevations_cameroon) = names(elevations_krkonose) = names(elevations_plowman) = c(1,2,3,4)

elevations_list <- list(
  elevations_plowman,
  elevations_krkonose,
  elevations_cameroon,
  elevations_cameroon
)
names(elevations_list) <- fwlist2

# Read in datasets----
datasets <- 
  fwlist %>% 
  map(~ read_csv(paste0("realdatasets/netsgrad1000noi_", ., "x1.csv")))
# Set names
names(datasets) <- fwlist2

# bind rows, pivot longer
datasets <-
  datasets %>%
  bind_rows(.id = "dataset_name") %>% 
  mutate(
    gr = map2_dbl(dataset_name, gr, ~ elevations_list[[.x]][as.character(.y)])
  )  %>%
  mutate(dataset_name = 
           factor(dataset_name, 
                  levels = fwlist2)) %>% 
  pivot_longer(NOconnections:vulnerability.LL,
               names_to = "index",
               values_to = "index_value") %>% 
  mutate(index = case_when(index == "NOconnections" ~ "Total~realised~links",
                           #index == "SA" ~ "Specialisation~asymmetry",
                           index == "generality.HL" ~ "Generality",
                           index == "vulnerability.LL" ~ "Vulnerability",
                           index == "H2" ~ "H[2]'",
                           index == "nestedness" ~ "NestednessOld",
                           index == "NODF" ~ "Nestedness",
                           index == "connectance" ~ "Connectance",
                           .default = index)) %>% 
  mutate(index = factor(index, 
                        levels = 
                          c("Total~realised~links",
                            "Connectance",
                            "Nestedness",
                            #"NODF",
                            #"Specialisation~asymmetry",
                            "H[2]'",  # Apply bquote() for H2 in fct_relevel
                            "Generality",
                            "Vulnerability"),
                        labels = 
                          c("Total~realised~links",
                            "Connectance",
                            "Nestedness",
                            #                          "NODF",
                            #"Specialisation~asymmetry",
                            bquote(italic("H"[2]*"'")),  # Apply bquote() for H2 in fct_relevel
                            "Generality",
                            "Vulnerability")))

datasets <- datasets %>% filter(!is.na(index))

summary_nulls_1 <- 
  datasets %>% 
  filter(null != 0) %>% 
  group_by(dataset_name, gr, index) %>% 
  summarise(lower_quantile = quantile(index_value, 0.025, na.rm = T),
            upper_quantile = quantile(index_value, 0.95, na.rm = T),
            sd = sd(index_value, na.rm = T),
            mean = mean(index_value, na.rm = T))

summary_nulls_2 <- 
  summary_nulls_1 %>% 
  left_join(datasets %>% filter(null != 0)) %>% 
  group_by(dataset_name, index, gr, null) %>% 
  summarise(ses_null_gr = (abs(index_value - mean)/abs(sd))) %>% 
  distinct() %>% 
  group_by(dataset_name, index, null) %>% 
  summarise(ses_null = sum(ses_null_gr)) %>% 
  summarise(upper_quantile_ses = quantile(ses_null, 0.95, na.rm = T),
            lower_quantile_ses = quantile(ses_null, 0.025, na.rm = T))

index_points <- datasets %>% 
  filter(null == 0) %>% 
  left_join(summary_nulls_1) %>% 
  group_by(dataset_name, index, gr, null) %>% 
  mutate(index_shape = case_when(index_value > upper_quantile |
                                   index_value < lower_quantile ~ "2",
                                 .default = "1"))

stars_true <-
  datasets %>% 
  filter(null == 0) %>% 
  left_join(summary_nulls_1) %>% 
  group_by(dataset_name, index, gr, null) %>% 
  summarise(ses_true_gr = (abs(index_value - mean)/abs(sd))) %>% 
  distinct() %>% 
  group_by(dataset_name, index, null) %>% 
  summarise(ses_true_gr = sum(ses_true_gr)) %>% 
  left_join(summary_nulls_2) %>% 
  mutate(sig = 
           case_when(ses_true_gr > upper_quantile_ses |
                       ses_true_gr < lower_quantile_ses ~ "*",
                     .default = "")) %>% 
  distinct(dataset_name, index, sig)




plot1 <- 
  datasets %>% 
  filter(null != 0) %>% 
  group_by(dataset_name, gr, index) %>% 
  summarise(lower_quantile = quantile(index_value, 0.025, na.rm = T),
            upper_quantile = quantile(index_value, 0.95, na.rm = T)) %>% 
  
  ggplot(aes(x = gr, fill = index, colour = index)) +
  
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_quantile, ymax = upper_quantile), alpha = 0.25, colour = NA) +
  geom_line(aes(y = (upper_quantile-(upper_quantile - lower_quantile)/2)),
            # alpha = 0.75,
            linewidth = 1,
            linetype = "dotted") +
  
  # True points
  geom_point(data = index_points,
             aes(x = gr, y = index_value, shape = index_shape), colour = "black",
             size = 3) +
  scale_shape_manual(values = c("1" = 1, "2" = 24)) +
  geom_line(data = index_points,
            aes(x = gr, y = index_value)) +
  
  # Significance stars
  geom_text(data = stars_true, aes(x = Inf, y = Inf, label = sig),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, size = 12,
            colour = "red") +
  
  # Faceting
  facet_grid2(index ~ dataset_name,
              labeller = labeller(
                index = label_parsed,  # Parse 'index' variable
                dataset_name = label_value),  # Leave 'dataset_name' as regular text
              scales = "free",
              independent = "y",
              switch = "y",
              axes = "all") +
  coord_cartesian(clip = 'off') +
  # Theming
  theme_classic(base_size = 15) +
  theme(
    plot.margin = unit(c(0, 1, 1, 0), "lines"),  # Widens the right margin
    legend.position = "none",
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0, "cm"),
    strip.background = element_blank(), 
    strip.text.y.left = element_text(angle = 0),
    panel.spacing.x = unit(0.5, "cm"),
    # panel.spacing.y = unit(0.5, "cm")
  ) +
  labs(x = "Elevation (m a.s.l.)", y = "")
plot1


ggsave(filename = paste0("realfw2_EL", ".pdf"), plot = plot1, device = NULL, width = 12, height = 10.69, dpi = 300)
