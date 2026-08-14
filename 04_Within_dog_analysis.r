library(here)
library(tidyverse)

# Read Bray-Curtis distance matrix generated in Python
dist_matrix <- read.csv(
  here("braycurtis_distance_matrix.csv"),
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

storage.mode(dist_matrix) <- "numeric"

# Continuous-infection comparisons used to define the reference range
continuous_pairs <- tribble(
  ~dog,     ~sample_time1,       ~sample_time2,
  "74FU1",  "74FU1_08_23",      "74FU1_08_24",
  "74JS11", "74JS11_08_23",     "74JS11_05_24",
  "74DA4",  "74DA4_10_23",      "74DA4_02_24",
  "74AN23", "74AN23_08_24",     "74AN23_08_25",
  "74GE9",  "74GE9_07_24",      "74GE9_06_25",
  "WSTI11", "WSTI11_08_24",     "WSTI11_08_25",
  "WSMI8",  "WSMI8_05_24",      "WSMI8_08_25"
)

# Reinfection comparisons
reinfection_pairs <- tribble(
  ~dog,     ~sample_time1,       ~sample_time2,
  "74FU1",  "74FU1_08_23",      "74FU1_09_25",
  "74FU1",  "74FU1_08_24",      "74FU1_09_25",

  "74DA4",  "74DA4_10_23",      "74DA4_11_25",
  "74DA4",  "74DA4_02_24",      "74DA4_11_25",

  "74JS11", "74JS11_08_23",     "74JS11_11_25",
  "74JS11", "74JS11_05_24",     "74JS11_11_25",

  "74TE56", "74TE56_07_24",     "74TE56_11_25",
  "74TU17", "74TU17_04_24",     "74TU17_08_25",
  "74WB19", "74WB19_04_24",     "74WB19_11_25",
  "74RI13", "74RI13_05_24",     "74RI13_11_25",
  "74JO10", "74JO10_05_24",     "74JO10_09_25"
)

# Check that all required samples are present
all_pairs <- bind_rows(
  continuous_pairs,
  reinfection_pairs
)

required_samples <- unique(
  c(
    all_pairs$sample_time1,
    all_pairs$sample_time2
  )
)

missing_samples <- setdiff(
  required_samples,
  rownames(dist_matrix)
)

if (length(missing_samples) > 0) {
  stop(
    paste0(
      "Samples not found in distance matrix:\n",
      paste(missing_samples, collapse = "\n")
    )
  )
}

# Extract Bray-Curtis distances for continuous-infection comparisons
continuous_distances <- continuous_pairs %>%
  rowwise() %>%
  mutate(
    distance = dist_matrix[sample_time1, sample_time2]
  ) %>%
  ungroup() %>%
  mutate(
    group = "Continuous infection"
  )

print(continuous_distances)

# Mean and SD reference ranges from continuous infections
continuous_summary <- continuous_distances %>%
  summarise(
    n = n(),
    mean_distance = mean(distance),
    sd_distance = sd(distance),
    lower_2sd = max(0, mean_distance - 2 * sd_distance),
    upper_2sd = min(1, mean_distance + 2 * sd_distance),
    lower_3sd = max(0, mean_distance - 3 * sd_distance),
    upper_3sd = min(1, mean_distance + 3 * sd_distance)
  )

print(continuous_summary)

mean_distance <- continuous_summary$mean_distance
lower_2sd <- continuous_summary$lower_2sd
upper_2sd <- continuous_summary$upper_2sd
lower_3sd <- continuous_summary$lower_3sd
upper_3sd <- continuous_summary$upper_3sd

# Extract Bray-Curtis distances for reinfection comparisons
reinfection_distances <- reinfection_pairs %>%
  rowwise() %>%
  mutate(
    distance = dist_matrix[sample_time1, sample_time2]
  ) %>%
  ungroup() %>%
  mutate(
    group = "Reinfection"
  )

print(reinfection_distances)

# Combine comparisons and create plotting labels
plot_distances <- bind_rows(
  reinfection_distances,
  continuous_distances
) %>%
  mutate(
    month1 = str_extract(sample_time1, "(?<=_)\\d{2}(?=_\\d{2}$)"),
    year1 = str_extract(sample_time1, "\\d{2}$"),
    month2 = str_extract(sample_time2, "(?<=_)\\d{2}(?=_\\d{2}$)"),
    year2 = str_extract(sample_time2, "\\d{2}$"),

    date1_order = as.numeric(paste0("20", year1, month1)),
    date2_order = as.numeric(paste0("20", year2, month2)),

    time1_label = paste0(month1, "/", year1),
    time2_label = paste0(month2, "/", year2),

    comparison = paste0(
      dog,
      "\n",
      time1_label,
      " vs ",
      time2_label
    ),

    group = factor(
      group,
      levels = c(
        "Reinfection",
        "Continuous infection"
      )
    )
  ) %>%
  arrange(
    group,
    dog,
    date1_order,
    date2_order
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = unique(comparison)
    )
  )

print(plot_distances)

# Save extracted distances and reference statistics
write.csv(
  continuous_distances,
  here("continuous_infection_distances.csv"),
  row.names = FALSE
)

write.csv(
  continuous_summary,
  here("continuous_infection_mean_SD.csv"),
  row.names = FALSE
)

write.csv(
  reinfection_distances,
  here("reinfection_all_year_comparisons.csv"),
  row.names = FALSE
)

write.csv(
  plot_distances,
  here("all_within_dog_comparisons_for_plot.csv"),
  row.names = FALSE
)

# Set y-axis range
largest_value <- max(
  plot_distances$distance,
  upper_3sd,
  na.rm = TRUE
)

y_upper <- min(
  1,
  max(
    0.7,
    ceiling(largest_value * 10) / 10
  )
)

# Plot within-dog Bray-Curtis distances
p <- ggplot() +

  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = lower_3sd,
    ymax = upper_3sd,
    fill = "grey90"
  ) +

  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = lower_2sd,
    ymax = upper_2sd,
    fill = "grey75"
  ) +

  geom_hline(
    yintercept = mean_distance,
    linewidth = 1,
    color = "black"
  ) +

  geom_point(
    data = plot_distances,
    aes(
      x = comparison,
      y = distance,
      color = group
    ),
    shape = 16,
    size = 5
  ) +

  scale_color_manual(
    values = c(
      "Continuous infection" = "black",
      "Reinfection" = "blue"
    ),
    name = NULL
  ) +

  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(
      mult = c(0.01, 0.05)
    )
  ) +

  coord_cartesian(
    ylim = c(0, y_upper)
  ) +

  labs(
    x = NULL,
    y = "Bray-Curtis distance"
  ) +

  theme_classic() +

  theme(
    legend.position = "right",

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 15,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 16,
      color = "black"
    ),

    axis.title.y = element_text(
      size = 16,
      color = "black"
    ),

    legend.text = element_text(
      size = 16
    )
  )

print(p)

# Save figure
ggsave(
  here("within_dog_distances_with_continuous_SD_ranges.pdf"),
  plot = p,
  width = 14,
  height = 6
)

ggsave(
  here("within_dog_distances_with_continuous_SD_ranges.png"),
  plot = p,
  width = 14,
  height = 6,
  dpi = 300
)