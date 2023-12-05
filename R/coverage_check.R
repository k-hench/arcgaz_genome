library(tidyverse)
library(here)
library(glue)
library(patchwork)
source(here("R/plot_defaults.R"))

scfs <- str_c("mscaf_a1_", c(str_pad(1:17, width = 2, pad = 0), "x"))

read_coverages <- \(scf, type){
  read_tsv(here(glue("results/neutral_tree/cov/by_{type}/combined/combined-{scf}.tsv.gz")))
}

summarize_cov <- \(dat){
  dat |>
    mutate(ota_min_2 = cov_ota > 2,
           pho_min_2 = cov_pho > 2) |>
    select(ends_with("_min_2")) |>
    pivot_longer(everything()) |>
    group_by(name, value) |>
    count() |>
    ungroup()
}

coverage_busco <- scfs |> map_dfr(read_coverages, type = "busco")
coverage_win <- scfs |> map_dfr(read_coverages, type = "win")

coverage_busco_summary <- coverage_busco |> summarize_cov()
coverage_win_summary <- coverage_win |> summarize_cov()

coverage_busco_summary |> pivot_wider(values_from = n, names_from = value)
coverage_win_summary |> pivot_wider(values_from = n, names_from = value)

p1 <- coverage_busco_summary |> ggplot() + labs(subtitle = "BUSCOS")
p2 <- coverage_win_summary |> ggplot() + labs(subtitle = "win")

p1 / p2 +
  plot_layout(guide = "collect") &
  geom_bar(stat = "identity",
           aes(x = value,
               y = n,
               fill = name),
           position = "dodge") &
  theme_ms(fontsize = 14)
