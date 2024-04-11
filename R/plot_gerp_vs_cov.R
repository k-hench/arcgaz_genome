library(tidyverse)
library(here)
library(glue)
library(prismatic)
library(ggtext)
library(ggrastr)
source(here("R/plot_defaults.R"))

read_cov <- \(scf){
  read_tsv(here(glue("results/neutral_tree/cov/by_win/combined/combined-{scf}.tsv.gz")))
}


scfs <- str_c("mscaf_a1_", c(str_pad(1:17, width = 2, pad = 0), "x"))
data_cov <- scfs |> map_dfr(read_cov) |> select(chr, start, end, cov_all)
data_win  <- read_tsv(here("results/pinniped/win_gerp_fst.tsv")) |>
  select(chr, start, end, gerp_rs_mean, fst_m_0)

data <- data_win |>
  left_join(data_cov, by = c("chr", "start", "end"))

p <- data |>
  ggplot(aes(x = cov_all, y = gerp_rs_mean)) +
  rasterise(geom_point(),dpi = 300) +
    labs(x = "mean coverage (50kb windows)",
       y = "mean GERP score (50kb windows)")

ggsave(filename = here("results/img/gerp_vs_cov.pdf"),
       plot = p,
       width = 6,
       height = 5,
       device = cairo_pdf)
