library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

tsv_in <- args[[1]]
tsv_out <- args[[2]]

read_tsv(tsv_in,
         col_names = c("chr", "start","end", "win_idx","neutral_n", "gerp_rs")) |>
  filter(gerp_rs != 0) |>
  group_by(chr, start, end ,win_idx) |>
  summarise(n_snps = n(),
            across(c(neutral_n, gerp_rs),
                   list(mean = mean, med = median, sd = sd),
                   .names = "{.col}_{.fn}")) |>
  ungroup() |>
  mutate(mid = (start + 1 + end )/2) |>
  write_tsv(file = tsv_out)
