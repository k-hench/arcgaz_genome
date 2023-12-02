library(tidyverse)
library(here)
library(glue)

mscaf <- args[[1]]

read_and_summarize <- \(file, type = "all"){
  miv_cov <- c(all = 4, pho = 2, ota = 2)

  read_table(here(file)) |>
    mutate(c_length = c_end - c_start ) |>
    group_by(chr, start, end, win_idx) |>
    summarise(average_cov = sum(cov * c_length) / sum(c_length),
              prcent_smaller_2 = sum(c_length[cov > (miv_cov[type] - 1)]) / sum(c_length)) |>
    ungroup() |>
    set_names(nm = c("chr", "start", "end", "win_idx",
                     str_c(c("cov_", "cov_min_"),
                           c("", miv_cov[type]),
                           c("", "_"),
                           type)))
}

data_all <- read_and_summarize(here(glue("results/neutral_tree/cov/by_win/{mscaf}.tsv.gz")), type = "all")
data_ota <- read_and_summarize(here(glue("results/neutral_tree/cov/by_win/fam/ota-{mscaf}.tsv.gz")), type = "ota")
data_pho <- read_and_summarize(here(glue("results/neutral_tree/cov/by_win/fam/pho-{mscaf}.tsv.gz")), type = "pho")

data_all |>
  left_join(data_ota |> select(-c(chr, start, end))) |>
  left_join(data_pho |> select(-c(chr, start, end))) |>
  mutate(all_min_2 = cov_ota > 2 & cov_pho > 2) |>
  write_tsv(here(glue("results/neutral_tree/cov/by_win/combined/combined-{mscaf}.tsv.gz")))
