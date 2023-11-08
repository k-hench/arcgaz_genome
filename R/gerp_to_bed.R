library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
tsv_in <- args[[1]]
bed_out <- args[[2]]
path_out <- args[[3]]

read_gerp <- \(file){
  read_tsv(file, col_names = c("neutral_rate_n", "rs_score")) |>
    mutate(seq = str_c("mscaf_a1_", str_replace(file, ".*pinniped_set_([0-9x]*).maf.*", "\\1")),
           pos = row_number(),
           start = pos -1)
}

dir.create(path = path_out,
           showWarnings = FALSE)

read_gerp(tsv_in) |>
  select(seq, start, pos, neutral_rate_n, rs_score) |>
  write_tsv(file = bed_out, col_names = FALSE)
