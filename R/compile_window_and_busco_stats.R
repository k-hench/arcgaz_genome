# input:
# - expand( "results/pinniped/gerp/tsv/gerp_snps_{scf}_summary.tsv.gz", scf = MSCFS )
# - expand( "results/pinniped/gerp/tsv/gerp_busco_{scf}_summary.tsv.gz", scf = MSCFS )
# - expand( "results/pinniped/fst/tsv/fst_snps_{scf}_summary.tsv.gz", scf = MSCFS )
# - expand( "results/pinniped/fst/tsv/fst_busco_{scf}_summary.tsv.gz", scf = MSCFS )
# - "data/genomes/arcgaz_anc_h1.genome"
# - "results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"
# - "results/fst_win.tsv.gz"
# output:
# - "results/pinniped/go_terms/thresholds.tsv"
# - "results/pinniped/busco_gerp_fst.tsv"
# - "results/pinniped/win_gerp_fst.tsv"
# - "results/pinniped/win_outlier_summary.tsv"
library(tidyverse)
library(glue)
library(here)

scfs <- str_c("mscaf_a1_",
              c(str_pad(1:17, width = 2, pad =  0), "x"))

read_win_gerp <- \(scf){read_tsv(here(glue("results/pinniped/gerp/tsv/gerp_snps_{scf}_summary.tsv.gz")))}
read_busco_gerp <- \(scf){read_tsv(here(glue("results/pinniped/gerp/tsv/gerp_busco_{scf}_summary.tsv.gz")))}
read_win_fst <- \(scf){read_tsv(here(glue("results/pinniped/fst/tsv/fst_snps_{scf}_summary.tsv.gz")))}
read_busco_fst <- \(scf){read_tsv(here(glue("results/pinniped/fst/tsv/fst_busco_{scf}_summary.tsv.gz")))}

genome <- read_tsv(here("data/genomes/arcgaz_anc_h1.genome"),
                   col_names = c("chr", "length")) |>
  mutate(end = cumsum(length),
         start = lag(end, default = 0),
         mid = (start + end)/2,
         eo = row_number() %% 2)

plain_busco <- read_tsv(here("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"), skip = 2) |>
  filter(Status == "Complete") |>
  mutate(bstart = if_else(`Gene Start` < `Gene End`, `Gene Start`, `Gene End`),
         bend = if_else(`Gene Start` < `Gene End`, `Gene End`, `Gene Start`)) |>
  select(chr = Sequence,
         bstart, bend,
         busco_id = `# Busco id`)  |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = bstart + scaf_start,
         gend = bend + scaf_start,
         gmid = (gstart + gend)/2)

data_win_gerp <- scfs |>
  map_dfr(read_win_gerp) |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = start + scaf_start,
         gend = end + scaf_start,
         gmid = (gstart + gend)/2)

data_win_fst_vt <- read_tsv(here("results/fst_win.tsv.gz")) |>
  left_join(genome |> select(CHROM = chr, start = start, eo)) |>
  mutate(gstart = BIN_START + start,
         gend = BIN_END + start,
         gmid = (gstart + gend)/2,
         fst_w_0 = if_else(WEIGHTED_FST > 0, WEIGHTED_FST, 0),
         fst_m_0 = if_else(MEAN_FST > 0, MEAN_FST, 0))

data_win_fst <- scfs |>
  map_dfr(read_win_fst) |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = start + scaf_start,
         gend = end + scaf_start,
         gmid = (gstart + gend)/2)


read_coverages <- \(scf, type){
  read_tsv(here(glue("results/neutral_tree/cov/by_{type}/combined/combined-{scf}.tsv.gz")))
}

coverage_win <- scfs |> map_dfr(read_coverages, type = "win") |>
  mutate(all_min_2_75 = cov_ota > 2 & cov_pho > 2 & cov_min_2_ota > 0.75 & cov_min_2_pho > 0.75 )

data_win <- data_win_gerp |>
  select(chr, start, end, eo, gstart,gend, gmid, gerp_rs_mean, gerp_rs_med) |>
  left_join(data_win_fst_vt |>
              select(chr = CHROM, end = BIN_END, n_snps = N_VARIANTS, fst_w_0, fst_m_0) |>
              left_join(data_win_fst |>
                          select(chr, end, fst_mean, fst_med),
                        by = c("chr", "end"))) |>
  left_join(coverage_win |> select(chr, start, all_min_2),
            by = c("chr", "start")) |>
  filter(n_snps >= 500,
         all_min_2)

data_win |> write_tsv(here("results/pinniped/win_gerp_fst.tsv"))

data_buso_gerp <- scfs |>
  map_dfr(read_busco_gerp) |>
  left_join(plain_busco)

data_buso_fst <- scfs |>
  map_dfr(read_busco_fst)

coverage_busco <- scfs |> map_dfr(read_coverages, type = "busco")

data_busco <- data_buso_gerp |>
  select(busco_id, chr,
         bstart,   bend,
         gstart, gend, gmid,
         neutral_n_mean, neutral_n_med,
         gerp_rs_mean, gerp_rs_med) |>
  left_join(data_buso_fst |>
              select(busco_id,
                     n_snps_fst = n_snps,
                     fst_mean,
                     fst_med)) |>
  left_join(coverage_busco |> select(busco_id, all_min_2)) |>
  filter(n_snps_fst >= (bend-bstart)*0.01,
         all_min_2)

# focal_gerp <- "gerp_rs_mean"
# gerp_lab <- "average RS score (GERP)"
# gerp_tag <- "mean"
# focal_fst <- "fst_m_0"
# fst_lab <- "mean"
# focal_fst_b <- "fst_mean"

# within BUSCO comparison
querry_thresholds <- \(stat, type, by = "BUSCO"){
  dat <- list(BUSCO  = data_busco,
              win = data_win)
  prbs <- list(gerp = 1 - c(.99, .95, .975, .995, .5, .025, .01,.0001),
               fst = c(.5, .95, .975, .99, .995, .999, .9995, .9999, .99995))
  trs <- quantile(dat[[by]][[stat]],
                  probs = prbs[[type]],
                  na.rm = TRUE)
  tibble(stat = stat,
         type = type,
         threshold = trs,
         prob = names(trs),
         by = by)
}

all_thresholds <- tibble(stat = rep(c("gerp_rs_mean", "fst_mean"), 2),
                         type = rep(c("gerp", "fst"),2),
                         by = rep(c("BUSCO", "win"), each = 2)) |>
  pmap_dfr(querry_thresholds)

all_thresholds |>
  write_tsv(here("results/pinniped/go_terms/thresholds.tsv"))

get_threshold <- \(tp = "gerp", prb = "2.5%", b = "BUSCO"){
  all_thresholds |> filter(type == tp, by == b, prob == prb) |>  pluck("threshold")
}

data_busco |>
  mutate(fst_top_10 = fst_mean > get_threshold("fst", "99%"),
         gerp_low_10 = gerp_rs_mean < get_threshold("gerp", "1%"),
         gerp_top_10 = gerp_rs_mean > get_threshold("gerp", "99%"),
         fst_top_25 = fst_mean > get_threshold("fst", "97.5%"),
         gerp_low_25 = gerp_rs_mean < get_threshold("gerp", "2.5%"),
         gerp_top_25 = gerp_rs_mean > get_threshold("gerp", "97.5%")) |>
  write_tsv(here("results/pinniped/busco_gerp_fst.tsv"))

data_win_tr <- data_win |>
  mutate(gerp_top = gerp_rs_mean >= get_threshold("gerp", "99.99%", "win"),
         fst_top = if_else(is.na(fst_mean),
                           FALSE,
                           fst_mean >= get_threshold("fst", "99.99%", "win"))) |>
  group_by(chr) |>
  mutate(gerp_ol_count = cumsum(gerp_top & !lag(gerp_top, default = FALSE)),
         gerp_ol_id = if_else(gerp_top, str_c(chr,"_gerp_", gerp_ol_count), NA),
         fst_ol_count = cumsum(fst_top & !lag(fst_top, default = FALSE)),
         fst_ol_id = if_else(fst_top, str_c(chr,"_fst_", fst_ol_count), NA)) |>
  ungroup() |>
  select(-c(gerp_top, gerp_ol_count, fst_ol_count))

data_win_tr |>
  write_tsv(here("results/pinniped/win_gerp_fst.tsv"))

outlier_summary_gerp <- data_win_tr |>
  filter(complete.cases(gerp_ol_id)) |>
  filter(all_min_2) |>
  group_by(gerp_ol_id) |>
  summarise(chr = chr[[1]],
            start = min(start),
            end = max(end)) |>
  ungroup() |>
  mutate(length = end - start + 1,
          outlier_type = "GERP",
         outlier_label = str_c("g", row_number())) |>
  rename(outlier_id = "gerp_ol_id")

outlier_summary_fst <- data_win_tr |>
  filter(complete.cases(fst_ol_id)) |>
  filter(all_min_2) |>
  group_by(fst_ol_id) |>
  summarise(chr = chr[[1]],
            start = min(start),
            end = max(end)) |>
  ungroup() |>
  mutate(length = end - start + 1,
         outlier_type = "fst",
         outlier_label = str_c("f", row_number())) |>
  rename(outlier_id = "fst_ol_id")

bind_rows(outlier_summary_gerp,
          outlier_summary_fst)  |>
  left_join(genome |> select(chr, c_start = start)) |>
  mutate(gmid = (start +  end) / 2 + c_start) |>
  select(- c_start) |>
  write_tsv(here("results/pinniped/win_outlier_summary.tsv"))
