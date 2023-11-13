library(tidyverse)
library(here)

genome <- read_tsv(here("data/genomes/arcgaz_anc_h1.genome"),
                   col_names = c("seq", "length")) |>
  mutate(end = cumsum(length),
         start = lag(end, default = 0),
         mid = (start + end)/2,
         eo = row_number() %% 2)

read_cov <- \(fl){
  read_tsv(str_c("~/Downloads/vcf_test/cov/", fl),
           col_names = c("seq", "cstart", "cend", "cov"),
           col_types = "ciii")
}

files_cov <- dir("~/Downloads/vcf_test/cov/")

data_cov <- files_cov |>
  map_dfr(read_cov) |>
  left_join(genome |> select(seq, start)) |>
  mutate(gstart = start + cstart,
         gend = start + cend,
         gmid = (gstart + gend)/2)

table(data_cov$cov)

cum_cov2 <- data_cov |>
  mutate(clength = (cend - cstart)) |>
  group_by(cov) |>
  summarise(len = as.double(sum(clength))) |>
  ungroup() |>
  mutate(cum_len = lag(cumsum(len),default = 0)) |>
  bind_rows(tibble(cov = 10, len = as.double(NA), cum_len = 2382865058))

cum_cov2 |>
  ggplot(aes(x = cum_len, y = cov))+
  geom_step() +
  scale_y_continuous(breaks = 1:10)

seqs <- genome$seq
w_sz <- 50000

win_cov <- data_cov |>
  group_by(seq) |>
  mutate(cmid = (cstart+cend)/2,
         win = (cmid %/% w_sz)) |>
  group_by(seq, win) |>
  summarise(avg_cov = mean(cov)) |>
  mutate(wmid = win * w_sz + .5 * w_sz) |>
  ungroup() |>
  left_join(genome |> select(seq, start)) |>
  mutate(gmid = start + wmid)

ggplot() +
  geom_rect(data = genome,
            aes(ymin = -Inf, ymax = Inf, xmin = start, xmax = end, fill = factor(eo))) +
  geom_point(data = win_cov,
             aes(x = gmid,
                 y = avg_cov),
             size = .2, color = rgb(0,0,0,.4)) +
  scale_fill_manual(values = c(`0` = "transparent",
                             `1` = rgb(.4,.4,.4,.4)),
                  guide = "none") +
  scale_x_continuous(labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(breaks = genome$mid,
                                         labels = str_remove(genome$seq, "mscaf_a1_"),
                                         trans = identity),
                     expand = c(0, 0)) +
  # coord_cartesian(ylim = c(-.005, .005)) +
  labs(y = "FST",
       x = "Genomic position (bp)") +
  theme_minimal() #+
#theme(axis.text.y = element_blank())
