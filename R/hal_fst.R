library(tidyverse)
library(here)

genome <- read_tsv(here("data/genomes/arcgaz_anc_h1.genome"),
                   col_names = c("seq", "length")) |>
  mutate(end = cumsum(length),
         start = lag(end, default = 0),
         mid = (start + end)/2,
         eo = row_number() %% 2)

data_fst <- read_tsv("~/Downloads/vcf_test/ota_pho_win.tsv.gz") |>
  left_join(genome |> select(CHROM = seq, start)) |>
  mutate(gstart = start + BIN_START,
         gend = start + BIN_END,
         gmid = (gstart + gend)/2)

ggplot() +
  # geom_rect(data = genome,
  #           aes(ymin = -Inf, ymax = Inf, xmin = start, xmax = end, fill = factor(eo))) +
   geom_point(data = data_fst,
             aes(x = (BIN_START+BIN_END)/2,#gmid,
                 y = WEIGHTED_FST),
             size = .2, color = rgb(0,0,0,.4)) +
  ggforce::facet_row(CHROM ~., scales = "free_x")
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
