library(tidyverse)
library(here)

genome <- read_tsv(here("data/genomes/arcgaz_anc_h1.genome"),
                   col_names = c("seq", "length")) |>
  mutate(end = cumsum(length),
         start = lag(end, default = 0),
         mid = (start + end)/2,
         eo = row_number() %% 2)

buscos <- read_tsv("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv",
                   skip = 2, col_types = "ccciicdicc") |>
  filter(Status == "Complete") |>
  left_join(genome |> select(Sequence = seq, start)) |>
  mutate(gstart = start + `Gene Start`,
         gend = start + `Gene End`,
         gmid = (gstart + gend)/2) |>
  arrange(gstart)

ggplot() +
  geom_rect(data = genome,
            aes(ymin = -Inf, ymax = Inf, xmin = start, xmax = end, fill = factor(eo))) +
  geom_density(data = buscos,
               trim = TRUE,
               adjust = .2,
             aes(x = gmid, #y = 0,
                 y = after_stat(count),
                 group = Sequence
                 ),
             fill = rgb(.3,.3,.3,.7),
             linewidth = .4)+
  scale_fill_manual(values = c(`0` = "transparent",
                               `1` = rgb(.4,.4,.4,.4)),
                    guide = "none") +
  scale_x_continuous(labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(breaks = genome$mid,
                                         labels = str_remove(genome$seq, "mscaf_a1_"),
                                          trans = identity),
                     expand = c(0, 0)) +
  # coord_cartesian(ylim = c(-.005, .005)) +
  labs(y = "BUSCO density",
       x = "Genomic position (bp)") +
  theme_minimal() +
  theme(axis.text.y = element_blank())
