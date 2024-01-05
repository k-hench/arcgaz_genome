library(tidyverse)
library(glue)
library(prismatic)
library(ggtext)
library(ggrastr)
library(here)
library(arcgazgen)
source("R/plot_defaults.R")

fs <- fnt_sz
clr_l_gray <- rgb(.85,.85,.85)

data_win  <- read_tsv(here("results/pinniped/win_gerp_fst.tsv"))

p <- ggplot() +
  geom_ag_genome()+
  rasterize(geom_point(data = data_win,
                       aes(x = gmid, y = cov_all + 1,
                           color = factor(eo),
                           group = chr),
                       size = .01,
                       shape = 19),
            dpi = 300, dev = "ragg") +
  geom_smooth(data = data_win,
              color = "white",
              linewidth = .2,
              span = 0.5,
              aes(x = gmid,
                  y = cov_all + 1,
                  group = chr),
              se = FALSE) +
  scale_color_manual(values = c("black", clr_lighten("black",.33)),
                     guide = "none") +
  scale_x_continuous(name = "Genomic Position (Gb)", labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis("Scaffold Id",
                                         trans = identity,
                                         breaks = arcgaz_genome$mid,
                                         labels = str_remove(arcgaz_genome$chr, "mscaf_a1_")))+
  scale_alpha_manual(values = c(`TRUE` = 6,
                                `FALSE` = 1),
                     guide = "none") +
  scale_y_continuous(breaks = (0:5)*2, labels = \(x){sprintf("%.0f", x)}) +
  labs(y = "Average Alignment Coverage (n genomes)") +
  coord_cartesian(xlim = c(0, max(arcgaz_genome$end)),
                  ylim = c(0, 11),
                  expand = 0)+
  theme_ms(fontsize = fs)+
  theme(axis.title.y = element_markdown())

ggsave(plot = p,
       filename = here("results/img/win_coverage.pdf"),
       width = 5.5, height = 2, device = cairo_pdf)
