library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(ggtext)
library(here)
library(arcgazgen)
source("R/plot_defaults.R")

read_fai <- \(file){
  read_tsv(file,
           col_names = c("Scaffold", "length", "offset", "linebases", "linewidth")) |>
    select(Scaffold, length)
}

g1 <- read_fai(here("data/genomes/arcgaz_anc_h1.fa.gz.fai"))
g2 <- read_fai(here("data/genomes/arcgaz_anc_h2.fa.gz.fai"))

data_h1 <- read_tsv("~/Downloads/david_mail/SNPs_Positions_Hap1.txt") |>
  filter(grepl("mscaf", Scaffold)) |>
  group_by(Scaffold) |>
  count() |>
  left_join(g1)

data_h2 <- read_tsv("~/Downloads/david_mail/SNPs_Positions_Hap2.txt") |>
  filter(grepl("mscaf", Scaffold)) |>
  group_by(Scaffold) |>
  count() |>
  left_join(g2)

clr_fg <- "gray45"
clr_bg <- "white"

p1 <- data_h1|>
  ggplot(aes(x = Scaffold, y = n,
             fill = Scaffold != "mscaf_a1_x"))+
  geom_bar(stat = 'identity',
           color = clr_darken(clr_fg),
           linewidth = plt_lwd,
           width = .8) +
  scale_x_discrete(labels = \(x){str_remove(x,"mscaf_a1_")}) +
  coord_cartesian(xlim = c(0.4,18.6),
                  ylim = c(0, 8200),
                  expand = 0) +
  labs(y = "SNPs per scaffold (n)") +
  theme_ms()

p2 <- data_h2 |>
  ggplot(aes(x = Scaffold, y = n,
             fill = Scaffold != "mscaf_a2_x"))+
  geom_bar(stat = 'identity',
           color = clr_darken(clr_fg),
           linewidth = plt_lwd,
           width = .8) +
  scale_x_discrete(labels = \(x){str_remove(x,"mscaf_a2_")}) +
  coord_cartesian(xlim = c(0.4,18.6),
                  ylim = c(0, 8200),
                  expand = 0) +
  theme_ms() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

xlblr <- \(x){sprintf("%.1f", x*1e-6)}

p3 <- data_h1 |>
  ggplot(aes(x = length, y = n,
             fill = Scaffold != "mscaf_a1_x"))+
  geom_smooth(method = "lm",
              color = clr_fg,
              # linetype = 3,
              linewidth = plt_lwd,
              fullrange = TRUE,
              se = F) +
  geom_point(shape = 21,
             size = .9,
             color = clr_fg)   +
  labs(y = "SNPs per scaffold (n)") +
  scale_x_continuous("Scaffold length (Mb)",
                     limits = c(.1e7, 2.8e8),
                     labels = xlblr) +
  coord_cartesian(xlim = c(4.5e7,2.25e8),
                  ylim = c(1600, 8200),
                  expand = 0) +
  theme_ms()

p4 <- data_h2 |>
  ggplot(aes(x = length, y = n,
             fill = Scaffold != "mscaf_a2_x"))+
  geom_smooth(method = "lm",
              color = clr_fg,
              # linetype = 3,
              linewidth = plt_lwd,
              fullrange = TRUE,
              se = F) +
  geom_point(shape = 21,
             size = .9,
             color = clr_fg)  +
  scale_x_continuous("Scaffold length (Mb)",
                     limits = c(.1e7, 2.8e8),
                     labels = xlblr) +
  coord_cartesian(xlim = c(4.5e7,2.25e8),
                  ylim = c(1600, 8200),
                  expand = 0) +
  theme_ms() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

pp1 <- p1 + p2 +
  p3 + p4 +
  plot_layout() +
  plot_annotation(tag_levels = "a") &
  scale_fill_manual(values = c(`TRUE` = clr_fg, `FALSE` = clr_bg),
                    guide = "none") &
  theme(plot.tag = element_text(family = fnt_sel))

ggsave(filename = here("results/img/snp_stats.pdf"),
       width = 8, height = 4, device = cairo_pdf)
