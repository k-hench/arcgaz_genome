library(tidyverse)
library(prismatic)
library(ggtext)
library(circlize)
library(glue)
library(patchwork)
library(here)
source(here("R/plot_defaults.R"))
source(here("R/cartesian_alignment_helpers.R"))

genomes <- c("arcgaz_dt_h1_hardmasked",
             "arcgaz_dt_h2_hardmasked")

sizes <- tibble(genome = genomes,
                y_base = seq_along(genomes)-1) |>
  pmap_dfr(read_size)

genome_summary <- sizes |>
  group_by(genome, y_base) |>
  summarise(n = n(),
            toal_length = max(end))

g_height <- .8
align_min_length <- 3e5

n_largest <- 37
sizes |>
  filter(size_idx < n_largest) |>
  ggplot() +
  geom_linerange(data = tibble(end = align_min_length),
                 aes(y = .9, xmin = 0, xmax = end ), size = 1.5) +
  geom_rect(data = genome_summary,
            aes(xmin = 0, xmax = toal_length,
                ymin = y_base, ymax = y_base + g_height),
            color = "gray50", fill = "gray90") +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = y_base, ymax = y_base + g_height,
                color = factor(eo),
                fill = after_scale(clr_lighten(color, .5)))) +
  geom_richtext(data = genome_summary,
                aes(y = y_base + .5 * g_height,
                    label = glue("{str_remove(genome, '_hardmasked')}<br>n: {n}"),
                    x = -.1e9),hjust = 1,
                fill = NA, label.color = NA, family = fnt_sel,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_y_continuous(NULL) +
  scale_x_continuous(NULL, breaks = c(0:3)*1e9,
                     labels = str_c(0:3, " Gb")) +
  scale_color_manual(values = clrs |> clr_lighten(.2) |> clr_desaturate(), guide = "none") +
  coord_cartesian(xlim = c(-4.5e8, 2.6e9),
                  ylim = c(-.1, length(genomes)),
                  expand = 0)+
  # coord_polar() +
  labs(title = "pre-anchoring genomes",
       subtitle = "(scaffolds sorted by size)") +
  theme_minimal(base_family = fnt_sel) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = .2),
        axis.ticks.length = unit(-4,"pt"),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

clr_genome = "gray50"
genome_width <- .05
skip <- .025

manual_order_18 <- list(c(14,9,12,4, 17,13,10,2,1,18,11,15,7,5,8,3,16,6),
                        c( 18, 28, 33, 19, 24,
                           26, 8, 15, 5, 1,
                           23, 29, 22, 13, 16,
                           36, 12, 27, # ----
                           6, 4,
                           21, 9, 25, 17, 34,
                           10, 20, 31, 30, 2,
                           14, 11, 32, 35, 7,
                           3),
                        c( 29, 33, 19, 24, 26, 1, 7, 9, 18, 12, 4, 8, 23, 31, 22,  5, 16, 17, #-----
                           3, 35, 21, 13, 10, 25, 28, 30, 34, 20, 2, 14, 11, 32, 36, 27,6, 15)) |>
  set_names(nm = c( "arcgaz_v3_hardmasked", "arcgaz_dt_h1_hardmasked", "arcgaz_dt_h2_hardmasked"));

all_alignments <- tibble(genomes = list(c("arcgaz_dt_h1_hardmasked", "arcgaz_dt_h2_hardmasked")),
                         n_first = list(c(36,36))) |>
  pmap_dfr(import_alignment,
           genome_width = genome_width,
           skip = skip,
           manual_order_18 = manual_order_18,
           n_longest = 25)
# p1 <-

all_alignments |>
  unnest(cols = psl_diag) |>
  ggplot() +
  facet_grid(querry ~ .,switch = "y") +
  geom_diagonal_wide(aes(x = y,
                         y = x,
                         group = str_c(querry, "_" ,group),
                         color = factor(dir),
                         fill = after_scale(clr_lighten(color))),
                     orientation = "y", linewidth = plt_lwd)  +
  geom_linerange(data = all_alignments |> unnest(genome_with_skips),
                 aes(xmin = start, xmax = end, y = y),
                 color = clr_genome,
                 linewidth = .5 * plt_lwd) +
  geom_rect(data = all_alignments |> unnest(sizes_18),
            aes(xmin = start, xmax = end,
                ymin = y_base + (label_sign * skip),
                ymax = y_base + (label_sign * (skip + genome_width))),
            color = clr_genome, fill = clr_lighten(clr_genome,.75),
            linewidth = .5 * plt_lwd) +
  # geom_text(data = all_alignments |> unnest(sizes_18),
  #           aes(x = mid, y = y_base, label = str_c(size_idx,"\n",pre_manual))) +
  # geom_ribbon(data = all_alignments |> unnest(data_cov_windowed),
  #             aes(x = gg_mid, group = seqnames,
  #                 ymin = (genome_idx -1) + 2 * skip *(genome_idx - 1.5) + genome_width *(genome_idx - 2),
  #                 ymax = (genome_idx -1) + 2 * skip *(genome_idx - 1.5) + genome_width *(genome_idx - 2) +  genome_width * avg_cov_25),
  #             linewidth = plt_lwd,
  #             fill = clr_genome) +
  geom_text(data = all_alignments |> unnest(sizes_18),
            aes(x = mid,
                y = y_base + (label_sign *  (4.3 * skip + genome_width)),
                label = simplify_names(chr)#,
                # hjust = c(1, 0)[as.numeric(factor(label_sign))]
            ),
            color = "black",
            family = fnt_sel,
            size = .35 * fnt_sz * 2,
            angle = 90)  +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.2, 1.25),
                  xlim = c(-1e7, max((all_alignments |> unnest(sizes_18))$end_with_skip)),
                  clip = "off") +
  scale_color_manual(values = clrs |> clr_alpha(.3)) +
  theme_void(base_family = fnt_sel) +
  theme(legend.position = "none",
        strip.text = element_text(angle = 90))


ggsave("img/initial_haplotype_alignment.pdf", width = 14, height = 3.5, device = cairo_pdf)
