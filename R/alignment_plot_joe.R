library(tidyverse)
library(prismatic)
library(ggtext)
library(glue)
library(patchwork)
library(here)
source(here("R/plot_defaults.R"))
source(here("R/cartesian_alignment_helpers.R"))


genomes <- c("arcgaz_anc_h1",
             "arcgaz_anc_h2",
             "arcgaz_v3_hardmasked",
             "zalcal_v1")

sizes <- tibble(genome = genomes,
                y_base = seq_along(genomes)-1) |>
  pmap_dfr(read_size)

genome_summary <- sizes |>
  group_by(genome, y_base) |>
  summarise(n = n(),
            toal_length = max(end))

g_height <- .8
align_min_length <- 3e5
n_largest <- 18
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
  labs(title = "anchored genomes",
       subtitle = "(scaffolds sorted by size)") +
  theme_minimal(base_family = fnt_sel) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = .2),
        axis.ticks.length = unit(-4,"pt"),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

genome_summary <- sizes |>
  group_by(genome, y_base) |>
  summarise(n = n(),
            toal_length = max(end))

g_height <- .8
align_min_length <- 3e5

clr_genome = "gray50"
genome_width <- .05
skip <- .025

all_alignments <- tibble(genomes = list(c("arcgaz_anc_h2", "zalcal_v1"),
                                        c("arcgaz_anc_h1", "zalcal_v1"))) |>
  pmap_dfr(import_alignment,
           genome_width = genome_width,
           skip = skip,
           manual_order_18 = manual_order_18,
           n_longest = 25)
# p1 <-
fnt_sel <- "Gill sans"
skip_y <- .06
clrs <- scales::colour_ramp(rcartocolor::carto_pal(12,"Prism"))( ((1:18) - 1 ) /17 )


all_alignments |>
  filter(querry == "zalcal_v1") |>
  unnest(cols = psl_diag) |>
  mutate(x = if_else(target == "arcgaz_anc_h2", x - 1 - skip_y, - x + 1 +skip_y)) |>
  ggplot() +
  # facet_grid(target + querry ~ .,switch = "y") +
  geom_diagonal_wide(aes(x = y,
                         y = x,
                         group = str_c(querry, "_", target,"_",group),
                         color = factor(dir),
                         fill = after_scale(clr_lighten(color))),
                     orientation = "y", linewidth = plt_lwd)  +
  geom_rect(data = all_alignments |>
              filter(querry == "zalcal_v1") |>
              unnest(sizes_18) |>
              filter(!(target == "arcgaz_anc_h2" & genome == "zalcal_v1")) |>
              mutate(y_base = case_when(
                genome == "arcgaz_anc_h1" ~ 1 + 2 * skip_y,
                genome == "arcgaz_anc_h2" ~ -1 - 2 * skip_y,
                genome == "zalcal_v1" ~ 0,
              )),
            aes(xmin = start, xmax = end,
                ymin = y_base - .5 *skip_y ,
                ymax = y_base + .5 *skip_y ),
            color = clr_genome, fill = clr_lighten(clr_genome,.75),
            linewidth = .5 * plt_lwd) +
  geom_text(data = all_alignments |>
              filter(querry == "zalcal_v1") |>
              unnest(sizes_18) |>
              filter(!(target == "arcgaz_anc_h2" & genome == "zalcal_v1")) |>
              mutate(y_base = case_when(
                genome == "arcgaz_anc_h1" ~ 1 + 2 * skip_y,
                genome == "arcgaz_anc_h2" ~ -1 - 2 * skip_y,
                genome == "zalcal_v1" ~ 0,
              )),
            aes(x = mid,
                y = y_base,
                label = simplify_names(chr)#,
                # hjust = c(1, 0)[as.numeric(factor(label_sign))]
            ),
            color = "black",
            family = fnt_sel,
            size = .35 * fnt_sz * 1,
            angle = 0)  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = -1:1 * (1+skip_y*2),
                     labels = c("haplotype 2", "*Z. californianus*", "haplotype 1")) +
  coord_cartesian(#ylim = c(-.2, 1.25),
                  xlim = c(-1e7, max((all_alignments |> unnest(sizes_18))$end_with_skip)),
                  clip = "off") +
  scale_color_manual(values = clrs) +
  theme_void(base_family = fnt_sel) +
  theme(legend.position = "none",
        strip.text = element_text(angle = 90),
        axis.text.y = element_markdown())

ggsave("~/Dropbox/joe_hoffman/figures/haplo_alignment_zalcal.png", width = 7, height = 5, dpi = 500, bg = "white")


