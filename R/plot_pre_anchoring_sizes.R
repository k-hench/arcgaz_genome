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
             "arcgaz_dt_h2_hardmasked",
             "arcgaz_v3_hardmasked")

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

( p_sizes <- sizes |>
  filter(size_idx < n_largest) |>
  mutate(l18 = size_idx < 19,
         eo_18 = str_c(eo, as.numeric(l18))) |>
  ggplot() +
  geom_linerange(data = tibble(end = align_min_length),
                 aes(y = .9, xmin = 0, xmax = end ), size = 1.5) +
  geom_rect(data = genome_summary,
            aes(xmin = 0, xmax = toal_length,
                ymin = y_base, ymax = y_base + g_height),
            color = "gray50", fill = "gray90",
            linewidth = plt_lwd) +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = y_base, ymax = y_base + g_height,
                color = factor(eo_18),
                fill = after_scale(clr_lighten(color, .5))),
            linewidth = plt_lwd) +
  geom_richtext(data = genome_summary,
                aes(y = y_base + .5 * g_height,
                    label = glue("{str_replace(str_remove(genome, '_hardmasked'),'arcgaz_dt_h', 'HiRise h')}<br>n: {n}"),
                    x = -.1e9),
                hjust = 1, size = fnt_sz/2,
                fill = NA, label.color = NA, family = fnt_sel,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_y_continuous(NULL) +
  scale_x_continuous(NULL, breaks = c(0:3)*1e9,
                     labels = str_c(0:3, " Gb")) +
  scale_color_manual(values = c(clrs |> clr_lighten(.1) |> clr_desaturate(),
                                clrs |> clr_lighten(.6) |> clr_desaturate())|>
                       set_names(nm = c("01", "11", "00", "10")),
                     guide = "none") +
  coord_cartesian(xlim = c(-4.5e8, 2.6e9),
                  ylim = c(length(genomes),-.1),
                  expand = 0)+
  # coord_polar() +
  theme_minimal(base_family = fnt_sel, base_size = fnt_sz) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = .2),
        axis.ticks.length = unit(-2, "pt"),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) )


ggsave(plot = p_sizes,
       filename = here("results/img/fig_s_hirise_sizes.pdf"),
       width = fwidth,  height = fwidth/4,
       device = cairo_pdf)
