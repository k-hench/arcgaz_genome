library(tidyverse)
library(prismatic)
library(patchwork)
library(here)
library(ggraph)
library(tidygraph)
library(glue)
library(ggtext)
library(ggnewscale)

source(here("R/plot_defaults.R"))
l_gerp_top <- readRDS(file = "~/Downloads/l_gerp.RDS")
l_fst <- readRDS(file = "~/Downloads/l_fst.RDS")

all_gerp_top_names <- l_gerp_top$tb |> as_tibble() |> filter(top_go) |> pluck("name")
all_gerp_names <- l_gerp_top$tb |> as_tibble() |> pluck("name")

all_fst_top_names <- l_fst$tb |> as_tibble() |> filter(top_go) |> pluck("name")
all_fst_names <- l_fst$tb |> as_tibble() |> pluck("name")

distances_to <- function(x){ function(...){ node_distance_to(x, mode = "all") } }

all_gerp_distances <- l_gerp_top$tb |>
  mutate(across(name,
                .fns = map(1:10, distances_to) |> set_names(nm = all_gerp_top_names),
                .names = "{.fn}")) |>
    as_tibble() |>
    filter(top_go) |>
    pivot_longer(starts_with("GO:"), names_to = "to", values_to = "distance") |>
  mutate(name = factor(name, levels = all_gerp_top_names),
         to = factor(to, levels = all_gerp_top_names))

all_fst_distances <- l_fst$tb |>
  mutate(across(name,
                .fns = map(1:10, distances_to) |> set_names(nm = all_fst_top_names),
                .names = "{.fn}")) |>
  as_tibble() |>
  filter(top_go) |>
  pivot_longer(starts_with("GO:"), names_to = "to", values_to = "distance") |>
  mutate(name = factor(name, levels = all_fst_top_names),
         to = factor(to, levels = all_fst_top_names))

p <- ggplot(mapping = aes(x = as.numeric(name),
                     y = as.numeric(to))) +
  geom_tile(data = all_gerp_distances |> filter(as.numeric(name) < as.numeric(to)),
            aes(fill = distance)) +
  geom_text(data = all_gerp_distances |>
              filter(as.numeric(name) < as.numeric(to),
                     distance <= 6),
            aes(label = distance), size = fnt_sz/ggplot2::.pt, color = "white", family = fnt_sel) +
  scale_fill_binned("Distance GERP GO",
                    limits = c(0, 12),
                    breaks = c(0,3,6,9,12),
                    type = \(...){scale_fill_stepsn(colors = c(clrs[1], "gray85"),...)},
                    guide = guide_colorsteps(title.position = "top",
                                             barheight = unit(5, "pt"))) +
  new_scale_fill() +
  geom_tile(data = all_fst_distances |> filter(as.numeric(name) > as.numeric(to)),
            aes(fill = distance)) +
  geom_text(data = all_fst_distances |> filter(as.numeric(name) > as.numeric(to),
                     distance <= 6),
            aes(label = distance), size = fnt_sz/ggplot2::.pt, color = "white", family = fnt_sel) +
  scale_fill_binned("Distance *F<sub>ST</sub>* GO",
                    limits = c(0, 12),
                    breaks = c(0,3,6,9,12),
                    type = \(...){scale_fill_stepsn(colors = c(clr_darken(clrs[2]),"gray85"),...)},
                    guide = guide_colorsteps(title.position = "top",
                                             barheight = unit(5, "pt"))) +
  scale_x_continuous(NULL,
                     breaks = 2:10, labels = all_fst_top_names[2:10],
                     sec.axis = sec_axis(trans = identity,
                                         breaks = 1:9,
                                         labels = all_gerp_top_names[1:9])) +
  scale_y_continuous(NULL,
                     breaks = 2:10, labels = all_gerp_top_names[2:10],
                     sec.axis = sec_axis(trans = identity,
                                         breaks = 1:9,
                                         labels = all_fst_top_names[1:9])) +
  coord_equal(expand = 0) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        axis.text.x.top = element_text(angle = 90, vjust = .5),
        legend.position = "bottom",
        legend.title = element_markdown())

ggsave(filename = here("results/img/go_distances.pdf"),
       width = 3, height = 3.5, device = cairo_pdf)
