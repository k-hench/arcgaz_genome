library(tidyverse)
library(prismatic)
library(cowplot)
library(ggfx)
source("R/plot_defaults.R")

clr_i <- 1

il_arcgaz <- hypoimg::hypo_read_svg("img/illustrations/arcgaz_ln.c.svg") |>
  hypoimg::hypo_recolor_svg(layer = 1, clr_lighten(clrs[clr_i], .75)) |>
  hypoimg::hypo_recolor_svg(layer = 2, clrs[clr_i])

ggdraw(with_shadow(il_arcgaz,
                   x_offset = 0, y_offset = 0,
                   sigma = 7, colour = clr_darken(clrs[clr_i])))

colorblindr::cvd_grid(il_arcgaz)
