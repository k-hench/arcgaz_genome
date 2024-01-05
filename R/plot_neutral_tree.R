# input:
# - "results/neutral_tree/rerooted.tree"
# output:
# - "results/img/neutral_tree.pdf"
library(tidyverse)
library(ggtree)
library(prismatic)
library(here)
library(ggtext)
source("R/plot_defaults.R")

tree <- read.tree(file = here("results/neutral_tree/rerooted.tree"))

# cumulative branch length (as reference for GERP RS scores)
sum(tree$edge.length)

spec_long <- c(odoros = "*Odobenus rosmarus*",
               calurs = "*Callorhinus ursinus*",
               arcgaz = "__*Arctocephalus gazella*__",
               eumjub = "*Eumetopias jubatus*",
               zalcal = "*Zalophus californianus*",
               lepwed = "*Leptonychotes weddellii*",
               mirleo = "*Mirounga leonina*",
               mirang = "*Mirounga angustirostris*",
               neosch = "*Neomonachus schauinslandi*",
               phovit = "*Phoca vitulina*",
               halgry = "*Halichoerus grypus*")

gens <- c("otariidae", "phocidae")

# clrs <- rcartocolor::carto_pal(5, "Prism")[c(1, 4)] |>
#   clr_lighten() |>
#   clr_alpha() |>
#   color()

specs_short <- str_replace(spec_long, "([A-Z])[a-z]* ([[a-z]*])", "\\1. \\2") |>
  set_names(nm = names(spec_long))

y_offfset <- .1
i_h <- 2.4/2

ott <- hypoimg::hypo_read_svg(here("results/img/illustrations/arcgaz.c.svg")) |>
  hypoimg::hypo_recolor_svg(color = "white")
pho <- hypoimg::hypo_read_svg(here("results/img/illustrations/halgry.c.svg")) |>
  hypoimg::hypo_recolor_svg(color = "white")
odo <- hypoimg::hypo_read_svg(here("results/img/illustrations/odoros.c.svg")) |>
  hypoimg::hypo_recolor_svg(color = clr_alpha("black"))

p <- ggtree(tree, color = NA) +
  geom_rect(inherit.aes = FALSE,
            data = tibble(x1 = c(.00725, .0232),
                          x2 = c(.03, .0465),
                          y1 = c(1.5, 5.5) + y_offfset/2,
                          y2 = c(5.5, 11.5) - y_offfset/2,
                          group = gens),
            aes(xmin = x1, xmax = x2,
                ymin = y1, ymax = y2,
                fill = group)) +
  annotation_custom(ott,
                    xmin = .0195,
                    xmax = .028,
                    ymin = 3.5 - i_h,
                    ymax = 3.5 + i_h) +
  annotation_custom(pho,
                    xmin = .0355,
                    xmax = .0465,
                    ymin = 6.2 - i_h,
                    ymax = 6.2 + i_h) +
  annotation_custom(odo,
                    xmin = .016,
                    xmax = .0225,
                    ymin = 1 - i_h/2.5,
                    ymax = 1 + i_h/2.5) +
  geom_text(data = tibble(x = c(.029, .0455),
                          y = c(3.5, 8.5),
                          label = gens),
            aes(label = label),
            angle = -90, size = 5, fontface = "bold", color = "white") +
  geom_tree(size = .4) +
  # geom_tiplab(aes(x = x + .0003,
  #                 label = specs_short[label]),
  #             fontface = "italic") +
  geom_richtext(aes(x = x + .0003,
                    label = specs_short[label]),
                hjust = 0,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_treescale(offset = -.4,
                 width = 0.005,
                 x = 0,
                 y = 8.5,
                 linesize = .4,
                 label = "subs. / site",
                 offset.label = .45)  +
  scale_fill_manual(values = clrs |>
                      clr_lighten(.2) |>
                      clr_alpha(.7),
                    guide = "none") +
  coord_cartesian(xlim = c(0, .0465),
                  ylim = c(0.75, 11.5),
                  expand = 0)

ggsave(here("results/img/neutral_tree.pdf"),
       plot = p,
       width = 8, height = 6,
       device = cairo_pdf)

