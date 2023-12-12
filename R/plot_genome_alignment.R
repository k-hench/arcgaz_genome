library(tidyverse)
library(ggforce)
library(prismatic)
library(here)
library(glue)
library(patchwork)
library(ggtext)
library(plyranges)
library(rlang)

source(here("R/plot_defaults.R"))
clr_genome = "gray50"
genomes <- c("arcgaz_anc_h1", "zalcal_v1")

scaf_order_18 <- list(str_c("mscaf_a1_", c(str_pad(1:17, width = 2, pad = 0), "x")),
                      rev(c("NC_045612.1", "NC_045599.1", "NC_045610.1", "NC_045604.1", "NC_045595.1", "NC_045601.1",
                            "NC_045597.1", "NC_045605.1", "NC_045600.1", "NC_045596.1", "NC_045608.1", "NC_045606.1",
                            "NC_045611.1", "NC_045607.1", "NC_045598.1", "NC_045609.1", "NC_045602.1", "NC_045603.1"))) |>
  set_names(nm = genomes)

read_size <- \(genome = "", y_base = 0, skip = 0){
  read_tsv(here::here("results", "genome", str_c(genome,".size")),
           col_names = c("chr", "size")) |>
    filter(chr %in% scaf_order_18[[genome]]) |>
    mutate(chr_ord = factor(chr, levels = scaf_order_18[[genome]])) |>
    arrange(as.numeric(chr_ord)) |>
    mutate(end_with_skip = cumsum(size + skip),
           start = lag(end_with_skip, default =  0),
           end = start + size,
           mid = (start + end) /2,
           eo = row_number() %% 2,
           y_base = y_base,
           genome = genome)
}

get_psl <- \(file){
  vroom::vroom(here("results","psl", file),
               delim = "\t",
               col_names = c("matches", "misMatches", "repMatches", "nCount",
                             "qNumInsert", "qBaseInsert", "tNumInsert",
                             "tBaseInsert", "strand", "qName", "qSize", "qStart",
                             "qEnd", "tName", "tSize", "tStart", "tEnd",
                             "blockCount")) |>
    select( tName, tStart, tEnd, qName, qStart, qEnd, strand ) |>
    mutate(tSize = tEnd - tStart,
           qSize = qEnd - qStart)
}

diag_to_wide <- \(x_start = 0,
                  x_end = 1,
                  ymin_start, ymax_start,
                  ymin_end, ymax_end,
                  dir = "+",
                  group){
  tibble(`1_x_min` = x_start,
         `2_x_max` = x_end,
         `3_x_max` = x_end,
         `4_x_min` = x_start,
         `1_y_min` = ymin_start,
         `2_y_min` = ymin_end,
         `3_y_max` = ymax_end,
         `4_y_max` = ymax_start,
         group = group,
         dir = dir #c('-','+')[1 + ((ymin_start < ymax_start) == (ymin_end < ymax_end))]
  ) |>
    pivot_longer(cols = -c(group, dir)) |>
    separate(name,
             into = c("ord", "axis", "role"),
             convert = TRUE) |>
    arrange(ord) |>
    select(-role) |>
    pivot_wider(names_from = axis,
                values_from = value) |>
    select(-ord)
}

sizes <- tibble(genome = genomes,
                y_base = seq_along(genomes)-1) |>
  pmap_dfr(read_size, skip = 1e7)  |>
  mutate(label_sign = 2 * ( .5 - (genome == genomes[[1]])))

genome_summary <- sizes |>
  group_by(genome, y_base) |>
  summarise(n = n(),
            toal_length = max(end))

psl <- get_psl(glue("slim_{genomes[[2]]}_on_{genomes[[1]]}.psl.gz"))

align_min_length <- 0.2e6

psl_filtered <- psl |>
  filter(tName %in% sizes$chr,
         qName %in% sizes$chr) |>
  arrange(-tSize) |>
  filter(tSize > align_min_length) |>
  left_join(sizes |> select(tName = chr, tg_start = start, t_eo = eo)) |>
  left_join(sizes |> select(qName = chr, qg_start = start, q_eo = eo)) |>
  mutate(ymin_start = tStart + tg_start,
         ymax_start = tEnd + tg_start,
         ymin_end = qStart + qg_start,
         ymax_end = qEnd + qg_start,
         group = row_number())

psl_diag <- psl_filtered |>
  select(ymin_start:group, dir = t_eo) |>
  pmap_dfr(diag_to_wide)

genome_width <- .05
skip <- .025

il_arcgaz <- hypoimg::hypo_read_svg("img/illustrations/arcgaz_ln.c.svg") |>
  hypoimg::hypo_recolor_svg(layer = 1, "transparent")

il_zalcal <- hypoimg::hypo_read_svg("img/illustrations/zalcal_ln.c.svg") |>
  hypoimg::hypo_recolor_svg(layer = 1, "transparent")

il_height <- .5

p1 <- psl_diag |>
  ggplot() +
  annotation_custom(grob = il_arcgaz, xmin = -2.5e8, xmax = 3e7,
                    ymax = -.15, ymin = -.15 -il_height) +
  annotation_custom(grob = il_zalcal, xmin = -2.5e8, xmax = 3e7,
                    ymax = 1.1 + il_height, ymin = 1.1 ) +
  geom_diagonal_wide(aes(x = y, y = x, group = group,
                         color = factor(dir),
                         fill = after_scale(clr_lighten(color))),
                     orientation = "y", linewidth = plt_lwd)  +
  geom_rect(data = sizes,
            aes(xmin = start, xmax = end,
                ymin = y_base + (label_sign * skip),
                ymax = y_base + (label_sign * (skip + genome_width))),
            color = clr_genome, fill = clr_lighten(clr_genome,.75),
            linewidth = .5 * plt_lwd) +
  geom_text(data = sizes,
            aes(x = mid,
                # y = y_base + (label_sign *  (skip + genome_width * .5)),
                y = y_base + (label_sign *  (2.2 * skip + genome_width)),
                label = chr,
                hjust = c(1, 0)[as.numeric(factor(label_sign))]),
            color = "grey30",#"black", # clr_genome,
            family = fnt_sel,
            size = .3 * fnt_sz,
            angle = 90)  +
  geom_text(data = tibble(x = -2.5e7,
                          y = 1*c(-genome_width,genome_width) +  c(0,1),
                          label = c("A. gazella", "Z. californianus")),
            aes(x = x,
                y = y,
                label = label),
                hjust = 1,
            vjust = .5,
            color = "grey30",
            family = fnt_sel,
            size = .3 * fnt_sz,
            fontface = "italic")  +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.1, 1.35),
                  xlim = c(-2e8, max(sizes$end_with_skip)),
                  clip = "off") +
  scale_color_manual(values = clrs |> clr_alpha(.3)) +
  theme_void(base_family = fnt_sel) +
  theme(legend.position = "none"
        )

p2 <- psl |>
  ggplot(aes(x = log10(tSize))) +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2]),
               linewidth = plt_lwd)+
  geom_vline(xintercept = log10(min(psl_filtered$tSize)),
             color = clrs[2], linetype = 3, linewidth = plt_lwd) +
  scale_x_continuous(glue("Alignment Size log<sub>10</sub>(bp)<br>{round(range(psl$tSize) / c(1, 1e6), digits = 2) %>% str_c(' ',.,c('bp ', 'Mb '), collapse = '–')}, n: {nrow(psl)}")) +
  coord_cartesian(expand = 0, clip = "off") +
  theme_minimal(base_family = fnt_sel) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = clr_genome, linewidth = plt_lwd),
        panel.grid = element_blank(),
        axis.title.x = element_markdown(lineheight = 1.5))

p3 <- psl_filtered |>
  ggplot(aes(x = tSize * 1e-6)) +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2]),
               linewidth = plt_lwd) +
  geom_vline(xintercept = min(psl_filtered$tSize) * 1e-6,
             color = clrs[2], linetype = 3, linewidth = plt_lwd) +
  scale_x_continuous(glue("Alignment Size (Mb)<br>{sprintf('%.2f',range(psl_filtered$tSize)/1e6) %>% str_c(' ',.,'Mb ', collapse = '–')}, n: {nrow(psl_filtered)}")) +
  coord_cartesian(expand = 0) +
  theme_minimal(base_family = fnt_sel) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = clr_genome, linewidth = plt_lwd),
        panel.grid = element_blank(),
        axis.title.x = element_markdown(lineheight = 1.5))

p <- p1 + theme(plot.margin = margin(t = unit(12,"pt"))) +
  ( p2/ p3 + plot_layout(heights = c(1,1))) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(widths = c(1, .2)) &
  theme(text = element_text(family = fnt_sel,
                            size = fnt_sz),
        plot.tag = element_text(family = fnt_sel,
                                margin = margin(t = -unit(15, "pt"))))

ggsave(plot = p,
       filename = here("results/img/arcgaz_a1_zalcal.pdf"),
       width = fwidth,
       height = fwidth * 0.3125,
       device = cairo_pdf)

