library(tidyverse)
library(ggforce)
library(prismatic)
library(here)
library(glue)
library(patchwork)
library(ggtext)
source(here("R/plot_defaults.R"))
clr_genome = "gray50"

genomes <- c("arcgaz_anc_h2", "zalcal_v1")

scaf_order_18 <- list(c(1:18),
                      rev(c(17, 18, 16, 15, 14, 12, 11, 10, 13, 9, 8, 5, 6, 7, 4, 3, 2, 1))) |>
  set_names(nm = genomes)

read_size <- \(genome = "", y_base = 0, skip = 0, order_by = "size"){
  meta_order <- list(size = c(1:18),
                     name = c( 14, 9, 12, 4, 13, 10, 17, 2, 1,
                               15, 11, 7, 18, 5, 8, 3, 6, 16))[[order_by]]

  read_tsv(here::here("results", "genome", str_c(genome,".size")),
           col_names = c("chr", "size")) |>
    arrange(-size) |>
    mutate(size_ord = row_number(),
           size_meta = c(scaf_order_18[[genome]] , 19:length(size_ord))) |>
    arrange(size_meta) |>
    mutate(size_idx = c( meta_order[size_meta[1:18]], 19:length(size_ord))) |>
    arrange(size_idx) |>
    mutate(end_with_skip = cumsum(size + skip),
           start = lag(end_with_skip, default =  0),
           end = start + size,
           mid = (start + end) /2,
           eo = size_idx %% 2,
           y_base = y_base,
           genome = genome)
}

get_psl <- \(file){
  vroom::vroom(here::here("results","psl", file),
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
  pmap_dfr(read_size, skip = 1e7, order_by = "name")  |>
  mutate(label_sign = 2 * ( .5 - (genome == genomes[[1]])))

genome_summary <- sizes |>
  group_by(genome, y_base) |>
  summarise(n = n(),
            toal_length = max(end))

sizes_18 <- sizes |>
  filter(size_idx <= 18 )

psl <- get_psl(glue("slim_{genomes[[2]]}_on_{genomes[[1]]}.psl.gz"))

n_alignments <- 20
align_min_length <- 0.2e6

psl_filtered <- psl |>
  filter(tName %in% sizes_18$chr,
         qName %in% sizes_18$chr) |>
  arrange(-tSize) |>
  filter(tSize > align_min_length) |>
  # group_by(tName) |>
  # mutate(in_top_n_t = row_number() < n_alignments) |>
  # ungroup() |>
  # group_by(qName) |>
  # mutate(in_top_n_q = row_number() < n_alignments) |>
  # ungroup() |>
  # filter(in_top_n_t | in_top_n_q) |>
  left_join(sizes_18 |> select(tName = chr, tg_start = start, t_eo = eo)) |>
  left_join(sizes_18 |> select(qName = chr, qg_start = start, q_eo = eo)) |>
  mutate(ymin_start = tStart + tg_start,
         ymax_start = tEnd + tg_start,
         ymin_end = qStart + qg_start,
         ymax_end = qEnd + qg_start,
         group = row_number())

psl_diag <- psl_filtered |>
  select(ymin_start:group, dir = t_eo) |>
  pmap_dfr(diag_to_wide)

genome_width <- .04
skip <- .025

p1 <- psl_diag |>
  ggplot() +
  geom_diagonal_wide(aes(x = y, y = x, group = group,
                         color = factor(dir),
                         fill = after_scale(clr_lighten(color))),
                     orientation = "y", linewidth = .5) +
  geom_rect(data = sizes_18,
            aes(xmin = start, xmax = end,
                ymin = y_base + (label_sign * skip),
                ymax = y_base + (label_sign * (skip + genome_width))),
            color = clr_genome, fill = clr_lighten(clr_genome)) +
  geom_text(data = sizes_18,
            aes(x = mid,
                # y = y_base + (label_sign *  (skip + genome_width * .5)),
                y = y_base + (label_sign *  (2.2 * skip + genome_width)),
                label = chr,
                hjust = c(1, 0)[as.numeric(factor(label_sign))]),
            color = "black", # clr_genome,
            family = fnt_sel,
            angle = 90)  +
  # coord_polar() +
  # scale_x_continuous(limits = c(0,2.59e9)) +
  # scale_y_continuous(limits = c(-1.2,1.2)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.15, 1.4),
                  xlim = c(-2e7, max(sizes_18$end_with_skip)),
                  clip = "off") +
  scale_color_manual(values = clrs |> clr_alpha(.3)) +
  theme_void(base_family = fnt_sel) +
  theme(legend.position = "none")

p2 <- psl |>
  ggplot(aes(x = log10(tSize))) +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2]))+
  geom_vline(xintercept = log10(min(psl_filtered$tSize)),
             color = clrs[2], linetype = 3) +
  scale_x_continuous(glue("alignment size log<sub>10</sub>(bp)<br>{round(range(psl$tSize) / c(1, 1e6), digits = 2) %>% str_c(' ',.,c('bp ', 'Mb '), collapse = '–')}, n: {nrow(psl)}")) +
  coord_cartesian(expand = 0) +
  theme_minimal(base_family = fnt_sel) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = clr_genome),
        panel.grid = element_blank(),
        axis.title.x = element_markdown(lineheight = 1.5))

p3 <- psl_filtered |>
  ggplot(aes(x = tSize * 1e-6)) +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2])) +
  geom_vline(xintercept = min(psl_filtered$tSize) * 1e-6,
             color = clrs[2], linetype = 3) +
  scale_x_continuous(glue("alignment size (Mb)<br>{sprintf('%.2f',range(psl_filtered$tSize)/1e6) %>% str_c(' ',.,'Mb ', collapse = '–')}, n: {nrow(psl_filtered)}")) +
  coord_cartesian(expand = 0) +
  theme_minimal(base_family = fnt_sel) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = clr_genome),
        panel.grid = element_blank(),
        axis.title.x = element_markdown(lineheight = 1.5))

p1 +
  ( p2/ p3 + plot_layout(heights = c(1,1))) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(widths = c(1, .2)) &
  theme(text = element_text(family = fnt_sel))

scl <- .9
ggsave("img/arcgaz_a2_zalcal.pdf",
       width = 16 * scl,
       height = 5 * scl,
       device = cairo_pdf)

# ==============================
# scaf_pick <- 12:13
#
# psl_diag_sub <- psl_filtered |>
#   filter(tName %in% sizes_18$chr[sizes_18$size_idx %in% scaf_pick]) |>
#   select(ymin_start:group, dir = t_eo) |>
#   pmap_dfr(diag_to_wide)
#
# psl_diag_sub |>
#   ggplot() +
#   geom_diagonal_wide(aes(x = y, y = x, group = group,
#                          color = factor(dir),
#                          fill = after_scale(clr_lighten(color))),
#                      orientation = "y", linewidth = .5) +
#   geom_rect(data = sizes_18 |> filter(chr %in% sizes_18$chr[sizes_18$size_idx %in% scaf_pick]),
#             aes(xmin = start, xmax = end,
#                 ymin = y_base + (label_sign * skip),
#                 ymax = y_base + (label_sign * (skip + genome_width))),
#             color = clr_genome, fill = clr_lighten(clr_genome)) +
#   geom_text(data = sizes_18  |> filter(chr %in% sizes_18$chr[sizes_18$size_idx %in% scaf_pick]),
#             aes(x = mid, y = y_base + (label_sign *  (skip + genome_width * .5)), label = chr),
#             color = clr_genome, family = fnt_sel)  +
#   scale_color_manual(values = clrs |> clr_alpha(.3)) +
#   theme_void(base_family = fnt_sel) +
#   theme(legend.position = "none")
