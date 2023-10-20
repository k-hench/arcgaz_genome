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

cov_from_psl <- \(prefix, genome){
  psl |>
    select(seqnames = str_c(prefix, "Name"),
           start =  str_c(prefix, "Start"),
           end =  str_c(prefix, "End"),
           strand) |>
    as_granges() |>
    compute_coverage() |>
    as_tibble() |>
    left_join(sizes_18 |> select(seqnames = chr, g_start = start, eo = eo))  |>
    mutate(gg_start = start + g_start,
           gg_end = end + g_start,
           genome = genome) |>
    arrange(gg_start)
}

sizes_to_granges <- \(genome_select){
  sizes |>
    filter(genome == genome_select) |>
    select(seqnames = chr, end = size) |>
    mutate(start = 1) |> as_granges()
}

genomes_granges <- c("arcgaz_anc_h2", "zalcal_v1") |> map(sizes_to_granges)
genome_windows <- genomes_granges |> map(\(g){slidingWindows(g, width = 10e4, step = 5e4) |> as_tibble() |>  mutate(window_id = row_number()) |> as_granges()})

data_coverage <- tibble(prefix = c("t", "q"), genome = c("arcgaz", "zalcal")) |>
  pmap_dfr(cov_from_psl)

  # bind_rows(cov_from_psl(tName, tStart, tEnd, "arcgaz"),
  #         cov_from_psl(qName, qStart, qEnd, "zalcal")) |

coverage_windows <- \(genome_idx){
  prefix = c("t", "q")[[genome_idx]]
  join_overlap_intersect(psl |>
                  select(seqnames = str_c(prefix, "Name"),
                         start =  str_c(prefix, "Start"),
                         end =  str_c(prefix, "End"),
                         strand) |>
                  as_granges() |>
                  compute_coverage(),
                as_granges(genome_windows[[genome_idx]])) |>
    as_tibble()  |>
    group_by(seqnames, window_id) |>
    summarise(start = min(start),
              end = max(end),
              avg_cov = sum(width * score) / sum(width)) |>
    ungroup() |>
    arrange(window_id) |>
    left_join(sizes_18 |> select(seqnames = chr, g_start = start, eo = eo))  |>
    mutate(gg_start = start + g_start,
           gg_end = end + g_start,
           gg_mid = (gg_start + gg_end) / 2,
           genome_idx = genome_idx,
           avg_cov_scl = avg_cov / max(avg_cov),
           avg_cov_25 = avg_cov / 5,
           avg_cov_25 = if_else(avg_cov_25 > 1, 1, avg_cov_25))
}

data_cov_windowed <- 1:2 |> map_dfr(coverage_windows)

genome_width <- .1
skip <- .025

il_arcgaz <- hypoimg::hypo_read_svg("img/illustrations/arcgaz_ln.c.svg") |>
  hypoimg::hypo_recolor_svg(layer = 1, "transparent")

il_zalcal <- hypoimg::hypo_read_svg("img/illustrations/zalcal_ln.c.svg") |>
  hypoimg::hypo_recolor_svg(layer = 1, "transparent")

# il_arcgaz <- grid::rasterGrob(png::readPNG("img/illustrations/argcaz_line.png"), interpolate = TRUE)
# il_zalcal <- grid::rasterGrob(png::readPNG("img/illustrations/zalcal_line.png"), interpolate = TRUE)
il_height <- .5

# data_cov_windowed |>
#   ggplot(aes()) +
#   geom_ribbon(data = data_cov_windowed,
#               aes(x = gg_mid, group = seqnames,
#                   ymin = 2 * (genome_idx -1),
#                   ymax = 2 * (genome_idx -1) + avg_cov_25),
#               linewidth = .3,
#               fill = "red") +
#   geom_hline(yintercept = c(0, 2) + (1/25))

p1 <- psl_diag |>
  ggplot() +
  annotation_custom(grob = il_arcgaz, xmin = -2.5e8, xmax = 3e7,
                    ymax = -.15, ymin = -.15 -il_height) +
  annotation_custom(grob = il_zalcal, xmin = -2.5e8, xmax = 3e7,
                    ymax = 1.02 + il_height, ymin = 1.02 ) +
  geom_diagonal_wide(aes(x = y, y = x, group = group,
                         color = factor(dir),
                         fill = after_scale(clr_lighten(color))),
                     orientation = "y", linewidth = plt_lwd)  +
  geom_linerange(data = tibble(y = c(0 - genome_width - skip,
                                     1 + skip) + (1/5) * genome_width,
                               start = -1e7,
                               end = c(max(sizes_18$end_with_skip[sizes_18$genome == "arcgaz_anc_h2"]),
                                       max(sizes_18$end_with_skip[sizes_18$genome == "zalcal_v1"]))),
                 aes(xmin = start, xmax = end,
                     y = y),
                 color = clr_genome,
                 linewidth = .5 * plt_lwd) +
  geom_rect(data = sizes_18,
            aes(xmin = start, xmax = end,
                ymin = y_base + (label_sign * skip),
                ymax = y_base + (label_sign * (skip + genome_width))),
            color = clr_genome, fill = clr_lighten(clr_genome,.75),
            linewidth = .5 * plt_lwd) +
    geom_ribbon(data = data_cov_windowed,
                aes(x = gg_mid, group = seqnames,
                    ymin = (genome_idx -1) + 2 * skip *(genome_idx - 1.5) + genome_width *(genome_idx - 2),
                    ymax = (genome_idx -1) + 2 * skip *(genome_idx - 1.5) + genome_width *(genome_idx - 2) +  genome_width * avg_cov_25),
                linewidth = plt_lwd,
                fill = clr_genome) +
    # geom_linerange(data = sizes_18,
    #           aes(xmin = start, xmax = end,
    #               y = y_base + (label_sign * (skip +  (label_sign == -1) * genome_width)) + genome_width * (1/5)),
    #           color = rgb(1,1,1,.5),
    #           linewidth = .5 * plt_lwd) +
  geom_text(data = sizes_18,
            aes(x = mid,
                # y = y_base + (label_sign *  (skip + genome_width * .5)),
                y = y_base + (label_sign *  (2.2 * skip + genome_width)),
                label = chr,
                hjust = c(1, 0)[as.numeric(factor(label_sign))]),
            color = "black", # clr_genome,
            family = fnt_sel,
            size = .35 * fnt_sz,
            angle = 90)  +
  # coord_polar() +
  # scale_x_continuous(limits = c(0,2.59e9)) +
  # scale_y_continuous(limits = c(-1.2,1.2)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.1, 1.35),
                  xlim = c(-2e8, max(sizes_18$end_with_skip)),
                  clip = "off") +
  scale_color_manual(values = clrs |> clr_alpha(.3)) +
  theme_void(base_family = fnt_sel) +
  theme(legend.position = "none")

p2 <- psl |>
  ggplot(aes(x = log10(tSize))) +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2]),
               linewidth = plt_lwd)+
  geom_vline(xintercept = log10(min(psl_filtered$tSize)),
             color = clrs[2], linetype = 3, linewidth = plt_lwd) +
  scale_x_continuous(glue("Alignment Size log<sub>10</sub>(bp)<br>{round(range(psl$tSize) / c(1, 1e6), digits = 2) %>% str_c(' ',.,c('bp ', 'Mb '), collapse = '–')}, n: {nrow(psl)}")) +
  coord_cartesian(expand = 0) +
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

p <- p1 +
  ( p2/ p3 + plot_layout(heights = c(1,1))) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1, .2)) &
  theme(text = element_text(family = fnt_sel,
                            size = fnt_sz))

# colorblindr::cvd_grid(p)
ggsave(plot = p,
       filename = "img/arcgaz_a2_zalcal.pdf",
       width = fwidth,
       height = fwidth * 0.3125,
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
