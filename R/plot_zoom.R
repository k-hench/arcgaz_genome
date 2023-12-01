# input:
# - "data/genomes/arcgaz_anc_h1.genome"
# - "results/pinniped/go_terms/thresholds.tsv"
# - "results/pinniped/busco_gerp_fst.tsv"
# - "results/pinniped/win_gerp_fst.tsv"
# - "results/pinniped/win_outlier_summary.tsv"
# output:
# -
library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(ggtext)
library(here)
library(arcgazgen)
source("R/plot_defaults.R")

fs <- 12#fnt_sz
clr_l_gray <- rgb(.85,.85,.85)

scfs <- str_c("mscaf_a1_",
              c(str_pad(1:17, width = 2, pad = 0), "x"))

read_10k <- \(scf, type = "gerp", win = "GERP"){
  read_tsv(here(glue("results/pinniped/{type}/tsv/fine/{type}_snps_{scf}_w10k_s5k_summary.tsv.gz"))) |>
    mutate(mid = (start + end)/2,
           window = win)
}

data_gerp_10k <- map_dfr(scfs, read_10k)
data_fst_10k <- map_dfr(scfs, read_10k, type = "fst", win = "*F<sub>ST</sub>*")
# data_fst_bp <- read_tsv(here("results/fst_bp.tsv.gz")) |>
#   filter(!is.na(WEIR_AND_COCKERHAM_FST))

data_cov <- scfs |>
  map(\(scf){read_tsv(here(glue("results/neutral_tree/cov/{scf}.collapsed.bed.gz")),
                                col_names = c('chr', "start", "end", "cov"))}) |>
  set_names(nm = scfs)

data_win  <- read_tsv(here("results/pinniped/win_gerp_fst.tsv"))
data_outlier <- read_tsv(here("results/pinniped/win_outlier_summary.tsv"))
data_busco_go_summary <- read_tsv(here("results/pinniped/go_terms/go_term_busco_stats.tsv")) |>
  group_by(busco_id) |>
  summarise(busco_type = case_when(
    any(gerp_top_rank <= 10) & any(fst_rank <= 10) ~ 3,
    any(fst_rank <= 10) ~ 1,
    any(gerp_top_rank <= 10) ~ 2,
    !(any(gerp_top_rank <= 10) | any(fst_rank <= 10)) ~ 4,
    .default = NA
  )) |>
  ungroup()

data_busco <- read_tsv(here("results/pinniped/busco_gerp_fst.tsv")) |>
  left_join(data_busco_go_summary) |>
  mutate(busco_type = replace_na(busco_type, 4),
         busco_color = c(clrs, "black", "lightgray")[busco_type])

zoom_buffer <- 6e4

plot_zoom <- \(chr, start, end,
               z_buff = zoom_buffer,
               n_tracks = 5,
               n_tracks_b = 5,
               focal_genes = c(),
               # b_track_height = .5,
               outlier_label = "",
               ...){
  z_start <- start - z_buff
  z_end <- end + z_buff
  chr_in <- chr

  data_gerp <- data_gerp_10k |>
    filter(chr == chr_in,
           start < z_end & end > z_start)

  data_fst <- data_fst_10k |>
    filter(chr == chr_in,
           start < z_end & end > z_start)

  data_fst_bp_z <- data_fst_bp |>
    filter(CHROM == chr_in,
           POS < z_end & POS > z_start) |>
    mutate(window = "*F<sub>ST</sub> (bp)*")


  data_cov_z <- data_cov[[chr_in]]|>
    filter(end > z_start & start < z_end) |>
    pivot_longer(start:end, values_to = "pos") |>
    mutate(window = "coverage")
    #   mutate(y_track = (dplyr::row_number() - 1) %% n_tracks_b) |>
    #   rowwise()
  # data_busco_z <- data_busco |>
  #   filter(chr == chr_in) |>
  #   filter(bend > z_start & bstart < z_end) |>
  #   mutate(y_track = (dplyr::row_number() - 1) %% n_tracks_b) |>
  #   rowwise() |>
  #   dplyr::mutate(bstart = max(bstart, z_start),
  #                 bend = min(bend, z_end)) |>
  #   ungroup() |>
  #   mutate(mid = (start  + end )/2,
  #          window = "BUSCO") |>
  #   select(busco_id, chr, bstart, bend, y_track, busco_type, busco_color, window) |>
  #   arrange(bstart)

  gerp_range <- range(data_gerp_10k$gerp_rs_mean)
  # fst_range <- range(data_fst$fst_mean)
  # fst_range <- range(data_fst_bp_z$WEIR_AND_COCKERHAM_FST)
  fst_range <- 0:1

  win_types <- c("genes",# "BUSCO",
                 "GERP", "*F<sub>ST</sub>*", "*F<sub>ST</sub> (bp)*", "coverage")
  p1 <- ggplot() +
    geom_blank(data = tibble(window = rep(win_types[1], each = 2),
                             y = c(c(-0.5, n_tracks - 0.5)#,
                                   # c(-0.5, n_tracks_b - 0.5)
                                   )),
                             aes(x = -Inf, y = y)) +
                 geom_rect(data = tibble(x1 = start, x2 = end),
                           aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
                           color = "transparent",
                           fill = clr_alpha("black", .1)) +
                 # geom_rect(data = data_busco_z,
                 #           aes(xmin = bstart,
                 #               xmax = bend,
                 #               ymin = y_track - b_track_height/2,
                 #               ymax = y_track + b_track_height/2,
                 #               color = busco_color,
                 #               fill =  busco_color,
                 #               group = busco_id)) +
                 scale_fill_identity() +
                 scale_color_identity() +
                 ggnewscale::new_scale_fill() +
                 ggnewscale::new_scale_color() +
                 ag_plot_zoom(z_chr = chr,
                              z_start = z_start, z_end = z_end,
                              n_tracks = n_tracks,
                              style = "simple",
                              focal_genes = focal_genes,
                              label_offset = -0.7,
                              label_size = .5*fs/ggplot2::.pt,
                              label_fun = \(lb){ if_else(str_length(lb) > 10,
                                                         str_c(str_sub(lb, 1, 8), ".."),
                                                         lb) }) +
                 facet_grid(factor(window, levels = win_types) ~ .,
                            scales = "free", switch = "y") +
                 coord_cartesian(xlim = c(z_start, z_end),
                                 expand = 0,
                                 clip = "off")+
                 labs(subtitle = outlier_label) +
                 theme_ms(fontsize = fs) +
                 theme(strip.text.y.left = element_markdown(),
                       axis.text = element_blank(),
                       axis.line = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title = element_blank())

               p2 <- ggplot() +
                 geom_blank(data = tibble(window = rep(win_types[2:5], each = 2),
                                          y = c(gerp_range + c(-.1,.1) * diff(gerp_range),
                                                fst_range + c(-.1,.1),
                                                fst_range + c(-.1,.1),
                                                c(-.5,11.5))),
                            aes(x = -Inf, y = y)) +
                 geom_rect(data = tibble(x1 = start, x2 = end),
                           aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
                           color = "transparent",
                           fill = clr_alpha("black", .1)) +
                 geom_line(data = data_gerp,
                           aes(x = mid, y = gerp_rs_mean),
                           linewidth = .3, color = clrs[[1]]) +
                 geom_point(data = data_fst_bp_z,
                           aes(x = POS, y = WEIR_AND_COCKERHAM_FST),
                           size = .3, color = clr_alpha(clrs[[2]],.15)) +
                 geom_line(data = data_fst,
                           aes(x = mid, y = fst_mean),
                           linewidth = .3, color = clrs[[2]]) +
                 geom_step(data = data_cov_z,
                           aes(x = pos, y = cov),
                           linewidth = .3) +
                 scale_x_continuous(labels = \(x){sprintf("%.1f", x *1e-6)}) +
                 facet_grid(factor(window, levels = win_types) ~ .,
                            scales = "free", switch = "y") +
                 coord_cartesian(xlim = c(z_start, z_end),
                                 expand = 0) +
                 labs(x = glue("Position (Mb, {str_remove(chr,'mscaf_a1_')})")) +
                 theme_ms(fontsize = fs) +
                 theme(strip.text.y.right = element_markdown(),
                       axis.title.y = element_blank())

               p <- p1 / p2 +
                 plot_layout(heights = c(.5, 1)) &
                 theme(strip.placement = "outside",
                       plot.subtitle = element_text(hjust = .5))

               list(p = p,
                    # z_busco = data_busco_z,
                    data_genes = ag_gene_data(z_chr = chr, z_start = z_start, z_end = z_end, n_tracks = n_tracks))
}

p_list <- data_outlier[c(4,8,9,10,11,15,16,17,19,20,
                         21,22,23,24,25,26,27,28,29,
                         30,31,32,33,34,35),] |>
  pmap(plot_zoom, focal_genes = c("HOXA5"))

# p_list[c(4,8,9,10,11,15,16,17,19,20,
#          21,22,23,24,25,26,27,28,29,
#          30,31,32,33,34,35)] |>
p_list |>
  map(\(x){x$p}) |>
  wrap_plots(nrow = 4, guides = "collect")

p_list[[1]]$p


random_buscos <- data_busco |>
  select(chr, start = bstart, end = bend, outlier_label = busco_id) |>
  slice_sample(n = 16) |>
  pmap(plot_zoom)

random_buscos |>
  map(\(x){x$p}) |>
  wrap_plots(nrow = 4, guides = "collect")

top_gerp_buscos <- data_busco |>
  arrange(-gerp_rs_mean) |>
  select(chr, start = bstart, end = bend, outlier_label = busco_id, gerp_rs_mean) |>
  slice_head(n = 16) |>
  pmap(plot_zoom)

top_gerp_buscos |>
  map(\(x){x$p}) |>
  wrap_plots(nrow = 4, guides = "collect")

top_fst_buscos <- data_busco |>
  arrange(-fst_mean) |>
  select(chr, start = bstart, end = bend, outlier_label = busco_id, fst_mean) |>
  slice_head(n = 16) |>
  pmap(plot_zoom)

top_fst_buscos |>
  map(\(x){x$p}) |>
  wrap_plots(nrow = 4, guides = "collect")


data_busco |>
  mutate(length = bend - bstart) |>
  select(length, fst_mean, gerp_rs_mean) |>
  pivot_longer(-length) |>
  filter(!is.na(value)) |>
  group_by(name) |>
  # mutate(value = rethinking::standardize(value)) |>
  ungroup() |>
  ggplot(aes(x = log10(length),#*1e-3,
             y = value)) +
  ggdensity::geom_hdr() +
  facet_wrap(name ~ ., scales = "free") +
  theme_ms(fontsize = 13)
