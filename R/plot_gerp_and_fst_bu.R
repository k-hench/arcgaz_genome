# input:
# - expand( "results/pinniped/gerp/tsv/gerp_snps_{scf}_summary.tsv.gz", scf = MSCFS )
# - expand( "results/pinniped/gerp/tsv/gerp_busco_{scf}_summary.tsv.gz", scf = MSCFS )
# - expand( "results/pinniped/fst/tsv/fst_snps_{scf}_summary.tsv.gz", scf = MSCFS )
# - expand( "results/pinniped/fst/tsv/fst_busco_{scf}_summary.tsv.gz", scf = MSCFS )
# - "data/genomes/arcgaz_anc_h1.genome"
# - "results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"
# - "results/fst_win.tsv.gz"
# output:
# - "results/pinniped/go_terms/thresholds.tsv"
# - "results/pinniped/busco_gerp_fst.tsv"
library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(ggtext)
library(ggrastr)
library(here)
source("R/plot_defaults.R")

fs <- 9
clr_l_gray <- rgb(.85,.85,.85)
scfs <- str_c("mscaf_a1_",
              c(str_pad(1:17, width = 2, pad =  0), "x"))

read_win_gerp <- \(scf){read_tsv(here(glue("results/pinniped/gerp/tsv/gerp_snps_{scf}_summary.tsv.gz")))}
read_busco_gerp <- \(scf){read_tsv(here(glue("results/pinniped/gerp/tsv/gerp_busco_{scf}_summary.tsv.gz")))}
read_win_fst <- \(scf){read_tsv(here(glue("results/pinniped/fst/tsv/fst_snps_{scf}_summary.tsv.gz")))}
read_busco_fst <- \(scf){read_tsv(here(glue("results/pinniped/fst/tsv/fst_busco_{scf}_summary.tsv.gz")))}

genome <- read_tsv(here("data/genomes/arcgaz_anc_h1.genome"),
                   col_names = c("chr", "length")) |>
  mutate(end = cumsum(length),
         start = lag(end, default = 0),
         mid = (start + end)/2,
         eo = row_number() %% 2)

plain_busco <- read_tsv(here("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"), skip = 2) |>
  filter(Status == "Complete") |>
  mutate(bstart = if_else(`Gene Start` < `Gene End`, `Gene Start`, `Gene End`),
         bend = if_else(`Gene Start` < `Gene End`, `Gene End`, `Gene Start`)) |>
  select(chr = Sequence,
         bstart, bend,
         busco_id = `# Busco id`)  |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = bstart + scaf_start,
         gend = bend + scaf_start,
         gmid = (gstart + gend)/2)

data_win_gerp <- scfs |>
  map_dfr(read_win_gerp) |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = start + scaf_start,
         gend = end + scaf_start,
         gmid = (gstart + gend)/2)

data_win_fst_vt <- read_tsv(here("results/fst_win.tsv.gz")) |>
  left_join(genome |> select(CHROM = chr, start = start, eo)) |>
  mutate(gstart = BIN_START + start,
         gend = BIN_END + start,
         gmid = (gstart + gend)/2,
         fst_w_0 = if_else(WEIGHTED_FST > 0, WEIGHTED_FST, 0),
         fst_m_0 = if_else(MEAN_FST > 0, MEAN_FST, 0))

data_win_fst <- scfs |>
  map_dfr(read_win_fst) |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = start + scaf_start,
         gend = end + scaf_start,
         gmid = (gstart + gend)/2)

data_win <- data_win_gerp |>
  select(chr, start, end, eo, gstart,gend, gmid, gerp_rs_mean, gerp_rs_med) |>
  left_join(data_win_fst_vt |>
              select(chr = CHROM, end = BIN_END, fst_w_0, fst_m_0) |>
              left_join(data_win_fst |>
                          select(chr, end, fst_mean, fst_med),
                        by = c("chr", "end")))

# fst_w_breaks <- quantile(data_win$fst_w_0,
#                          probs = c(.005, .995),
#                          na.rm = TRUE)

# ggplot(data_win) +
#   geom_point(aes(x = fst_m_0, y = fst_mean)) +
#   coord_equal()

data_buso_gerp <- scfs |>
  map_dfr(read_busco_gerp) |>
  left_join(plain_busco)

data_buso_fst <- scfs |>
  map_dfr(read_busco_fst)

data_busco <- data_buso_gerp |>
  select(busco_id, chr,
         bstart,   bend,
         gstart, gend, gmid,
         neutral_n_mean, neutral_n_med,
         gerp_rs_mean, gerp_rs_med) |>
  left_join(data_buso_fst |>
              select(busco_id,
                     n_snps_fst = n_snps,
                     fst_mean,
                     fst_med))

focal_gerp <- "gerp_rs_mean"
gerp_lab <- "average RS score (GERP)"
gerp_tag <- "mean"
focal_fst <- "fst_m_0"
fst_lab <- "mean"
focal_fst_b <- "fst_mean"

# focal_gerp <- "gerp_rs_mean"
# gerp_lab <- "average RS score (GERP)"
# gerp_tag <- "mean"
# focal_fst <- "fst_w_0"
# fst_lab <- "weighted_mean"
# focal_fst_b <- "fst_mean"

# focal_gerp <- "gerp_rs_med"
# gerp_lab <- "median RS score (GERP)"
# gerp_tag <- "median"
# focal_fst <- "fst_med"
# fst_lab <- "median"
# focal_fst_b <- "fst_med"

p1 <- ggplot() +
  geom_rect(data = genome,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = factor(eo))) +
  rasterize(geom_point(data = data_win,
                       aes(x = gmid, y = .data[[focal_gerp]],
                           color = factor(eo),
                           group = chr),
                       size = .15),
            dpi = 300) +
  geom_smooth(data = data_win,
              color = "white",
              linewidth = .5,
              span = 0.5,
              aes(x = gmid, y = .data[[focal_gerp]],
                  group = chr),
              se = FALSE) +
  scale_fill_manual(values = c(`1` = "white",
                               `0` = clr_l_gray),
                    guide = "none") +
  scale_color_manual(values = clr_alpha(clrs),
                     guide = "none") +
  scale_x_continuous(name = "Genomic Position",
                     labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(name = "Scaffold Id",
                                         trans = identity,
                                         breaks = genome$mid,
                                         labels = str_remove(genome$chr, "mscaf_a1_")))+
  labs(y = gerp_lab) +
  coord_cartesian(xlim = c(0, max(genome$end)),
                  expand = 0)+
  theme_ms(fontsize = fs)

p2 <- ggplot() +
  geom_rect(data = genome,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = factor(eo))) +
  rasterize(geom_point(data = data_win,
                       aes(x = gmid, y = .data[[focal_fst]],
                           color = factor(eo),
                           group = chr),
                       size = .15),
            dpi = 300) +
  geom_smooth(data = data_win,
              color = "white",
              linewidth = .5,
              span = 0.5,
              aes(x = gmid,
                  y = .data[[focal_fst]],
                  group = chr),
              se = FALSE) +
  scale_fill_manual(values = c(`1` = "white",
                               `0` = clr_l_gray),
                    guide = "none") +
  scale_color_manual(values = clr_alpha(clrs),
                     guide = "none") +
  scale_x_continuous(name = "Genomic Position", labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(trans = identity,
                                         breaks = genome$mid,
                                         labels = str_remove(genome$chr, "mscaf_a1_")))+
  scale_alpha_manual(values = c(`TRUE` = 6,
                                `FALSE` = 1),
                     guide = "none") +
  labs(y = glue("*F<sub>ST</sub>* ({fst_lab})")) +
  coord_cartesian(xlim = c(0, max(genome$end)),
                  expand = 0)+
  theme_ms(fontsize = fs)+
  theme(axis.title.y = element_markdown())

pb <- ggplot() +
  geom_rect(data = genome,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = factor(eo))) +
  rasterize(geom_point(data = plain_busco,
                       aes(x = gmid, y = 0, color = factor(eo)),
                       shape = "|", alpha = .3),
            dpi = 300) +
  scale_fill_manual(values = c(`1` = "white",
                               `0` = clr_l_gray),
                    guide = "none") +
  scale_color_manual(values = clr_alpha(clrs),
                     guide = "none") +
  scale_x_continuous(name = "Genomic Position (Gb)", labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(trans = identity,
                                         breaks = genome$mid,
                                         labels = str_remove(genome$chr, "mscaf_a1_")))+
  scale_alpha_manual(values = c(`TRUE` = 6,
                                `FALSE` = 1),
                     guide = "none") +
  labs(y = "BUSCOs") +
  coord_cartesian(xlim = c(0, max(genome$end)),
                  ylim = c(-1,1),
                  expand = 0)+
  theme_void(base_family = fnt_sel, base_size = fs)

pp1 <- p1 +
  theme(axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()) +
  pb +
  p2 +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) +
  plot_layout(ncol = 1, heights = c(1,.1,1))

p3 <- data_win |>
  ggplot(aes(x = .data[[focal_gerp]],
             y = .data[[focal_fst]])) +
  geom_hex(aes(fill = after_stat(log10(count)),
               color = after_scale(clr_darken(fill,.2))),
           linewidth = .2,
           bins = 50) +
  scale_fill_gradientn(colours = c(clr_l_gray,
                                   clrs[2],
                                   clrs[1]))+
  labs(x = gerp_lab,
       y = glue("*F<sub>ST</sub>* ({fst_lab})")) +
  guides(fill = guide_colorbar(title.position = "left",
                               barwidth = unit(5, "pt"),
                               barheight = unit(100, "pt"))) +
  theme_ms(fontsize = fs) +
  theme(legend.position = c(0.01,0.01),
        legend.justification = c(0,0),
        legend.title = element_text(angle = 90, hjust = .5),
        axis.title.y = element_markdown())

pp2 <- p3 + pp1 +
  plot_layout(widths = c(.5, 1)) +
  plot_annotation(#tag_levels = "a",
    tag_levels = list(c('a', 'b', "", "c")),
    theme = theme(plot.tag = element_text(family = fnt_sel),
                  plot.tag.position = c(-.1,.5))) &
  theme(plot.tag = element_text(family = fnt_sel),
        plot.tag.position = c(-.01, 1))

ggsave(plot = pp2,
       filename = glue("~/Downloads/arcgaz_genome_gerp_{gerp_tag}_fst_{fst_lab}_whg.pdf"),
       width = 10, height = 5, device = cairo_pdf)

# within BUSCO comparison
querry_thresholds <- \(stat, type, by = "BUSCO"){
  dat <- list(BUSCO  = data_busco,
              win = data_win)
  prbs <- list(gerp = 1 - c(.95, .975, .995, .001),
               fst = c(.95, .975, .995))
  trs <- quantile(dat[[by]][[stat]],
           probs = prbs[[type]],
           na.rm = TRUE)
  tibble(stat = stat,
         type = type,
         threshold = trs,
         prob = names(trs),
         by = by)
}

all_thresholds <- tibble(stat = rep(c("gerp_rs_mean", "fst_mean"), 2),
       type = rep(c("gerp", "fst"),2),
       by = rep(c("BUSCO", "win"), each = 2)) |>
  pmap_dfr(querry_thresholds)

# all_thresholds |>
#   write_tsv(here("results/pinniped/go_terms/thresholds.tsv"))

quantile(data_busco[[focal_gerp]],
         probs = 1 - c(.95, .975, .995, .00000000000001),
         na.rm = TRUE)


gerp_threshold <- quantile(data_win[[focal_gerp]],
                          probs = 1 - c(.95, .975, .995, .00000000000001),
                          na.rm = TRUE)

fst_threshold <- quantile(data_win[[focal_fst]],
                         probs = c(.95, .975, .995),
                         na.rm = TRUE)

range(data_busco$gerp_rs_mean)

p4 <- data_busco |>
  ggplot(aes(x = .data[[focal_gerp]])) +
  geom_vline(xintercept = gerp_threshold,
             color = "darkgray",
             linewidth = .3,
             linetype = c(1,2,3,3)) +
  geom_histogram(boundary = 0,
                 color = clrs[[2]],
                 fill = clr_alpha(clrs[[2]], .6),
                 linewidth = plt_lwd,
                 binwidth = .001) +
  geom_vline(xintercept = mean(data_busco[[focal_gerp]]),
             color = "black", linetype = 3) +
  coord_cartesian(expand = 0,
                  xlim = c(1.1, 1.05) * range(data_busco[[focal_gerp]])) +
  theme_ms(fontsize = fs) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

p5 <- data_busco |>
  ggplot(aes(x = .data[[focal_fst_b]])) +
  geom_vline(xintercept = fst_threshold,
             color = "darkgray",
             linewidth = .3,
             linetype = c(1,2,3)) +
  geom_histogram(boundary = 0,
                 color = clrs[[2]],
                 fill = clr_alpha(clrs[[2]], .6),
                 linewidth = plt_lwd,
                 binwidth = .015) +
  geom_vline(xintercept = mean(data_busco[[focal_fst_b]], na.rm = TRUE),
             color = "black", linetype = 3)  +
  coord_flip(expand = 0,
             xlim = c(-.025,1.025)) +
  theme_ms(fontsize = fs) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

p6 <- data_busco |>
  ggplot(aes(x = .data[[focal_gerp]],
             y = .data[[focal_fst_b]])) +
  geom_vline(xintercept = gerp_threshold,
             color = "darkgray",
             linewidth = .3,
             linetype = c(1,2,3, 3)) +
  geom_hline(yintercept = fst_threshold,
             color = "darkgray",
             linewidth = .3,
             linetype = c(1,2,3)) +
  geom_hex(aes(fill = after_stat(log10(count)),
               color = after_scale(clr_darken(fill,.2))),
           linewidth = .2,
           bins = 100) +
  ggdensity::geom_hdr_lines(probs = c(0.95, 0.66),
                            aes(linewidth = after_stat(probs),
                                linetype = after_stat(probs)),
                            alpha = 1) +
  scale_fill_gradientn(colours = c(clr_l_gray,
                                   clrs[2],
                                   clrs[1])) +
  scale_linetype_manual(values = c(`95%` = 3,
                                   `66%` = 1),
                        guide = "none") +
  scale_linewidth_manual(values = c(`95%` = .5,
                                    `66%` = .25),
                         guide = "none") +
  labs(x = gerp_lab,
       y = glue("*F<sub>ST</sub>* ({fst_lab})")) +
  guides(fill = guide_colorbar(title.position = "left",
                               barwidth = unit(5, "pt"),
                               barheight = unit(100, "pt"))) +
  coord_cartesian(xlim = c(1.1, 1.05) * range(data_busco[[focal_gerp]]),
                  ylim = c(-.025, 1.025),
                  expand = 0) +
  theme_ms(fontsize = fs) +
  theme(legend.position = c(0.01,0.01),
        legend.justification = c(0,0),
        legend.box = "horizontal",
        legend.title = element_text(angle = 90, hjust = .5),
        axis.title.y = element_markdown())

pp3 <- p4 + plot_spacer() +
  p6 + p5 +
  plot_layout(widths = c(1,.2), heights = c(.2, 1)) +
  plot_annotation(subtitle = "Within BUSCO comparison") &
  theme(plot.subtitle = element_markdown(family = fnt_sel, hjust = .5),
        plot.tag = element_text(family = fnt_sel))

ggsave(plot = pp3,
       filename = glue("~/Downloads/arcgaz_genome_busco_gerp_{gerp_tag}_vs_fst_{fst_lab}.pdf"),
       width = 5, height = 5, device = cairo_pdf)

# data_busco |>
#    mutate(fst_top_25 = fst_mean > fst_threshold[["97.5%"]],
#           gerp_low_25 = gerp_rs_mean < gerp_threshold[["2.5%"]],
#           gerp_top_005 = gerp_rs_mean < gerp_threshold[["99.995%"]]) |>
#   write_tsv(here("results/pinniped/busco_gerp_fst.tsv"))

# Fst median mean comparison

p7 <- data_win |>
  mutate(fst_med_0 = if_else(fst_med > 0, fst_med, 0)) |>
  ggplot(aes(x = fst_m_0,
             y = fst_med_0))  +
  labs(subtitle = "50kb windows")

p8 <- data_busco|>
  ggplot(aes(x = fst_mean,
             y = fst_med)) +
  labs(subtitle = "BUSCOS")

pp4 <- p7 + p8 +
  plot_annotation(tag_levels = "a",
                  subtitle = "*F<sub>ST</sub>* comparison of mean and median") &
  geom_abline(slope = 1, intercept = 0,
              linewidth = .3, color = "lightgray") &
  geom_hex(aes(fill = after_stat(log10(count)),
               color = after_scale(clr_darken(fill,.2))),
           linewidth = .2,
           bins = 100) &
    ggdensity::geom_hdr_lines(probs = c(0.95, 0.66),
                            aes(linewidth = after_stat(probs),
                                linetype = after_stat(probs)),
                            alpha = 1) &
  scale_fill_gradientn(colours = c(clr_l_gray,
                                   clrs[2],
                                   clrs[1])) &
  scale_linetype_manual(values = c(`95%` = 3,
                                   `66%` = 1),
                        guide = "none") &
  scale_linewidth_manual(values = c(`95%` = .3,
                                    `66%` = .25),
                         guide = "none") &
  guides(fill = guide_colorbar(title.position = "left",
                               barwidth = unit(5, "pt"),
                               barheight = unit(100, "pt"))) &
  coord_equal() &
  theme_ms(fontsize = fs) &
  theme(legend.position = c(0.01,0.99),
        legend.justification = c(0,1),
        legend.box = "horizontal",
        legend.title = element_text(angle = 90, hjust = .5),
        axis.title.y = element_markdown())  &
  theme(plot.subtitle = element_markdown(family = fnt_sel, hjust = .5),
        plot.tag = element_text(family = fnt_sel))

ggsave(plot = pp4,
       filename = glue("~/Downloads/arcgaz_genome_fst_median_vs_mean.pdf"),
       width = 9, height = 4.5, device = cairo_pdf)
