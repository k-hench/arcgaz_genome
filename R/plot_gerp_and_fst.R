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
library(ggrastr)
library(here)
library(arcgazgen)
source("R/plot_defaults.R")

fs <- fnt_sz
clr_l_gray <- rgb(.85,.85,.85)

data_win  <- read_tsv(here("results/pinniped/win_gerp_fst.tsv"))
data_outlier <- read_tsv(here("results/pinniped/win_outlier_summary.tsv"))
data_busco <- read_tsv(here("results/pinniped/busco_gerp_fst.tsv"))
data_busco_go_stats <- read_tsv(here("results/pinniped/go_terms/go_term_busco_stats.tsv"))

gerp_lab <- "Average GERP RS Score"
fst_lab <- "Average *F<sub>ST</sub>*"

p1 <- ggplot() +
  geom_ag_genome()+
  geom_linerange(data = data_outlier |> filter(outlier_type == "GERP"),
                 aes(x = gmid, ymin = .0775, ymax = Inf),
                 linewidth = plt_lwd,
                 color = "black") +
  rasterize(geom_point(data = data_win,
                       aes(x = gmid, y = gerp_rs_mean,
                           color = factor(eo),
                           group = chr),
                       size = .01,
                       shape = 19),
            dpi = 300, dev = "ragg") +
    geom_text(data = data_outlier |> filter(outlier_type == "GERP"),
                 aes(x = gmid -3e7*c(1,-1,1,1,1,-1,1,1),
                     y = .0775,
                     label = outlier_label),
                 color = "black",
                 size = .9 * fs / ggplot2::.pt,
                 family = fnt_sel) +
  geom_smooth(data = data_win,
              color = "white",
              linewidth = .2,
              span = 0.5,
              aes(x = gmid, y = gerp_rs_mean,
                  group = chr),
              se = FALSE) +
  scale_color_manual(values = c(clrs[[1]], clr_lighten(clrs[[1]],.33)),
                     guide = "none") +
  scale_x_continuous(name = "Genomic Position",
                     labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(name = "Scaffold Id",
                                         trans = identity,
                                         breaks = arcgaz_genome$mid,
                                         labels = str_remove(arcgaz_genome$chr, "mscaf_a1_")))+
  labs(y = gerp_lab) +
  coord_cartesian(xlim = c(0, max(arcgaz_genome$end)),
                  ylim = c(-.01, .0825),
                  expand = 0)+
  theme_ms(fontsize = fs)

p2 <- ggplot() +
  geom_ag_genome()+
  geom_linerange(data = data_outlier |>
                   filter(outlier_type == "fst"),
                 aes(x = gmid, ymin = 1, ymax = Inf),
                 linewidth = plt_lwd,
                 color = "black") +
  rasterize(geom_point(data = data_win,
                       aes(x = gmid, y = fst_m_0,
                           color = factor(eo),
                           group = chr),
                       size = .01,
                       shape = 19),
            dpi = 300, dev = "ragg") +
  geom_text(data = data_outlier |>
              filter(outlier_type == "fst") |>
              filter(row_number() != 6),
            aes(x = gmid -3e7  * rep(c(1,-1), c(5,1)),
                y = 1,
                label = outlier_label),
            color = "black",
            size = .9 * fs / ggplot2::.pt,
            family = fnt_sel) +
  geom_smooth(data = data_win,
              color = "white",
              linewidth = .2,
              span = 0.5,
              aes(x = gmid,
                  y = fst_mean,
                  group = chr),
              se = FALSE) +
  scale_color_manual(values = c(clrs[[2]], clr_lighten(clrs[[2]],.33)),
                     guide = "none") +
  scale_x_continuous(name = "Genomic Position", labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(trans = identity,
                                         breaks = arcgaz_genome$mid,
                                         labels = str_remove(arcgaz_genome$chr, "mscaf_a1_")))+
  scale_alpha_manual(values = c(`TRUE` = 6,
                                `FALSE` = 1),
                     guide = "none") +
  labs(y = fst_lab) +
  coord_cartesian(xlim = c(0, max(arcgaz_genome$end)),
                  ylim = c(0, 1.05),
                  expand = 0)+
  theme_ms(fontsize = fs)+
  theme(axis.title.y = element_markdown())

l_hight <- .75
pb <- ggplot() +
  geom_ag_genome()+
  rasterize(geom_linerange(data = data_busco,
                       aes(x = gmid, ymin = -.5*l_hight, ymax = .5 * l_hight),
                       color = "darkgray",
                       linewidth = .1,
                       alpha = .3),
            dpi = 300, dev = "ragg") +
  geom_linerange(data = data_busco_go_stats |>
                   filter(gerp_top_rank <= 10) |>
                   filter(!duplicated(busco_id)),
                 aes(x = gmid,
                     ymin = 1 -.5*l_hight,
                     ymax = 1 +.5 * l_hight),
                 color = clrs[[1]],
                 linewidth = .2) +
  geom_linerange(data = data_busco_go_stats |>
                   filter(fst_rank <= 10) |>
                   filter(!duplicated(busco_id)),
                 aes(x = gmid,
                     ymin = -1 -.5*l_hight,
                     ymax = -1 +.5 * l_hight),
                 color = clr_darken(clrs[[2]],.3),
                 linewidth = .2) +
  scale_fill_manual(values = c(`1` = "white",
                               `0` = clr_l_gray),
                    guide = "none") +
  scale_x_continuous(name = "Genomic Position (Gb)", labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(trans = identity,
                                         breaks = arcgaz_genome$mid,
                                         labels = str_remove(arcgaz_genome$chr, "mscaf_a1_")))+
  scale_alpha_manual(values = c(`TRUE` = 6,
                                `FALSE` = 1),
                     guide = "none") +
  labs(y = "BUSCOs") +
  coord_cartesian(xlim = c(0, max(arcgaz_genome$end)),
                  ylim = c(-1.75,1.75),
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
  plot_layout(ncol = 1, heights = c(1,.15,1)) +
  plot_annotation(
    tag_levels = list(c('b', 'c', 'd')),
    theme = theme(plot.tag = element_text(family = fnt_sel),
                  plot.tag.position = c(-.1,.5)))

p3 <- data_win |>
  ggplot(aes(x = gerp_rs_mean,
             y = fst_m_0)) +
  geom_hex(aes(fill = after_stat(log10(count)),
               color = after_scale(clr_darken(fill,.2))),
           linewidth = .2,
           bins = 50) +
  geom_vline(xintercept = median(data_win$gerp_rs_mean),
             color = "black",
             linewidth = .3,
             linetype = 3) +
  geom_hline(yintercept = median(data_win$fst_m_0, na.rm = TRUE),
             color = "black",
             linewidth = .3,
             linetype = 3) +
  scale_fill_gradientn(colours = c(clr_l_gray,
                                   clrs[2],
                                   clrs[1]))+
  labs(x = "Average GERP RS Score per 50kb Window",
       y = "Average *F<sub>ST</sub>* per 50kb Window") +
  guides(fill = guide_colorbar(title.position = "left",
                               barwidth = unit(2.5, "pt"),
                               barheight = unit(50, "pt"))) +
  coord_cartesian(expand = 0,
                  xlim = c(1.1, 1.1) * range(data_win$gerp_rs_mean),
                  ylim = c(-.025,1.025)) +
  theme_ms(fontsize = fs) +
  theme(legend.position = c(0.01,0.01),
        legend.justification = c(0,0),
        legend.title = element_text(angle = 90, hjust = .5),
        axis.title.y = element_markdown())

p4 <- data_win |>
  ggplot(aes(x = gerp_rs_mean)) +
  geom_text(data = tibble(gerp_rs_mean = (c(1.1, 1.1) * range(data_win$gerp_rs_mean) +
                                            median(data_win$gerp_rs_mean) )/2 +
                            c(0,.005),
                          lab = c("less conserved", "highly\nconserved")),
            aes(label = lab, y = 5),
            color = clrs[[1]],
            size = fs / ggplot2::.pt,
            family = fnt_sel) +
  geom_histogram(boundary = 0,
                 color = clrs[[1]],
                 aes(y = after_stat(count *1e-3)),
                 fill = clr_alpha(clrs[[1]], .6),
                 linewidth = plt_lwd,
                 binwidth = .001) +
  geom_vline(xintercept = median(data_win$gerp_rs_mean),
             color = "black",
             linewidth = .3,
             linetype = 3) +
  labs(y = "Count (10<sup>3</sup>)") +
  coord_cartesian(expand = 0,
                  xlim = c(1.1, 1.1) * range(data_win$gerp_rs_mean),
                  clip = "off") +
  theme_ms(fontsize = fs) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_markdown())

p5 <- data_win |>
  ggplot(aes(x = fst_m_0)) +
  geom_text(data = tibble(fst_m_0 = (c(-.025,1.025) +
                                        median(data_win$fst_m_0, na.rm = TRUE) )/2,
                          lab = c("less divergent", "more divergent")),
            aes(label = lab, y = 13),
            color = clr_darken(clrs[[2]], .2),
            size = fs / ggplot2::.pt,
            family = fnt_sel,
            angle = -90) +
  geom_histogram(boundary = 0,
                 color = clrs[[2]],
                 aes(y = after_stat(count *1e-3)),
                 fill = clr_alpha(clrs[[2]], .6),
                 linewidth = plt_lwd,
                 binwidth = .015) +
  geom_vline(xintercept = median(data_win$fst_m_0, na.rm = TRUE),
             color = "black",
             linewidth = .3,
             linetype = 3) +
  labs(y = "Count (10<sup>3</sup>)") +
  coord_flip(expand = 0,
             xlim = c(-.025,1.025),
             clip = "off") +
  theme_ms(fontsize = fs) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_markdown())

pp2 <- p4 + plot_spacer() +
  p3 + p5 +
  plot_layout(widths = c(1,.2), heights = c(.2, 1)) +
  plot_annotation(
    tag_levels = list(c('a', '', '', '')),
    theme = theme(plot.tag = element_text(family = fnt_sel),
                  plot.tag.position = c(-.1,.5)))

pp_out1 <- cowplot::plot_grid(pp2, pp1, rel_widths = c(.66, 1))

ggsave(plot = pp_out1,
       filename = here("results/img/win_gerp_fst.pdf"),
       width = 9, height = 4, device = cairo_pdf)

# ===========================================================
# ===========================================================

pp3 <- (pp2) + (pp1) +
  # plot_layout(widths = c(.5, 1)) +
  plot_annotation(#tag_levels = "a",
    tag_levels = list(c('a', 'b', "", "c")),
    theme = theme(plot.tag = element_text(family = fnt_sel),
                  plot.tag.position = c(-.1,.5))) &
  theme(plot.tag = element_text(family = fnt_sel),
        plot.tag.position = c(-.01, 1))

ggsave(plot = pp2,
       filename = glue("~/Downloads/arcgaz_genome_gerp_{gerp_tag}_fst_{fst_lab}_whg.pdf"),
       width = 10, height = 5, device = cairo_pdf)
# ================
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
