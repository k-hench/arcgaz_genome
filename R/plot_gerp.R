library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(ggtext)
source("R/plot_defaults.R")


fs <- 9
clr_l_gray <- rgb(.85,.85,.85)
scfs <- str_c("mscaf_a1_",
              c(str_pad(1:17, width = 2, pad =  0), "x"))

read_win <- \(scf){read_tsv(glue("results/pinniped/gerp/tsv/gerp_snps_{scf}_summary.tsv.gz"))}
read_busco <- \(scf){read_tsv(glue("results/pinniped/gerp/tsv/gerp_busco_{scf}_summary.tsv.gz"))}

genome <- read_tsv("data/genomes/arcgaz_anc_h1.genome",
                   col_names = c("chr", "length")) |>
  mutate(end = cumsum(length),
         start = lag(end, default = 0),
         mid = (start + end)/2,
         eo = row_number() %% 2)

busco_meta <- read_tsv("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv", skip = 2) |>
  filter(Status == "Complete") |>
  select(buco_id = `# Busco id`,
         bstart = `Gene Start`,
         bend = `Gene End`)

data_win <- scfs |>
  map_dfr(read_win) |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = start + scaf_start,
         gend = end + scaf_start,
         gmid = (gstart + gend)/2)

data_buso <- scfs |>
  map_dfr(read_busco) |>
  left_join(busco_meta) |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = bstart + scaf_start,
         gend = bend + scaf_start,
         gmid = (gstart + gend)/2)

p1 <- ggplot() +
  geom_rect(data = genome,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = factor(eo))) +
  geom_point(data = data_win,
             aes(x = gmid, y = gerp_rs_mean,
                 color = factor(eo),
                 group = chr),
             size = .3) +
  geom_smooth(data = data_win,
              color = "white",
              linewidth = .5,
              span = 0.5,
             aes(x = gmid, y = gerp_rs_mean,
                 group = chr#,
                 # fill = factor(eo),
                 # color = after_scale(fill)
                 ),
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
  labs(y = "average RS score (GERP)") +
  coord_cartesian(xlim = c(0, max(genome$end)),
                  expand = 0)+
  theme_ms(fontsize = fs)

data_buso |>
  ggplot(aes(x = gerp_rs_mean)) +
  geom_vline(xintercept = mean(data_buso$gerp_rs_mean),
             color = clrs[[1]], linetype = 3)+
  geom_histogram(boundary = 0,
                 color = clrs[[2]],
                 fill = clr_alpha(clrs[[2]], .6),
                 size = plt_lwd,
                 binwidth = .001) +
  coord_cartesian(expand = 0)+
  theme_ms(fontsize = fs)

data_fst_win <- read_tsv("results/fst_win.tsv.gz") |>
  left_join(genome |> select(CHROM = chr, start = start, eo)) |>
  mutate(gstart = BIN_START + start,
         gend = BIN_END + start,
         gmid = (gstart + gend)/2,
         fst_w_0 = if_else(WEIGHTED_FST > 0, WEIGHTED_FST, 0),
         fst_m_0 = if_else(MEAN_FST > 0, MEAN_FST, 0))

fst_w_breaks <- quantile(data_fst_win$fst_w_0,
                         probs = c(.005, .995))

p2 <- ggplot() +
  geom_rect(data = genome,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = factor(eo))) +
  geom_point(data = data_fst_win,
             aes(x = gmid, y = fst_w_0,
                 color = factor(eo),
                 group = CHROM),
             size = .3) +
  geom_smooth(data = data_fst_win,
              color = "white",
              linewidth = .5,
              span = 0.5,
              aes(x = gmid,
                  y = fst_w_0,
                  group = CHROM#,
                  # fill = factor(eo),
                  # color = after_scale(fill)
              ),
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
  labs(y = "*F<sub>ST</sub>* (weighted mean)") +
  coord_cartesian(xlim = c(0, max(genome$end)),
                  expand = 0)+
  theme_ms(fontsize = fs)+
  theme(axis.title.y = element_markdown())


plain_busco <- read_tsv("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv", skip = 2) |>
  filter(Status == "Complete") |>
  mutate(start = if_else(`Gene Start` < `Gene End`, `Gene Start`, `Gene End`),
         end = if_else(`Gene Start` < `Gene End`, `Gene End`, `Gene Start`)) |>
  select(chr = Sequence,
         start, end,
         busco_id = `# Busco id`)  |>
  left_join(genome |> select(chr, scaf_start = start, eo)) |>
  mutate(gstart = start + scaf_start,
         gend = end + scaf_start,
         gmid = (gstart + gend)/2)

pb <- ggplot() +
  geom_rect(data = genome,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = factor(eo))) +
  geom_point(data = plain_busco,
             aes(x = gmid, y = 0, color = factor(eo)),
             shape = "|", alpha = .3) +
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
  theme_void(base_family = fnt_sel, base_size = fs)#+
  # theme(axis.title.y = element_text(angle = 90))

pp1 <- p1 +
  theme(axis.title.x.bottom = element_blank(),
           axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()) +
  pb +
  p2 +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) +
  plot_layout(ncol = 1, heights = c(1,.1,1))

data_win_combined <- data_fst_win |>
  select(chr = CHROM, end = BIN_END, fst_w_0) |>
  left_join(data_win |> select(chr, start, end, gmid, gerp_rs_mean))

p3 <- data_win_combined |>
  ggplot(aes(x = gerp_rs_mean, y = fst_w_0)) +
  geom_hex(aes(fill = after_stat(log10(count)),
               color = after_scale(clr_darken(fill,.2))),
           linewidth = .2,
           bins = 50) +
  scale_fill_gradientn(colours = c(clr_l_gray,
                                   clrs[2],
                                   clrs[1]))+
  labs(x = "average RS score (GERP)",
       y = "*F<sub>ST</sub>* (weighted mean)") +
  guides(fill = guide_colorbar(title.position = "left",
                               barwidth = unit(5, "pt"),
                               barheight = unit(100, "pt"))) +
  theme_ms(fontsize = fs) +
  theme(legend.position = c(0.01,0.01),
        legend.justification = c(0,0),
        legend.title = element_text(angle = 90, hjust = .5),
        axis.title.y = element_markdown())

p3 + pp1 +
  plot_layout(widths = c(.5, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = fnt_sel)))

