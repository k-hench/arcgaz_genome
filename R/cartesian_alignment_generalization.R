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
source(here("R/cartesian_alignment_helpers.R"))

clr_genome = "gray50"
manual_order_18 <- list(c( 14, 9, 12, 4, 13, 10, 17, 2, 1, 15, 11, 7, 18, 5, 8, 3, 6, 16),#c(1:18),
                        c( 14, 9, 12, 4, 13, 10, 17, 1, 2, 11, 15, 7, 18, 5, 8, 3, 16, 6),
                      c(14,9,12,4, 17,13,10,2,1,18,11,15,7,5,8,3,16,6),# c(1:18),
                      c(14, 9, 12, 4, 17, 10, 13, 2, 1, 18, 15, 11, 7, 5, 8, 3, 16, 6)) |>
  set_names(nm = c("arcgaz_anc_h2", "arcgaz_anc_h1", "arcgaz_v3_hardmasked", "zalcal_v1"))

genome_width <- .05
skip <- .025

all_alignments <- list(c("arcgaz_anc_h2", "arcgaz_anc_h1"),
     c("arcgaz_anc_h2", "arcgaz_v3_hardmasked"),
     c("arcgaz_anc_h2", "zalcal_v1")) |>
  map_dfr(import_alignment,
          manual_order_18 = manual_order_18,
          genome_width = genome_width,
          skip = skip,
          n_longest = 50)

p1 <- all_alignments |>
  unnest(cols = psl_diag) |>
  ggplot() +
  facet_grid(querry ~ .,switch = "y") +
  geom_diagonal_wide(aes(x = y,
                         y = x,
                         group = str_c(querry, "_" ,group),
                         color = factor(dir),
                         fill = after_scale(clr_lighten(color))),
                     orientation = "y", linewidth = plt_lwd)  +
  geom_linerange(data = all_alignments |> unnest(genome_with_skips),
                 aes(xmin = start, xmax = end, y = y),
                 color = clr_genome,
                 linewidth = .5 * plt_lwd) +
  geom_rect(data = all_alignments |> unnest(sizes_18),
            aes(xmin = start, xmax = end,
                ymin = y_base + (label_sign * skip),
                ymax = y_base + (label_sign * (skip + genome_width))),
            color = clr_genome, fill = clr_lighten(clr_genome,.75),
            linewidth = .5 * plt_lwd) +
  # geom_text(data = all_alignments |> unnest(sizes_18),
  #           aes(x = mid, y = y_base, label = str_c(size_idx,"\n",pre_manual))) +
  # geom_ribbon(data = all_alignments |> unnest(data_cov_windowed),
  #             aes(x = gg_mid, group = seqnames,
  #                 ymin = (genome_idx -1) + 2 * skip *(genome_idx - 1.5) + genome_width *(genome_idx - 2),
  #                 ymax = (genome_idx -1) + 2 * skip *(genome_idx - 1.5) + genome_width *(genome_idx - 2) +  genome_width * avg_cov_25),
  #             linewidth = plt_lwd,
  #             fill = clr_genome) +
  geom_text(data = all_alignments |> unnest(sizes_18),
            aes(x = mid,
                y = y_base + (label_sign *  (4.3 * skip + genome_width)),
                label = chr |> str_remove("mscaf_") |> str_remove("NC_045") |> str_remove("CAAAJK01000")#,
                # hjust = c(1, 0)[as.numeric(factor(label_sign))]
                ),
            color = "black",
            family = fnt_sel,
            size = .35 * fnt_sz * 2,
            angle = 0)  +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.2, 1.25),
                  xlim = c(-1e7, max((all_alignments |> unnest(sizes_18))$end_with_skip)),
                  clip = "off") +
  scale_color_manual(values = clrs |> clr_alpha(.3)) +
  theme_void(base_family = fnt_sel) +
  theme(legend.position = "none",
        strip.text = element_text(angle = 90))

p2 <- all_alignments |> unnest(psl) |>
  ggplot(aes(x = log10(tSize))) +
  facet_grid(querry ~ .,switch = "y") +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2]),
               linewidth = plt_lwd)+
  geom_vline(data = all_alignments |>
               unnest(psl_filtered) |>
               group_by(querry) |>
               filter(tSize == min(tSize)) |>
               ungroup(),
             aes(xintercept = log10(tSize)),
             color = clrs[2], linetype = 3, linewidth = plt_lwd * 4) +
  scale_x_continuous(glue("Alignment Size<br>{round(range((all_alignments |> unnest(psl))$tSize) / c(1, 1e6), digits = 2) %>% str_c(' ',.,c('bp ', 'Mb '), collapse = '–')}, n: {nrow(all_alignments |> unnest(psl))}"),
                     breaks = 2:7, labels = c("100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb")) +
  coord_cartesian(expand = 0) +
  theme_minimal(base_family = fnt_sel) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = clr_genome, linewidth = plt_lwd),
        panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(lineheight = 1.5))

p3 <- all_alignments |>
  unnest(psl_filtered) |>
  ggplot(aes(x = log10(tSize * 1e-5))) +
  facet_grid(querry ~ .,switch = "y") +
  geom_density(color = clrs[2],
               fill = clr_alpha(clrs[2]),
               linewidth = plt_lwd) +
  geom_vline(data = all_alignments |>
               unnest(psl_filtered) |>
               group_by(querry) |>
               filter(tSize == min(tSize)) |>
               ungroup(),
             aes(xintercept = log10(tSize* 1e-5)),
             color = clrs[2], linetype = 3, linewidth = plt_lwd * 4) +
  scale_x_continuous(glue("Alignment Size<br>{sprintf('%.2f',range((all_alignments |> unnest(psl_filtered))$tSize)/1e6) %>% str_c(' ',.,'Mb ', collapse = '–')}"),
                     breaks = c(0,1,2),labels = c("0.1 MB", "1 Mb", "10 Mb")) +
  coord_cartesian(expand = 0) +
  theme_minimal(base_family = fnt_sel) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = clr_genome, linewidth = plt_lwd),
        panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(lineheight = 1.5))

p <- p1 + p2 +# p3+
  plot_layout(nrow = 1, widths = c(1,.2)) +
  plot_annotation(tag_levels = 'A') +
  theme(text = element_text(family = fnt_sel,
                            size = fnt_sz))
