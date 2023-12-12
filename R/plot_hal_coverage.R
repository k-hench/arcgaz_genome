library(tidyverse)
library(here)
library(glue)
library(prismatic)
library(ggtext)
source(here("R/plot_defaults.R"))

read_cov <- \(scf){
  read_tsv(here(glue("results/neutral_tree/cov/{scf}.collapsed.bed.gz")),
           col_names = c("chr", "start", "end", "cov"))
}

read_cov_fam <- \(scf, fam){
  read_tsv(here(glue("results/neutral_tree/cov/fam/{fam}-{scf}.collapsed.bed.gz")),
           col_names = c("chr", "start", "end", "cov"))
}

scfs <- str_c("mscaf_a1_", c(str_pad(1:17, width = 2, pad = 0), "x"))

data <- scfs |> map_dfr(read_cov)
data_ota <- scfs |> map_dfr(read_cov_fam, fam = "ota")
data_pho <- scfs |> map_dfr(read_cov_fam, fam = "pho")

cov_summary <- \(dat, grp){
  dat |>
    group_by(coverage = cov + 1) |>
    summarise(bp = sum(end - start)) |>
    ungroup() |>
    mutate(z_bp = cumsum(bp),
           a_bp = lag(z_bp, default = 0),
           percent = sprintf("%.1f", bp / sum(bp) * 100),
           group = grp)
}

data_all_summary <- data |> cov_summary(grp = "all")
data_ota_summary <- data_ota |> cov_summary(grp = "Otariidae")
data_pho_summary <- data_pho |> cov_summary(grp = "Phocidae")

data_summary <- data_all_summary |>
  bind_rows(data_ota_summary) |>
  bind_rows(data_pho_summary)

clr_grp <- c(all = "black", Otariidae = clrs[[1]], Phocidae = clrs[[2]])

clr_label <- \(grp){glue("<span style='color:{clr_grp[grp]}'>{grp}</span>")}

p <- data_summary |>
  pivot_longer(c(z_bp,a_bp),
               values_to = "position") |>
  arrange(group, -coverage, -position) |>
  ggplot(aes(x = max(data_summary$z_bp) - position, y = coverage)) +
  geom_area(aes(color = group,
                fill = after_scale(clr_alpha(color))),
            linewidth = plt_lwd) +
  geom_text(data = data_summary,
            mapping = aes(x = max(z_bp) - a_bp + 8e7,
                          y = coverage - .5,
                          label = coverage,
                          color = group),
            family = fnt_sel, size = fnt_sz / ggplot2::.pt,
            fontface = "bold") +
  scale_x_continuous("Genome Sequence (Gb)",
                     labels = \(x){sprintf("%.1f", x*1e-9)},
                     sec.axis = sec_axis("Genome Share (%)",
                                         trans = identity,
                                         breaks = seq(0,max(data_summary$z_bp),
                                                      length.out = 11),
                                         labels = 0:10 * 10)) +
  scale_y_continuous("Alignment Depth (n genomes)",
                     breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_color_manual(values = clr_grp, guide = "none") +
  ggforce::facet_col(factor(clr_label(group), levels = clr_label(names(clr_grp))) ~ .,
                     scales = "free_y",
                     space = "free",
                     strip.position = "right") +
  coord_cartesian(#ylim = c(0, 11.25),
                  xlim = c(0, 1.04 * max(data_summary$z_bp)),
                  expand = 0) +
  theme_ms() +
  theme(strip.text.y.right = element_markdown(size = fnt_sz),
        panel.spacing.y = unit(5, "pt"))

ggsave(filename = here("results/img/hal_coverage.pdf"),
       plot = p,
       width = 4, height = 4, device = cairo_pdf)

data_summary |>
  select(coverage, percent, group) |>
  pivot_wider(names_from = group, values_from = percent)
