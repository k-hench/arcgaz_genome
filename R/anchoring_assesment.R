# library(tidyverse)
# library(prismatic)
# library(patchwork)
# library(glue)
# library(here)

source(here("R/plot_defaults.R"))

new_prefix <- "mscaf_a" # beware of alphabetical sorting (< "S")
read_fai <- \(file){
  read_tsv(file,
           col_names = c("name", "length", "offset", "linebases", "linewidth"))
}

compile_bed <- \(ht = 1){
  scaf_x = str_c(new_prefix, ht,"_x")
  full_bed <- read_tsv(here(glue("results/anchoring/arcgaz_dt_h{ht}_hardmasked/anchored_arcgaz_dt_h{ht}.lifted.bed")),
                       col_names = c("name", "start", "end", "n1", "n2"))
  data <- full_bed |>
    mutate(org_scaffold_v3 = str_remove(n1, ":.*"),
           org_scaffold_dt = str_remove(n2, ":.*"),
           is_x = org_scaffold_v3 == "CAAAJK-010000043.1",
           is_single_diget = str_detect(name, "chr[0-9]$"),
           name_prep = case_when(
             is_x ~ scaf_x,
             !is_x & str_sub(name, 1, 3) == "chr" ~ str_c("chr", str_pad(str_sub(name, 4, -1), pad = 0, width = 4)),
             .default = name)) |>
    arrange(name_prep, start) |>
    mutate( new_chr = !(name_prep == lag(name_prep, default = "")),
            chr_nr = cumsum(new_chr),
            name_new = case_when(
              is_x ~ scaf_x,
              .default = str_c(new_prefix, ht,"_",str_pad(chr_nr, width = 2, pad = 0)))) |>
    select(name, name_new, org_scaffold_v3) |>
    group_by(name_new) |>
    mutate(n_scaf = length(name_new) / 2) |>
    ungroup() |>
    filter(!duplicated(name))

  faidx <- read_fai(here(glue("results/anchoring/arcgaz_dt_h{ht}_hardmasked/anchored_arcgaz_dt_h{ht}.fasta.gz.fai")))

  data_bed <- faidx |>
    left_join(data) |>
    mutate(name_new = if_else(is.na(name_new), name, name_new),
           start = 0) |>
    arrange(str_to_lower(name_new)) |>
    select(name, start, end = length, name_new, org_scaffold_v3, n_scaf) |>
    mutate(gend = cumsum(end),
           gstart = lag(gend, default = 0),
           gmid= (gstart + gend) / 2,
           lg_group = case_when(
             name_new == glue("{new_prefix}{ht}_x") ~ "x",
             row_number() > 18 ~ "unplaced",
             .default = c("even", "odd")[1 + (row_number() %% 2)] ) |>
             factor(levels = c("odd", "even", "x", "unplaced")))

  bed <- data_bed |>
    select(name, start, end, name_new)

  data_scaf <- full_bed |>
    mutate(org_scaf = str_remove(n2, ":.*"),
           scaf_abbref = str_remove(org_scaf, "^.*_") |> str_remove(";.*")) |>
    group_by(org_scaf) |>
    summarise(name = name[1],
              scaf_abbref = scaf_abbref[1],
              start = min(start) + 1,
              end = max(end)) |>
    ungroup() |>
    left_join(data_bed |>
                select(name, gstart, gmid)) |>
    mutate(gstart_s = gstart + start,
           gend_s = gstart + end,
           gmid_s = (gstart_s + gend_s) / 2) |>
    arrange(gstart_s) |>
    mutate(scaf_type = row_number() %% 2) |>
    select(name, org_scaf, scaf_abbref, gstart:scaf_type)

  p <- data_bed |>
    filter(grepl( new_prefix, name_new )) |>
    ggplot() +
    geom_rect(data = tibble(end = max(data_bed$gend)),
              aes(xmin = 0, xmax = end * 1e-9,
                  ymin = 0, ymax = 1),
              color = "gray80", fill = "gray95",
              linewidth = plt_lwd) +
    geom_rect(aes(xmin = gstart * 1e-9,
                  xmax = gend * 1e-9,
                  ymin = 0, ymax = 1,
                  color = lg_group,
                  fill = after_scale(clr_lighten(color))),
              linewidth = plt_lwd) +
    geom_rect(data = data_scaf,
              aes(xmin = gstart_s * 1e-9,
                  xmax = gend_s * 1e-9,
                  ymin = .25 + scaf_type * .25,
                  ymax = .5 + scaf_type * .25,
                  fill = factor(scaf_type),
                  color = after_scale(clr_darken(fill))),
              alpha = .3,
              linewidth = plt_lwd) +
    geom_point(aes(x = gmid * 1e-9,
                   y = .5),
               shape = 21, size = 4.5, stroke = plt_lwd,
               color = "black", fill = "white") +
    geom_text(data = data_bed |>
                select(name, name_new, gmid, gstart, gend, lg_group, n_scaf) |>
                mutate(placement_type = lg_group == "unplaced",
                       plot_group = if_else(!placement_type, name_new, "unplaced"),
                       n_scaf = replace_na(n_scaf, 1)) |>
                group_by(plot_group) |>
                summarize(n_scaf = sum(n_scaf),
                          gmid = (min(gstart) + max(gend))/2),
              aes(x = gmid * 1e-9,
                   y = .5,
                  label = n_scaf),
               color = "black",
              family = fnt_sel, size = fnt_sz/ggplot2::.pt) +
    geom_text(data = data_scaf,
              aes(x = gmid_s * 1e-9,
                  y = .2 + scaf_type * .6,
                  label = scaf_abbref),
              family = fnt_sel, size = fnt_sz/ggplot2::.pt) +
    scale_color_manual(values = c(clr_darken(clrs[1], .25), clrs[1], "gray25", "gray90")) +
    scale_fill_manual(values = c("gray65", "black")) +
    scale_x_continuous("genomic position (GB)",
                       sec.axis = sec_axis(breaks = c(data_bed$gmid[1:18],
                                                      (data_bed$gstart[19] + max(data_bed$gend)) / 2 )* 1e-9,
                                           labels = c(str_c(
                                             data_bed$name_new[1:18] |> str_remove("mscaf_"),
                                                            "\n", str_remove(data$org_scaffold_v3[1:18], "CAAAJK-01000") |> str_remove("^0*") %>% str_c("v3_",.)),
                                                      " unplaced"),
                                           trans = identity)) +
    scale_y_continuous(name = glue("haplotype {ht}"),
                       limits = c(-.3, 1.1)) +
    coord_cartesian(ylim = c(-.01, 1.01),
                    expand = 0) +
    # facet_wrap(name ~ . ) +
    theme_minimal(base_family = fnt_sel, base_size = fnt_sz) +
    theme(axis.text.y = element_blank())

  tibble(ht = ht,
         plot = list(p),
         bed_export = list(bed),
         data_bed = list(data_bed),
         data_scaf = list(data_scaf))
}

both_beds <- 1:2 |>
  map_dfr(compile_bed)
# both_beds$plot |>
#   wrap_plots(ncol = 1) &
#   theme(legend.position = "none")

