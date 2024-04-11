library(tidyverse)
library(here)
library(glue)
library(jsonlite)
library(ggtext)

source(here("R/plot_defaults.R"))

seq_levels <- c("scaffold",
                "contig",
                "gap")

read_assembly_stats <- \(gn){
  file <- here("results", "fa_stats", str_c(gn, "_stats_a.json"))

  js <- read_json(file)

  bind_rows(as_tibble(js$`Contig Stats`) |> mutate(seq_type = "contig"),
            as_tibble(js$`Scaffold Stats`) |> mutate(seq_type = "scaffold")) |>
    mutate(genome = gn) |>
    set_names(nm = str_to_lower) |>
    mutate(seq_type = factor(seq_type, levels = seq_levels))
}

read_busco_summary <- \(gn){
  file <- here("results", "busco", gn, glue("short_summary.specific.carnivora_odb10.{gn}.txt"))

  read_tsv(file,
           skip = 9, n_max = 5, col_names = c("n1", "count", "stat","n2","n3","n4")) |>
    dplyr::select(count:stat) |>
    mutate(stat = str_replace(stat, " \\(", "_") |> str_remove("\\)")) |>
    separate(stat, into = c("label", "abbrev"), sep = '_') |>
    mutate(abbrev = factor(abbrev, levels = c("C","S", "D", "F", "M"))) |>
    arrange(-as.numeric(abbrev)) |>
    mutate(genome = gn,
           prct = count / 14502,
           lab_y = c(1:4,0) * .2,#lag(cumsum(prct), default = 0) + .5 * prct,
           x_lab = c(1.6,1.7,1.8,1,.5))
}

genomes <- c("arcgaz", "calurs", "eumjub", "halgry", "lepwed",
             "mirang", "mirleo", "neosch", "odoros", "phovit",
             "zalcal")

data_n50 <- genomes |> map_dfr(read_assembly_stats)
data_busco <- genomes |> map_dfr(read_busco_summary)

fs <- fnt_sz

p <- data_n50 |>
  left_join(data_busco |>  filter(abbrev == "C")) |>
  filter(seq_type == "scaffold") |>
  ggplot(aes(x = log10(n50),
             y = prct)) +
  geom_point(aes(color = genome == "arcgaz"),
             show.legend = FALSE,
             alpha = .6,
             stroke = .4,
             size = .5) +
  geom_text(aes(label = genome,
                color = genome == "arcgaz",
                x = log10(n50) + c(1,0,1,1,1,  1,1,1,1,1,  1) * .16,
                y = prct +  c(.3,.6,0,0,0,  0,0,-.3,0,.6,  0) * .01),
            family = fnt_sel,
            size = fs / ggplot2::.pt,
            show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = clrs[[2]], `FALSE` = "black")) +
  scale_x_continuous(name = "N50",
                     breaks = 6:8,
                     #   labels = \(x){glue("10<sup>{x}</sup>")})+
                     labels = \(x){sprintf("%.0fMb", 10 ^ (x - 6) )}) +
  scale_y_continuous("BUSCO score (%)", labels = \(x){sprintf("%.0f", x * 100)}) +
  theme_ms(fontsize = fs)

ggsave(filename = here("results/img/pinni_busco_n50.pdf"),
       plot = p,
       width = fhalfwidth,
       height = fhalfwidth * .75,
       device = cairo_pdf)
