library(tidyverse)

read_fai <- \(file, h = 1){
  read_tsv(file,
           col_names = c("name", "length", "offset", "linebases", "linewidth")) |>
    mutate(hap = h)
}

# scaffold_summary

data <- c(here("data/genomes/arcgaz_anc_h1.fa.gz.fai"),
          here("data/genomes/arcgaz_anc_h2.fa.gz.fai")) |>
  map2_dfr(1:2, read_fai)


data |>
  group_by(hap) |>
  summarise(total_bp = sum(length),
            anchored_bp = sum(length[grepl("mscaf", name)])) |>
  ungroup() |>
  mutate(anchored_frct = sprintf("%.1f", 100 * anchored_bp / total_bp))

# busco summary
read_busco <- \(h){
  read_tsv(here(glue::glue("results/busco/arcgaz_anc_h{h}/run_carnivora_odb10/full_table.tsv")),
           skip = 2) |>
    arrange(`# Busco id`, Score) |>
    group_by(`# Busco id`) |>
    filter(!duplicated(`# Busco id`)) |>
    ungroup() |>
    mutate(hap = h)
  }

data_busco <- 1:2 |>
   map_dfr(read_busco )

data_busco  |>
  group_by(hap) |>
  summarise(n = n(),
            perc_ds = sum(Status == "Complete" | Status == "Duplicated") / length(Status),
            perc_s = sum(Status == "Complete") / length(Status),
            perc_d = sum(Status == "Duplicated") / length(Status),
            perc_f = sum(Status == "Fragmented") / length(Status),
            perc_m = sum(Status == "Missing") / length(Status),
            c_in_18_by_c = sum(Status == "Complete" & grepl("mscaf", Sequence)) /  sum(Status == "Complete"),
            c_in_18_by_total = sum(Status == "Complete" & grepl("mscaf", Sequence)) /  length(Status))


# annotation summary
data_anno <- read_tsv("~/work/hoffman_lab/genomes/arcgaz_anc/annotation/genes.bed.gz",
                      col_names = c("seq", "start", "end", "name"))

# percentage of genome covered by coding region
sprintf("%.2f", 100 * 32136216 / sum(data$length[data$hap == 1]) )

# percentage of gene predictions in to 18 scaffolds
sprintf("%.1f", 100 * sum(grepl("mscaf", data_anno$seq)) / length(data_anno$seq))

