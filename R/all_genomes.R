library(tidyverse)
library(stringr)

source("R/plot_defaults.R")

read_fai <- \(file){
  genome <- str_remove(file, ".*/") |> str_remove("\\.fa.*")
  if(grepl(pattern = "hap", genome)){
    genome <- str_extract(genome, "hap[1-2]")
    }
  genome
  read_tsv(file,
           col_names = c("name", "length", "offset", "linebases", "linewidth")) |>
    mutate(lc = str_length(name),
           genome = genome,
           o_org = row_number(),
           o_end = cumsum(length),
           o_start = lag(o_end, default = 0),
           o_mid = (o_start + o_end) / 2 ,
           o_type = o_org %% 2) |>
    arrange(-length) |>
    mutate(
      g_ord = row_number(),
      g_end = cumsum(length),
      g_start = lag(g_end, default = 0),
      g_mid = (g_start + g_end) / 2 ,
      g_type = g_ord %% 2,
      in_n50 = g_start < (max(g_end) * .5),
      in_n90 = g_start < (max(g_end) * .9)
    )
}

genome_dir <- "~/work/hoffman_lab/genomes/arcgaz/"
genomes <- c("arcGaz3/arcGaz3.fa.gz.fai",
             "v1/arcgaz_v1.fa.gz.fai",
             "humble_2016_arcgaz_1_2/arcgaz_v1.2.fa.gz.fai",
             "humble_2018_arcgaz_1_4/arcgaz_v1.4.fa.gz.fai",
             "dovetail/Hap1/jordan-zhang-dtg-bie3448-hap1-mb-hirise-ngfa1__09-13-2022__hic_output.fasta.gz.fai",
             "dovetail/Hap2/jordan-zhang-dtg-bie3448-hap2-mb-hirise-x0ohf__08-20-2022__hic_output.fasta.gz.fai",
             "anc/arcgaz_anc_h1.fa.gz.fai", "anc/arcgaz_anc_h2.fa.gz.fai")

genomes_sources <- tibble(genome = c("arcgaz_v1.2", "arcgaz_v1.4", "arcGaz3", "arcgaz_v1",
                                     "hap1", "hap2", "arcgaz_anc_h1", "arcgaz_anc_h2"),
                          source = c("Humble *et al.* 2016", "Humble *et al.* 2018",
                                     "Peart *et al.* 2021", "CeBiTec", "dovetail", "dovetail",
                                     "ALLMAPS", "ALLMAPS"),
                          accession = c("10.5061/dryad.8kn8c", "GCA_900500725.1",
                                        "GCA_900642305.1", "/prj/furseal-genome/Seals/Seal_genome",
                                        NA, NA, NA, NA)) |>
  mutate(genome = factor(genome, levels = genome))

data <- str_c(genome_dir, genomes) |>
  map_dfr(read_fai) |>
  group_by(genome) |>
  mutate(genome = factor(genome, levels = genomes_sources$genome )) |>
  nest() |>
  mutate(n_scaff = map_dbl(data, nrow),
         max_scaff = map_dbl(data, \(df){max(df$length)}),
         total_legth = map_dbl(data, \(df){max(df$g_end)}),
         n50 = map_dbl(data, \(df){sum(df$in_n50)}),
         l50 = map_dbl(data, \(df){min(df$length[df$in_n50])}),
         n90 = map_dbl(data, \(df){sum(df$in_n90)}),
         l90 = map_dbl(data, \(df){min(df$length[df$in_n90])})) |>
  left_join(genomes_sources) |>
  arrange(genome)


data |>
  unnest(cols = data) |>
  ggplot(aes(x = g_ord,
             y = g_end * 1e-9,
             color = genome)) +
  geom_line() +
  facet_grid(.~genome, scales = "free_x") +
  scale_color_manual(values = clrs_n(8))+
  # coord_cartesian(xlim = c(0, 1000)) +
  labs(y = "Cummulative Genome Size (Gb)",
       x = "Scaffold Number") +
  theme_minimal(base_family = fnt_sel) +
  theme(legend.position = "none")

# =====================================
