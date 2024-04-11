library(tidyverse)
library(here)
library(knitr)
library(glue)

read_pops <- \(fam){
  read_tsv(here(glue("data/{fam}.pop")), col_names = "Label") |>
    mutate(Family = str_to_sentence(fam))
}

fams <- dir(here('data'), pattern = ".pop") |>
  str_remove(".pop") |>
  map_dfr(read_pops) |>
  filter((!(Family == "Otarioidea")) |
           ( Family == "Otarioidea" & Label == "odoros"))

table_prep <- read_tsv(here("data/pinniped_genome_and_timetree.tsv")) |>
  select(Species = organism_name,
         Label = spec,
         `NCBI accession` = assembly_accession) |>
  mutate(Species = glue::glue("\\textit{{{str_replace_all(str_to_sentence(Species), '_', ' ')}}}")) |>
  left_join(fams) |>
  arrange(Family, Species) |>
  knitr::kable(format = "latex")

table_prep |>
  as.character() |>
  str_remove_all("\\\\hline\\n") |>
  str_replace_all("\\\\textbackslash\\{\\}textit","\\\\textit") |>
  write_lines(file = here("results/tab/aligned_genomes.tex"))
