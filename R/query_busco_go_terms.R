# input:
# - "results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"
# output:
# - expand( "results/pinniped/go_terms/busco_id_to_go_term_{type}.tsv", type = [ "bp", "cc", "mf" ] )
library(tidyverse)
library(jsonlite)
library(glue)
library(here)

busco_tab <- read_tsv(here("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"), skip = 2) |>
  filter(Status == "Complete")

z <- 0
n_buscos <- length(busco_tab$`# Busco id`)

l <- list()
go_from_busco_id <- \(busco_id){
  Sys.sleep(.025)
  z <<- z +1
  cat("\r")
  cat(glue("--- ({str_pad(round(z/n_buscos * 100), width = 3, pad = 0)}%) ---"))
  raw_data <- fromJSON(glue("https://v10-1.orthodb.org/group?id={busco_id}"))$data
  l_out <- list(bp = raw_data[["biological_process"]]$id,
                cc = raw_data[["cellular_component"]]$id,
                mf = raw_data[["molecular_function"]]$id)
  l <<- c(l, l_out)
  l_out
}

go_terms <- busco_tab$`# Busco id` |>
  map(go_from_busco_id) |>
  set_names(nm = busco_tab$`# Busco id`)

write_lines(x = glue("Download report for BUSCO related GO terms\nAccessed at {Sys.time()} from https://www.orthodb.org/ (Version 10.1) \nUsing the API in the form of https://v10-1.orthodb.org/group?id=<busco_id>"),
            file = here(glue("results/pinniped/go_terms/busco_id_to_go_term_info.txt")))

pick_sub <- \(L, sub = "bp"){L[[sub]]}

go_export <- \(type){
  go_prep <- go_terms |>
    map(pick_sub, sub = type) |>
    set_names(nm = names(go_terms))

  go_prep |>
    map2_dfr(names(go_prep),
             \(x, y){tibble(busco_id = y,
                            go_terms = str_c(x, collapse = ",") |>
                              str_replace(",", ", "))}) |>
    filter(go_terms != "") |>
    write_tsv(here(glue("results/pinniped/go_terms/busco_id_to_go_term_{type}.tsv")),
              col_names = FALSE)
}

c("bp", "cc", "mf") |>
  walk(go_export)
