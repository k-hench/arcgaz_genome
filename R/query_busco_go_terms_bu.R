# input:
# - "results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"
# output:
# - "results/pinniped/go_terms/thresholds.tsv"
library(tidyverse)
library(jsonlite)
library(glue)
library(here)

busco_tab <- read_tsv(here("results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"), skip = 2) |>
  filter(Status == "Complete")

# go_from_busco_id <- \(busco_id, go_type = "bp"){
#   if(!(go_type %in% c("bp", "cc", "mf"))){
#     stop(glue("go_type needs to be one of\n  - 'bp' (biological process)\n  - 'cc' (cellular component)\n  - 'mf' (molecular function)\n<currently it's '{go_type}'>"))
#   }
#   Sys.sleep(.005)
#   type <- c(bp = "biological_process",
#             cc = "cellular_component",
#             mf = "molecular_function")[go_type]
#
#   raw_data <- fromJSON(glue("https://v10.orthodb.org/group?id={busco_id}"))$data
#   raw_data[[type]]$id
#   }

z <- 0
n_buscos <- length(busco_tab$`# Busco id`)

l <- list()
go_from_busco_id <- \(busco_id){
  Sys.sleep(.025)
  z <<- z +1
  cat("\r")
  cat(glue("--- ({str_pad(round(z/n_buscos * 100), width = 3, pad = 0)}%) ---"))
  raw_data <- fromJSON(glue("https://v10.orthodb.org/group?id={busco_id}"))$data
  l_out <- list(bp = raw_data[["biological_process"]]$id,
                cc = raw_data[["cellular_component"]]$id,
                mf = raw_data[["molecular_function"]]$id)
  l <<- c(l, l_out)
  l_out
}

go_terms <- busco_tab$`# Busco id` |>
  map(go_from_busco_id) |>
  set_names(nm = busco_tab$`# Busco id`)

# saveRDS(go_terms, "~/Downloads/go_terms.Rds")


# go_terms <- readRDS("~/Downloads/go_terms.Rds")
# pick_sub <- \(L, sub = "bp"){L[[sub]]}
#
# go_bp_prep <- go_terms |>
#   map(pick_sub) |>
#   set_names(nm = names(go_terms))
#
# go_bp_prep |>
#   map2_dfr(names(go_bp_prep),
#            \(x, y){tibble(busco_id = y,
#                           go_terms = str_c(x, collapse = ",") |>
#                             str_replace(",", ", "))}) |>
#   filter(go_terms != "") |>
#   write_tsv("results/pinniped/busco_id_to_go_term_bp.tsv", col_names = FALSE)
#
