library(tidyverse)
library(jsonlite)
library(glue)

busco_tab <- read_tsv("results/busco/zalcal_v1/run_carnivora_odb10/full_table.tsv", skip = 2) |>
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
head(go_terms)

reduce_go_list <- \(by){
  go_terms |> reduce(.f = \(l1,l2){c(l1, l2[by])}, .init = list()) |> set_names(nm = names(go_terms))
}

go_terms_bp <- reduce_go_list("bp")
go_terms_cc <- reduce_go_list("cc")
go_terms_mf <- reduce_go_list("mf")

go_terms_bp |>
  head() |>
  map2_dfr(.y = names(go_terms_bp |>
                        head()),
           \(l,n){tibble(busco_id = n, go_terms = str_remove(str_c(" ", l, collapse = ","), "^ "))})

str(go_terms_bp |>
      head())
