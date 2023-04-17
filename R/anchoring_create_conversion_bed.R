library(tidyverse)
library(prismatic)
library(patchwork)
library(glue)
library(here)
source(here("R/anchoring_assesment.R"))

export_bed <- \(ht){
  write_tsv(
    x = both_beds$bed_export[[ht]],
    file = here(glue("results/anchoring/arcgaz_dt_h{ht}_hardmasked/conversion_h{ht}.bed")),
    col_names = FALSE)
}

export_bed(as.numeric(str_remove(snakemake@params[['ht']],"h")))
