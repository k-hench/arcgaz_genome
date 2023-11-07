library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
tsv_in <- args[[1]]
vcf_out <- args[[2]]

data <- read_tsv(tsv_in)

ref <- names(data)[[3]]
select_alleles <- \(x, a){replace_na(as.character(a[x]), ".")}

data |> 
  mutate(rn = row_number()) |> 
  rowwise() |> 
  nest(data = -c(refSequence, refPosition, rn)) |> 
  mutate(REF = map_chr(data, \(df){df[[ref]]}),
         ALT_LIST = map2(data, REF,
                             \(df, r){
                               a <- unique(unlist(df[1,]))
                               sort(a[a != r])}),
         ALT = map_chr(ALT_LIST, \(a){str_c(a, collapse = ",")}),
         ALLELE_NR = map2(REF, ALT_LIST, \(r,a){
           nm = c(r, a)
           (seq_along(nm) -1) |> set_names(nm = nm)
           })) |> 
  unnest(data) |> 
  ungroup() |>
  mutate(across(-c(refSequence, refPosition, rn, REF, ALT_LIST, ALT, ALLELE_NR),
               ~ map2_chr(.x, .y = ALLELE_NR, select_alleles))) |> 
  mutate(ID = ".",
         QUAL = ".",
         FILTER = ".",
         INFO = ".",
         FORMAT = "GT") |> 
  select(-c(rn, ALT_LIST, ALLELE_NR)) |> 
  select(`#CHROM` = refSequence, POS = refPosition, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, everything()) |> 
  write_tsv(file = vcf_out)
