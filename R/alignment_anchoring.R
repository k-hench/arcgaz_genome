library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)

read_size <- \(genome = "", y_base = 0){
  read_tsv(here::here("results", "genome", str_c(genome,".size")),
           col_names = c("chr", "size")) |>
    arrange(-size) |>
    mutate(size_idx = row_number(),
           end = cumsum(size),
           start = lag(end,default =  0),
           eo = size_idx %% 2,
           y_base = y_base,
           genome = genome)
}

genomes <- c("arcgaz_dt_h1_hardmasked",
             "arcgaz_dt_h2_hardmasked",
             "arcgaz_v1_hardmasked",
             "arcgaz_bp_hardmasked",
             "zalcal_v1_hardmasked")

sizes <- tibble(genome = genomes,
                y_base = 0:4) |>
  pmap_dfr(read_size)

genome_summary <- sizes |>
  group_by(genome, y_base) |>
  summarise(n = n(),
            toal_length = max(end))

align_min_length <- 3e5

get_psl <- \(file){
  vroom::vroom(here::here("results","psl", file),
               delim = "\t",
               col_names = c("matches", "misMatches", "repMatches", "nCount",
                             "qNumInsert", "qBaseInsert", "tNumInsert",
                             "tBaseInsert", "strand", "qName", "qSize", "qStart",
                             "qEnd", "tName", "tSize", "tStart", "tEnd",
                             "blockCount")) |>
    select( tName, tStart, tEnd, qName, qStart, qEnd, strand ) |>
    mutate(tSize = tEnd - tStart,
           qSize = qEnd - qStart)
}

simplify_names <- \(str){
  str_remove_all(str, ";.*") |>
    str_replace("ScDpJTZ_", "A.") |>
    str_replace("Sc3htSU_", "HT2.") |>
    str_replace("NC_045", "Z.") |>
    str_replace("NW_0233655", "Z.") |>
    str_replace("ScWAj4l_", "V1.")|>
    str_replace("ptg000", "PB.")
}

n_largest <- 45

inspect_alignment <- \(aln_name = "zalcal_v1", n_alignments = 1e3){
  z_psl <- get_psl(glue("{aln_name}-18.psl.gz"))

  sizes_az <- sizes |>
    filter(size_idx < n_largest,
           genome %in% c("arcgaz_dt_h1_hardmasked", aln_name)) |>
    mutate(l_start = 0,
           chr = simplify_names(chr)) |>
    select(chr, start = l_start, end = size, genome)

  z_psl_above <- z_psl |>
    arrange(-tSize) |>
    # filter(tSize > align_min_length) |>
    mutate(tName = simplify_names(tName),
           qName = simplify_names(qName)) |>
    filter(tName %in% sizes_az$chr,
           qName %in% sizes_az$chr) |>
    head(n_alignments)

  chr_sel <- sizes_az |>
    filter(genome != genomes[[1]]) |>
    pluck("chr")

  chr_sel_a <- sizes_az |>
    filter(genome == genomes[[1]]) |>
    pluck("chr")
  chr_f <- factor(chr_sel, levels = chr_sel)
  n_chr <- length(chr_sel)

  z_psl_in_top_length <- z_psl |>
    mutate(tName = simplify_names(tName),
           qName = simplify_names(qName)) |>
    filter(tName %in% sizes_az$chr,
           qName %in% sizes_az$chr)

  z_psl_matrix <- z_psl_in_top_length |>
    group_by(tName, qName) |>
    summarize(n = n(),
              alignment_length = sum(tSize)) |>
    ungroup() |>
    left_join(sizes_az |> select(qName = chr, qSize = end)) |>
    left_join(sizes_az |> select(tName = chr, tSize = end)) |>
    mutate(tName = factor(tName, chr_sel_a),
           qName = factor(qName, chr_sel),
           t_aling_fract = alignment_length / tSize,
           q_aling_fract = alignment_length / qSize)

  z_psl_matrix <- z_psl_in_top_length |>
    group_by(tName, qName) |>
    summarize(n = n(),
              alignment_length = sum(tSize)) |>
    ungroup() |>
    left_join(sizes_az |> select(qName = chr, qSize = end)) |>
    left_join(sizes_az |> select(tName = chr, tSize = end)) |>
    mutate(tName = factor(tName, chr_sel_a),
           qName = factor(qName, chr_sel),
           t_aling_fract = alignment_length / tSize,
           q_aling_fract = alignment_length / qSize)


  tibble(alignment = aln_name,
         psl = list(z_psl),
         sizes_df = list(sizes_az),
         matrix = list(z_psl_matrix))
}

inspections <- genomes[3] |>
  map_dfr(inspect_alignment)

alignment <- inspections$psl[[1]] |>
  mutate(across(contains("Name"), simplify_names,.names = "{col}_simple"))

matrix_az <- inspections$matrix[[1]]
sizes_az <- inspections$sizes_df[[1]]

best_hits <- matrix_az |>
  group_by(tName) |>
  filter(t_aling_fract == max(t_aling_fract)) |>
  ungroup() |>
  arrange(qName) |>
  mutate(tName = fct_reorder(tName, row_number())) |>
  filter(t_aling_fract > 0.3)
  # filter(t_aling_fract > 0.5)

matrix_az |>
  mutate(tName = factor(tName, levels = levels(best_hits$tName))) |>
  ggplot(aes(x = tName, y = qName)) +
  geom_tile(aes(fill = t_aling_fract,
                width = (.5 + t_aling_fract) / 2,
                height = after_stat(width),
                color = after_scale(clr_darken(fill)))) +
  geom_point(data = best_hits, shape = 1, size = 7) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")

best_hits |>
  ggplot(aes(x = t_aling_fract)) +
  geom_histogram(color = "#F89441FF", fill = clr_lighten("#F89441FF")) +
  theme_minimal()

best_hits

matrix_az |>
  group_by(tName, qName) |>
  filter(t_aling_fract == max(t_aling_fract)) |>
  ungroup()

primary_alignments <- alignment |>
  left_join(best_hits |> dplyr::select(tName_simple = tName, qName_simple = qName) |> mutate(check = TRUE)) |>
  mutate(check = replace_na(check, FALSE)) |>
  filter(check) |>
  group_by(qName) |>
  mutate(tIdx = as.numeric(factor(tName))) |>
  ungroup()

primary_alignments |>
  ggplot(aes(x = abs(tEnd - tStart))) +
  geom_histogram()  +
  facet_wrap(qName_simple ~ .) +
  coord_cartesian(xlim = c(0, 1e4))

primary_borders <- primary_alignments |>
  filter(abs(tEnd - tStart) > 1500) |>
  group_by(tName, qName) |>
  filter(qStart == min(qStart) | qEnd == max(qEnd) | tStart == min(tStart) | tEnd == max(tEnd) ) |>
  ungroup() |>
  mutate(nr = row_number()) |>
  filter(nr %in% c(2, 7, 1, 5, #6, 3,
                   11, 12,
                   33, 35, 37, 38, 39, 40,
                   42, 43, 44, 46,
                   15, 16, 18, 19,
                   145, 148, 146, 149, 147, 150,
                   54, 65, 63, 64,
                   74, 75, 70, 71,
                   81, 80, 77, 76,
                   83, 82,
                   91, 97, 100, 95, 101, 87, 96, 89,
                   156, 152, 154, 155, 157, 159,
                   108, 107, 111, 109,
                   115, 116, 118, 117,
                   121, 119, 123, 120,
                   160, 161, #165, 163,
                   24, 22, 25, 21, 26, 29,
                   126, 125, 127, 128,
                   133, 132, 134, 136))

# primary_borders_o <- primary_borders
# primary_borders <- primary_borders_o

cc <- c( "V1.1", "V1.17", "V1.18", "V1.2", "V1.2441", "V1.25", "V1.28",
         "V1.33", "V1.37", "V1.43", "V1.4319", "V1.44", "V1.46", "V1.49",
         "V1.5180", "V1.8", "V1.82", "V1.84"  )
ccc <- cc[[3]]

p1 <- primary_alignments |>
  filter(abs(tEnd - tStart) > 1500) |>
  # filter(qName_simple == ccc) |>
  ggplot() +
  geom_segment(aes(x = tStart, xend = tEnd,
                   y = qStart, yend = qEnd,
                   color = factor(tIdx))) +
  # geom_point(data = primary_borders,
  #            aes(x = tStart, y = qStart, color = factor(tIdx)),
  #            shape = 1, size = 1) +
  geom_text(data = primary_borders,# |> filter(qName_simple == ccc) ,
             aes(x = tStart, y = qStart, color = factor(tIdx),
                 label = nr), size = 5, fontface ="bold") +
  geom_point(data = primary_borders,# |> filter(qName_simple == ccc) ,
             aes(x = tEnd, y = qEnd, color = factor(tIdx)),
             shape = 1, size = 3) +
  facet_wrap(qName_simple ~ ., scales = 'free') +
  theme(legend.position = "none")

library(plyranges)
p2 <- primary_borders |>
  group_by(tName_simple, qName_simple) |>
  summarise(qstart = min(qStart),
            qend = max(qEnd),
            tstart = min(tStart),
            tend = max(tEnd)) |>
  ungroup() |>
  select(seqnames = qName_simple, start = qstart, end = qend, everything()) |>
  as_granges() |>
  compute_coverage() |>
  as_tibble() |>
  mutate(seqnames = as.character(seqnames)) |>
  ggplot() +
  # geom_area(aes(x = (start + end)/2, y = score))+
  geom_linerange(aes(xmin = start, xmax = end, y = score),
                 linewidth = 2, color = "red") +
  facet_wrap(seqnames ~ ., scales = 'free_x')

p1 + p2

sizes_ht1 <- sizes |> filter(genome == "arcgaz_dt_h1_hardmasked") |> select(tName = chr, tSize = size)
sizes_v1 <- sizes |> filter(genome == "arcgaz_v1_hardmasked") |>  select(qName = chr, qSize = size)

write_lines("ScWAj4l\t1",file = "results/anchoring/weights.txt")

primary_borders |>
  arrange(tName, tStart) |>
  mutate(tPos = floor((tStart + tEnd)/2),
         qPos = floor((qStart + qEnd)/2)) |>
  dplyr::select(tName, qName, tPos, qPos) |>
  left_join(sizes_ht1) |>
  left_join(sizes_v1) |>
  arrange(qName, qPos) |>
  mutate(target_start = tPos - 1,
         querry_label = str_c(str_remove(qName, ";.*") |> str_replace("_","-"),":", sprintf("%.5f", qPos/qSize)),
         # querry_label = str_c(str_remove(qName, ";.*") |> str_replace("_","-"),":",querry),
         target_label = str_c(tName,":", tPos)) |>
  select(tName, target_start, tPos, querry_label, target_label) |>
  write_tsv("results/anchoring/alignment_anchors.bed", col_names = FALSE)

