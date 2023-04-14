library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(plyranges)
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
  summarise(n = dplyr::n(),
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

inspect_alignment <- \(aln_name = "arcgaz_dt_h2_hardmasked",
                       n_alignments = 1e3,
                       ref = "arcgaz_v1_hardmasked"){
  z_psl <- get_psl(glue("slim_{aln_name}_on_arcgaz_v1_hardmasked.psl.gz"))

  sizes_az <- sizes |>
    mutate(in_larges_scaffold = if_else(genome == ref, 25, n_largest)) |>
    filter(size_idx < in_larges_scaffold,
           genome %in% c("arcgaz_v1_hardmasked", aln_name)) |>
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
    filter(genome != genomes[[3]]) |>
    pluck("chr")

  chr_sel_a <- sizes_az |>
    filter(genome == genomes[[3]]) |>
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
    summarize(n = dplyr::n(),
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
    summarize(n = dplyr::n(),
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

get_best_hits <- \(alignment, psl, matrix, sizes_df){

  alignment_df <- psl |>
    mutate(across(contains("Name"), simplify_names,.names = "{col}_simple"))

  sizes_az <- sizes_df
  matrix_az <- matrix

  best_hits <- matrix_az |>
    group_by(qName) |>
    filter(q_aling_fract == max(q_aling_fract)) |>
    ungroup() |>
    arrange(tName) |>
    mutate(qName = fct_reorder(qName, row_number())) |>
    filter(q_aling_fract > 0.33)

  p_fract <- matrix_az |>
    mutate(qName = factor(qName, levels = levels(best_hits$qName))) |>
    ggplot(aes(x = as.numeric(tName), y = qName)) +
    geom_tile(aes(fill = q_aling_fract,
                  width = (.5 + q_aling_fract) / 2,
                  height = after_stat(width),
                  color = after_scale(clr_darken(fill))),
              alpha = .5) +
    geom_point(data = best_hits, shape = 1, size = 7, alpha = .5) +
    scale_fill_viridis_c(option = "C", direction = -1) +
    scale_x_continuous(breaks = 1:24,
                       labels = levels(matrix_az$tName),
                       sec.axis = sec_axis(breaks = (1:12) * 2,
                                           trans = ~ .)) +
    coord_equal(xlim = c(0.5, 24.5),
                ylim = c(0.5, n_largest),
                expand = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90),
          axis.text.x.top = element_text(angle = 0),
          legend.position = "bottom")

  p_bar <- best_hits |>
    ggplot(aes(x = t_aling_fract)) +
    geom_histogram(color = "#F89441FF", fill = clr_lighten("#F89441FF")) +
    theme_minimal()

  # best_hits
  # matrix_az |>
  #   group_by(tName, qName) |>
  #   filter(t_aling_fract == max(t_aling_fract)) |>
  #   ungroup()

  primary_alignments <- alignment_df |>
    left_join(best_hits |> dplyr::select(tName_simple = tName, qName_simple = qName) |> mutate(check = TRUE)) |>
    mutate(check = replace_na(check, FALSE)) |>
    filter(check) |>
    group_by(tName) |>
    mutate(qIdx = as.numeric(factor(qName))) |>
    ungroup()

  p_hist <- primary_alignments |>
    ggplot(aes(x = abs(tEnd - tStart))) +
    geom_histogram()  +
    facet_wrap(qName_simple ~ .) +
    coord_cartesian(xlim = c(0, 1e4))

  tibble(alignment = alignment,
         best_hits = list(best_hits),
         primary_alignments = list(primary_alignments),
         p_frct = list(p_fract),
         p_bar = list(p_bar),
         p_hist = list(p_hist))
}

inspections <- genomes[1:2] |>
  map_dfr(inspect_alignment)

best_hits_df <- inspections |>
  pmap_dfr(get_best_hits)

# wrap_plots(best_hits_df$p_frct)
# wrap_plots(best_hits_df$p_hist)

plot_primary_borders <- \(alignment, primary_alignments, ...){
  size_cuttoff <- 3500
  primary_borders <- primary_alignments |>
    filter(abs(tEnd - tStart) > size_cuttoff) |>
    group_by(tName, qName) |>
    filter(qStart == min(qStart) | qEnd == max(qEnd) | tStart == min(tStart) | tEnd == max(tEnd) ) |>
    ungroup() |>
    mutate(nr = row_number())

  # primary_borders_o <- primary_borders
  # primary_borders <- primary_borders_o

  cc <- c( "V1.1", "V1.17", "V1.18", "V1.2", "V1.2441", "V1.25", "V1.28",
           "V1.33", "V1.37", "V1.43", "V1.4319", "V1.44", "V1.46", "V1.49",
           "V1.5180", "V1.8", "V1.82", "V1.84"  )
  # ccc <- cc[[18]]

  border_ids <- tibble(
    genome = genomes[1:2],
    ids = list(
      c(23, 22, 114, 115, 70, 71, 100, 99, # 01 ht 1
        80, 81, 136, 137, 28, 29,          # 02
        69, 68, 5, 6,                      # 03
        66, 65, 63, 64,                    # 04
        46, 47, 150, 151, 111, 112,        # 05
        32, 33, 62, 61,                    # 06
        13, 14, 107, 105,                  # 07
        1, 2, 56, 57,                      # 08
        9, 10,                             # 09
        54, 52, 133, 134, 73, 74,          # 10
        146, 147, 7, 8, 126, 125,          # 11
        24, 25, 144, 145,                  # 12
        37, 38, 77, 76,                    # 13
        139, 138, 40, 39,                  # 14
        27, 26,                            # 15
        21, 20, 148, 149, 104, 102,        # 16
        108, 109, 79, 78,                  # 17
        142, 141, 43, 42),                 # 18
      c(31, 32, 108, 107,                  # 01  ht2
        13, 14, 75, 76,                    # 02
        40, 41,                            # 03
        86, 85, 43, 42,                    # 04
        123, 120, 129, 130, 104, 105,      # 05
        5, 6, 57, 56,                      # 06
        94, 95, 36, 35,                    # 07
        49, 48, 137, 136, 90, 93, 15, 17,  # 08
        23, 25, 45, 44, 8, 7,              # 09
        60, 58, 53, 52, 99, 100, 88, 89,   # 10
        81, 82, 96, 97, 71, 63,            # 11
        46, 47, 133, 134,                  # 12
        33, 34, 51, 50,                    # 13
        110, 109, 28, 27,                  # 14
        2, 1,                              # 15
        30, 29, 138, 139, 111, 112,        # 16
        37, 38, 12, 11,                    # 17
        102, 103, 79, 77, 20, 19)          # 18
    ))

  p1 <- primary_alignments |>
    filter(abs(qEnd - qStart) > size_cuttoff/3) |>
    # filter(tName_simple == ccc) |>
    ggplot() +
    geom_segment(aes(x = tStart, xend = tEnd,
                     y = qStart, yend = qEnd,
                     color = factor(qIdx))) +
    geom_text(data = primary_borders |>
                filter(#tName_simple == ccc,
                  nr %in% (border_ids |> filter(genome == alignment) |> pluck("ids") |> unlist())
                ) ,
              aes(x = tStart, y = qStart, color = factor(qIdx),
                  label = nr), size = 5, fontface ="bold") +
    geom_point(data = primary_borders,# |> filter(tName_simple == ccc) ,
               aes(x = tEnd, y = qEnd, color = factor(qIdx)),
               shape = 1, size = 3) +
    facet_wrap(tName_simple ~ ., scales = 'free') +
    theme(legend.position = "none")

  p2 <- primary_borders |>
    filter(nr %in% (border_ids |> filter(genome == alignment) |> pluck("ids") |> unlist())) |>
    group_by(qName_simple, tName_simple) |>
    summarise(tstart = min(tStart),
              tend = max(tEnd),
              qstart = min(qStart),
              qend = max(qEnd)) |>
    ungroup() |>
    select(seqnames = tName_simple, start = tstart, end = tend, everything()) |>
    as_granges() |>
    compute_coverage() |>
    as_tibble() |>
    mutate(seqnames = as.character(seqnames)) |>
    ggplot() +
    # geom_area(aes(x = (start + end)/2, y = score))+
    geom_linerange(aes(xmin = start, xmax = end, y = score),
                   linewidth = 2, color = "red") +
    facet_wrap(seqnames ~ ., scales = 'free_x')

  sizes_t <- sizes |> filter(genome == "arcgaz_v1_hardmasked") |> select(tName = chr, tSize = size)
  sizes_q <- sizes |> filter(genome == alignment) |>  select(qName = chr, qSize = size)
  dir.create(glue("results/anchoring/{alignment}"), showWarnings = FALSE)
  write_lines("ScWAj4l\t1",
              file = glue("results/anchoring/{alignment}/weights.txt"))

  export_bed <- primary_borders  |>
    filter(nr %in% (border_ids |> filter(genome == alignment) |> pluck("ids") |> unlist())) |>
    arrange(qName, qStart) |>
    mutate(tPos = floor((tStart + tEnd)/2),
           qPos = floor((qStart + qEnd)/2)) |>
    dplyr::select(tName, qName, tPos, qPos, nr) |>
    left_join(sizes_t) |>
    left_join(sizes_q) |>
    arrange(tName, tPos) |>
    mutate(target_start = qPos - 1,
           querry_label = str_c(str_remove(tName, ";.*") |> str_replace("_","-"),":", sprintf("%.5f", tPos/tSize)),
           # querry_label = str_c(str_remove(qName, ";.*") |> str_replace("_","-"),":",querry),
           target_label = str_c(qName,":", qPos)) |>
    select(qName, target_start, qPos, querry_label, target_label)

  export_bed |>
   write_tsv(glue("results/anchoring/{alignment}/anchored_{str_remove(alignment,'_hardmasked')}.bed"), col_names = FALSE)

  tibble(alignment = alignment,
         p1 = list(p1),
         p2 = list(p2),
         export_bed = list(export_bed))
}

bed_exports <- best_hits_df |>
  pmap_dfr(plot_primary_borders)

# bed_exports$p1 |> wrap_plots()
# c(bed_exports$p1[2], bed_exports$p2[2]) |> wrap_plots()
# bed_exports$export_bed[[2]]
