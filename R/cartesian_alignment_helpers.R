library(tidyverse)
library(ggforce)
library(prismatic)
library(here)
library(glue)
library(patchwork)
library(ggtext)
library(plyranges)
library(rlang)

# funs ---------------
#' read_size imports the size table for a given genome
#'
#'  This is an import function for the `.size` files
#'  created by `faSize --detailed` (in the `conda`
#'  environment `msa_align`).
#'
#' @param genome the file name prefix (without the `.size`)
#' @param y_base <numeric scalar> to attach an aribitrary y value to table
#' @param skip <numeric scalar> adds a bp offset between the concatinated scaffolds
#' @param order_by <"size"/any other string> whether to concatinate scaffolds in the order of decreasing lenght
#' @param manual_order_18 <named list, contains a vector of lenght 18 whos name matches `genome`> a manual re-ordering of the first 18 scaffolds
#'
read_size <- \(genome = "", y_base = 0, skip = 0, order_by = "size", manual_order_18 = NULL, n_first = 18){
  if(is.null(manual_order_18)){manual_order_18 <- list(1:n_first) |> set_names(nm = genome)}

  read_tsv(here::here("results", "genome", str_c(genome,".size")),
           col_names = c("chr", "size")) |>
    mutate(org_pos = row_number()) |>
    arrange(desc(size)) |>
    mutate(pre_manual = rank(-size),
           size_idx = if(order_by == "size"){
             row_number()
           } else {c(manual_order_18[[genome]], (n_first+1):length(chr))}) |>
    arrange(size_idx) |>
    mutate(end_with_skip = cumsum(size + skip),
           start = lag(end_with_skip, default =  0),
           end = start + size,
           mid = (start + end) /2,
           eo = row_number() %% 2,
           y_base = y_base,
           genome = genome,
           in_top = size_idx <= n_first)
}

#' Import a psl file
#'
#' @param file name of the `.psl` file (located within `results/psl`)
#'
get_psl <- \(file){
  vroom::vroom(here::here("results","psl", file),
               delim = "\t",
               col_names = c("matches", "misMatches", "repMatches", "nCount",
                             "qNumInsert", "qBaseInsert", "tNumInsert",
                             "tBaseInsert", "strand", "qName", "qSize", "qStart",
                             "qEnd", "tName", "tSize", "tStart", "tEnd",
                             "blockCount")) |>
    select( tName, tStart, tEnd, qName, qStart, qEnd, strand ) |>
    mutate(tSize = abs(tEnd - tStart),
           qSize = abs(qEnd - qStart))
}

#' Convert a table into wide format for links
#'
#' Converts link data into the needed format
#' for usage with `ggforce::geom_diagonal_wide()`.
#'
#' @param x_start x value for first genome (assuming vertical alignment)
#' @param x_end  x value for second genome (assuming vertical alignment)
#' @param ymin_start starting coordignate of the link on the first genome
#' @param ymax_start ending coordignate of the link on the first genome
#' @param ymin_end starting coordignate of the link on the second genome
#' @param ymax_end ending coordignate of the link on the secong genome
#' @param dir optional value for the strand direction
#' @param group the link ID
#'
diag_to_wide <- \(x_start = 0,
                  x_end = 1,
                  ymin_start, ymax_start,
                  ymin_end, ymax_end,
                  dir = "+",
                  group = ""){
  tibble(`1_x_min` = x_start,
         `2_x_max` = x_end,
         `3_x_max` = x_end,
         `4_x_min` = x_start,
         `1_y_min` = ymin_start,
         `2_y_min` = ymin_end,
         `3_y_max` = ymax_end,
         `4_y_max` = ymax_start,
         group = group,
         dir = dir) |>
    pivot_longer(cols = -c(group, dir)) |>
    separate(name,
             into = c("ord", "axis", "role"),
             convert = TRUE) |>
    arrange(ord) |>
    select(-role) |>
    pivot_wider(names_from = axis,
                values_from = value) |>
    select(-ord)
}

#' Converting a genome size table into a granges object
#'
#' @param sizes the genome size table (can contain several genomes)
#' @param genome_select the genome id of the target genome
#'
sizes_to_granges <- \(sizes, genome_select){
  sizes |>
    filter(genome == genome_select) |>
    select(seqnames = chr, end = size) |>
    mutate(start = 1) |> as_granges()
}

#' Import function to parse genome sizes and links
#'
#' @param genomes
#' @param n_alignments
#' @param align_min_length
#' @param manual_order_18
#' @param genome_width
#' @param skip
#' @param n_longest
#' @param order_by
#' @param coverage
#'
import_alignment <- \(genomes,
                      n_alignments = 20,
                      align_min_length = 0.2e6,
                      manual_order_18 = NULL,
                      genome_width = .1,
                      skip = .025,
                      n_longest = 50,
                      order_by = "name",
                      coverage = FALSE,
                      n_first = 18){

  ## helper funs ----
  cov_from_psl <- \(prefix, genome){
    psl |>
      select(seqnames = str_c(prefix, "Name"),
             start =  str_c(prefix, "Start"),
             end =  str_c(prefix, "End"),
             strand) |>
      as_granges() |>
      compute_coverage() |>
      as_tibble() |>
      left_join(sizes_18 |> select(seqnames = chr, g_start = start, eo = eo))  |>
      mutate(gg_start = start + g_start,
             gg_end = end + g_start,
             genome = genome) |>
      arrange(gg_start)
  }


  coverage_windows <- \(genome_idx){
    prefix = c("t", "q")[[genome_idx]]
    join_overlap_intersect(psl |>
                             select(seqnames = str_c(prefix, "Name"),
                                    start =  str_c(prefix, "Start"),
                                    end =  str_c(prefix, "End"),
                                    strand) |>
                             as_granges() |>
                             compute_coverage(),
                           as_granges(genome_windows[[genome_idx]])) |>
      as_tibble()  |>
      group_by(seqnames, window_id) |>
      summarise(start = min(start),
                end = max(end),
                avg_cov = sum(width * score) / sum(width)) |>
      ungroup() |>
      arrange(window_id) |>
      left_join(sizes_18 |> select(seqnames = chr, g_start = start, eo = eo))  |>
      mutate(gg_start = start + g_start,
             gg_end = end + g_start,
             gg_mid = (gg_start + gg_end) / 2,
             genome_idx = genome_idx,
             avg_cov_scl = avg_cov / max(avg_cov),
             avg_cov_25 = avg_cov / 5,
             avg_cov_25 = if_else(avg_cov_25 > 1, 1, avg_cov_25))
  }
  # -----------
  sizes <- tibble(genome = genomes,
                  y_base = seq_along(genomes)-1,
                  n_first = n_first) |>
    pmap_dfr(read_size, skip = 1e7,
             order_by = order_by,
             manual_order_18 = manual_order_18)  |>
    mutate(label_sign = 2 * ( .5 - (genome == genomes[[1]])))

  genome_summary <- sizes |>
    group_by(genome, y_base) |>
    summarise(n = n(),
              toal_length = max(end))

  sizes_18 <- sizes |>
    filter(in_top)

  psl <- get_psl(glue("slim_{genomes[[2]]}_on_{genomes[[1]]}.psl.gz"))

  psl_filtered <- psl |>
    filter(tName %in% sizes_18$chr,
           qName %in% sizes_18$chr) |>
    # arrange(desc(tSize)) |>
    group_by(tName) |>
    slice_max(order_by = tSize, n = n_longest) |>
    ungroup() |>
    # filter(tSize > align_min_length) |>
    left_join(sizes_18 |> select(tName = chr, tg_start = start, t_eo = eo)) |>
    left_join(sizes_18 |> select(qName = chr, qg_start = start, q_eo = eo)) |>
    mutate(ymin_start = tStart + tg_start,
           ymax_start = tEnd + tg_start,
           ymin_end = qStart + qg_start,
           ymax_end = qEnd + qg_start,
           group = row_number())

  psl_diag <- psl_filtered |>
    select(ymin_start:group, dir = q_eo) |>
    pmap_dfr(diag_to_wide)

  genome_with_skips <- tibble(y = c(0 - genome_width - skip,
                                    1 + skip) + (1/5) * genome_width,
                              start = -1e7,
                              end = c(max(sizes_18$end_with_skip[sizes_18$genome == genomes[[1]]]),
                                      max(sizes_18$end_with_skip[sizes_18$genome == genomes[[2]]])))

  if(coverage){
    genomes_granges <- genomes |> map(sizes_to_granges, sizes = sizes)
    genome_windows <- genomes_granges |>
      map(\(g){
        slidingWindows(g, width = 10e4, step = 5e4) |>
          as_tibble() |>
          mutate(window_id = row_number()) |>
          as_granges()})

    data_coverage <- tibble(prefix = c("t", "q"), genome = genomes) |>
      pmap_dfr(cov_from_psl)

    data_cov_windowed <- 1:2 |> map_dfr(coverage_windows)

    return(
      tibble(target = genomes[1],
             querry = genomes[2],
             sizes = list(sizes),
             sizes_18 = list(sizes_18),
             genome_summary = list(genome_summary),
             psl = list(psl),
             psl_filtered = list(psl_filtered),
             psl_diag = list(psl_diag),
             data_cov_windowed = list(data_cov_windowed),
             genome_with_skips = list(genome_with_skips) )
    )
  } else {
    return(
      tibble(target = genomes[1],
             querry = genomes[2],
             sizes = list(sizes),
             sizes_18 = list(sizes_18),
             genome_summary = list(genome_summary),
             psl = list(psl),
             psl_filtered = list(psl_filtered),
             psl_diag = list(psl_diag),
             genome_with_skips = list(genome_with_skips) )
    )
  }
}


psl_matrix <- \(target, querry, psl, sizes, sizes_18, ...){
  chr_sel <- sizes |> filter(genome == querry) |> pluck("chr")
  chr_sel_a <- sizes |> filter(genome == target) |> pluck("chr")

  psl_matrix <- psl |>
    filter(tName %in% sizes_18$chr & qName %in% sizes_18$chr) |>
    group_by(tName, qName) |>
    summarize(n = length(tName),
              alignment_length_t = sum(abs(tSize)),
              alignment_length_q = sum(abs(qSize))) |>
    ungroup() |>
    left_join(sizes_18 |> select(qName = chr, qSize = end)) |>
    left_join(sizes_18 |> select(tName = chr, tSize = end)) |>
    mutate(tName = factor(tName, chr_sel_a),
           qName = factor(qName, chr_sel),
           t_align_fract = alignment_length_t / abs(tSize),
           q_align_fract = alignment_length_q / abs(qSize))

  tibble(target = target,
         querry = querry,
         psl_matrix = list(psl_matrix))
}

simplify_names <- \(chr){
  chr |>
    str_remove(";.*") |>
    str_remove("ScDpJTZ_") |>
    str_remove("Sc3htSU_") |>
    str_remove("CAAAJK01000") |>
    str_remove("mscaf_") |>
    str_remove("NC_045")
}

plot_matrix <- \(target, querry, psl_matrix,
                 tile_min_sz = 0, tile_scl = 1,
                 cmin = 0, cmax = 1){
      p1 <- psl_matrix |>
        mutate(across(ends_with("Name"), simplify_names)) |>
        ggplot(aes(y = tName, x = qName)) +
        geom_tile(aes(fill = t_align_fract,
                      width = tile_min_sz + q_align_fract * tile_scl,
                      height = after_stat(width),
                      color = after_scale(clr_darken(fill)))) +
        scale_fill_gradientn(limits = c(cmin, cmax),
                             colours = c(#"gray90",
                                         clrs[2] , "black"), na.value = "red") +
        guides(fill = guide_colorbar("Alginment Fraction (Target)",
                                     title.position = "top",
                                     barwidth = unit(.4,"npc"),
                                     barheight = unit(5,"pt"))) +
        labs(x = glue("{querry} (q)"),
             y = glue("{target} (t)"))

      p2 <- psl_matrix |>
        mutate(across(ends_with("Name"), simplify_names)) |>
        ggplot(aes(y = tName, x = qName)) +
        geom_tile(aes(fill = q_align_fract,
                      width = tile_min_sz + q_align_fract * tile_scl,
                      height = after_stat(width),
                      color = after_scale(clr_darken(fill)))) +
        scale_fill_gradientn(limits = c(cmin, cmax),
                             colours = c(#"gray90",
                                         clrs[1] , "black"), na.value = "red") +
        guides(fill = guide_colorbar("Alginment Fraction (Querry)",
                                     title.position = "top",
                                     barwidth = unit(.4,"npc"),
                                     barheight = unit(5,"pt"))) +
        labs(x = glue("{querry} (q)"),
             y = glue("{target} (t)"))

      p <- p1 / p2 &
        coord_equal() &
        theme_minimal(base_family = fnt_sel) &
        theme(axis.text.x = element_text(angle = 90),
              legend.position = "bottom")

      tibble(p1 = list(p1),
             p2 = list(p2),
             p = list(p))
}
