library(tidyverse)
library(glue)
library(patchwork)
fnt_sel <- "Josefin sans"

get_size <- \(x){
  read_tsv(glue("results/genome/genome_{x}.size"),
           col_names = c("chr", "size")) |> 
    mutate(end = cumsum(size),
           start = lag(end,default = 0),
           seq = glue("seq{x}"))
}

clrs <- c("gray30",
          "#f46d43",
           "#66c2a5") 

plot_lastz <- \(seq1, seq2, bp_scle = 1e-3){
  data <- read_tsv(glue("results/lastz/dplot/genome_{seq2}.tsv")) |> 
    filter(complete.cases(seq1)) |> 
    mutate(segment_idx = (row_number() - 1)%/%2 + 1,
           segment_role = c("start", "end")[(row_number() - 1)%%2 + 1]) |> 
    pivot_wider(values_from = starts_with("seq"),
                id_cols = segment_idx,
                names_from = segment_role)
  
  sizes <- c(1, 2) |> 
    map_dfr(get_size) |> 
    pivot_wider(values_from = size:start,
                id_cols = chr,
                names_from = seq)
  
  data |> 
    ggplot() +
    geom_linerange(data = sizes,
                   aes(y = -Inf,
                       xmin = start_seq1 * bp_scle,
                       xmax = end_seq1 * bp_scle,
                       color = chr),
                   linewidth = 3) +
    geom_linerange(data = sizes,
                   aes(x = -Inf,
                       ymin = start_seq2 * bp_scle, 
                       ymax = end_seq2 * bp_scle,
                       color = chr),
                   linewidth = 3) +
    geom_vline(data = sizes, aes(xintercept = end_seq1 * bp_scle), linetype = 3)+
    geom_hline(data = sizes, aes(yintercept = end_seq2 * bp_scle), linetype = 3)+
    geom_segment(aes(x = seq1_start * bp_scle, xend = seq1_end * bp_scle,
                     y = seq2_start * bp_scle, yend = seq2_end * bp_scle),
                 linewidth = 2, color = clrs[seq2]) +
    scale_color_manual(values = c("gray20", "gray75"), guide = "none") +
    labs(x = glue("seq{seq1} (kb)"), y = glue("seq{seq2} (kb)")) +
    theme_minimal()
}

p0 <- plot_lastz(1, 2) +
plot_lastz(1,3) +
  plot_annotation(subtitle = "LASTZ alignments") &
  theme_minimal(base_family = fnt_sel) &
  theme(plot.subtitle = element_text( hjust= .5 ))

ggsave(plot = p0, 
       "img/lastz_alignments.svg",
       width = 7, height = 3.5)
