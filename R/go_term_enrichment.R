# input:
# - "results/pinniped/go_terms/busco_id_to_go_term_bp.tsv"
# - "results/pinniped/busco_gerp_fst.tsv"
# - "results/pinniped/go_terms/thresholds.tsv"
# output:
# - "results/pinniped/go_terms/go_term_info.tsv"
# - "results/pinniped/go_terms/go_term_busco_stats.tsv"
# - "results/img/go_sub_graph.pdf"
# - "results/img/busco_go_term_2d_dens.pdf"
library(tidyverse)
library(prismatic)
library(patchwork)
library(topGO)
library(here)
library(ggraph)
library(tidygraph)
library(jsonlite)
library(glue)
library(ggdensity)
library(ggtext)

source(here("R/plot_defaults.R"))

go_bp <- readMappings(here("results/pinniped/go_terms/busco_id_to_go_term_bp.tsv"))
busco_data <- read_tsv(here("results/pinniped/busco_gerp_fst.tsv")) |>
  filter(all_min_2)

gerp_top_scores <- factor(as.numeric(busco_data$gerp_top_10)) |> set_names(busco_data$busco_id)
fst_scores <- factor(as.numeric(busco_data$fst_top_10)) |> set_names(busco_data$busco_id)
fst_scores <- fst_scores[!is.na(fst_scores)]

gene_intersect <- intersect(busco_data$busco_id, names(go_bp))

bp_genelist <- factor(rep(TRUE, length(names(go_bp))),
                      levels = c(FALSE, TRUE)) |>
  set_names(nm = names(go_bp))

genelist_gerp_top <- gerp_top_scores[gene_intersect]

gene_fst_intersect <- intersect(gene_intersect, names(fst_scores))
genelist_fst <- fst_scores[gene_fst_intersect]

all_thresholds <- read_tsv(here("results/pinniped/go_terms/thresholds.tsv"))

go_score_gerp_h <- function(allScore) {
  tr <- all_thresholds |> filter(type == "gerp", by == "BUSCO", prob == "99%") |>  pluck("threshold")
  return(allScore > tr)
}

go_score_fst_h <- function(allScore) {
  tr <- all_thresholds |> filter(type == "fst", by == "BUSCO", prob == "99%") |>  pluck("threshold")
  return(allScore > tr)
}

min_node_size <- 5
GOdata_gerp_top <- new("topGOdata",
                       ontology = "BP",
                       allGenes = genelist_gerp_top,
                       annot = annFUN.gene2GO,
                       nodeSize = min_node_size,
                       gene2GO = go_bp)

GOdata_fst <- new("topGOdata",
                   ontology = "BP",
                   allGenes = genelist_fst,
                   annot = annFUN.gene2GO,
                   nodeSize = min_node_size,
                   gene2GO = go_bp)


test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
result_gerp_top <- getSigGroups(GOdata_gerp_top, test.stat)
result_fst <- getSigGroups(GOdata_fst, test.stat)

check_scores <- \(obj, ttl = NULL){
  ggplot(tibble(score = obj@score),
         aes(x = score)) +
    geom_histogram(color = "darkgray",
                   fill = clr_alpha("darkgray"),
                   binwidth = .025) +
    labs(subtitle = ttl)
}

check_scores(result_gerp_top, "GERP") +
  check_scores(result_fst, "FST") &
  theme_minimal()

top_n <- 10
grph_gerp_top <- showSigOfNodes(GOdata_gerp_top, score(result_gerp_top), firstSigNodes = top_n, useInfo = 'all', swPlot = FALSE)
grph_fst <- showSigOfNodes(GOdata_fst, score(result_fst), firstSigNodes = top_n, useInfo = 'all', swPlot = FALSE)

plt_gprh <- \(grph, resultKS, ttl = NULL){
  test_graph <- igraph::graph_from_graphnel(grph$dag)

  tbl_grph <- tidygraph::as_tbl_graph(test_graph) |>
    left_join(tibble(name = names(resultKS@score),
                     score = resultKS@score)) |>
    arrange(score) |>
    mutate(signif = score < .01,
           top_go = row_number() <= top_n,
           go_rank = row_number())
  p <- tbl_grph |>
    ggraph(layout = "sugiyama") +
    geom_edge_diagonal2(aes(#alpha = after_stat(index),
                            edge_width = after_stat(index),
                            colour = node.score),
                        strength = 1.5,
                        alpha = .4
                        # color = "gray"
    ) +
    geom_node_point(aes(fill = score,
                        color = after_scale(clr_darken(fill)),
                        size = top_go,
                        shape = top_go),
                    stroke = .2) +
    geom_node_text(aes(label = go_rank,
                       filter = top_go, y = y -.35),
                   family = fnt_sel,
                   size = fnt_sz / ggplot2::.pt) +
    scale_edge_width_continuous(range = c(0,.5), guide = "none") +
    scale_fill_gradientn(colours = c(rev(clrs), "lightgray", "darkgray"),
                         values = c(0,.05,.1,1), limits = c(0, 1),
                         breaks = c(0,.05,.1, 1)) +
    scale_edge_color_gradientn(colours = c(rev(clrs), "lightgray", "darkgray"),
                               values = c(0,.05,.1,1),
                               breaks = c(0,.05,.1, 1),
                               guide = "none", limits = c(0, 1)) +
    scale_shape_manual(values = c(`TRUE` = 22, `FALSE` = 21), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1), guide = "none") +
    theme_ms() +
    labs(subtitle = ttl)

  list(tb = tbl_grph, p = p)
}

l_gerp_top <- plt_gprh(grph_gerp_top, result_gerp_top, "Sub-Ontology (GERP)")
l_fst <- plt_gprh(grph_fst, result_fst, "Sub-Ontology (*F<sub>ST<sub>*)")

p0 <- l_gerp_top$p +
  l_fst$p +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  guides(fill = guide_colorbar(title = "*p* value",
                               title.position = "top",
                               barheight = unit(3, "pt"),
                               barwidth = unit(.5, "npc"))) &
  theme(legend.position = "bottom",
        legend.title = element_markdown(hjust = .5),
        plot.subtitle = element_markdown(hjust = .5),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.tag = element_text(family = fnt_sel))

ggsave(filename = here("results/img/go_sub_graph.pdf"),
        plot = p0,
        width = 10,
        height = 5,
        device = cairo_pdf)

top_go_terms <- l_gerp_top$tb |>
  as_tibble() |>
  filter(top_go) |>
  select(go_term = name, gerp_top_p_val = score, gerp_top_rank = go_rank) |>
  full_join(
    l_fst$tb |>
      as_tibble() |>
      filter(top_go) |>
      select(go_term = name, fst_p_val = score, fst_rank = go_rank)
    ) |>
  select(go_term, ends_with("rank"), ends_with("val"))

get_go_details <- \(go_term){
  Sys.sleep(.005)
  tibble(fromJSON(glue("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_term}/complete"))$results)
  }

go_details <- top_go_terms |>
  left_join(map_dfr(top_go_terms$go_term, get_go_details), by = c(go_term = "id"))

go_details_filled <- go_details |>
  select(go_term, isObsolete:xRelations) |>
  mutate(across(go_term,
                .fns = list(gerp_top_rank = \(x){rank(score(result_gerp_top), ties.method = "min")[x]},
                            gerp_top_p_val = \(x){score(result_gerp_top)[x]},
                            fst_rank = \(x){rank(score(result_fst), ties.method = "min")[x]},
                            fst_p_val = \(x){score(result_fst)[x]}),
                .names = "{.fn}")) |>
  dplyr::select(go_term, ends_with("rank"), ends_with("p_val"), everything()) |>
  mutate(busco_id = genesInTerm(GOdata_gerp_top, go_term))

write_lines(x = glue("# Download report for GO term descriptions\n# Accessed at {Sys.time()} from https://www.ebi.ac.uk/QuickGO \n# Using the API in the form of https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/<go_term>/complete"),
            file = here("results/pinniped/go_terms/go_term_info.tsv"))

go_details_filled |>
  unnest(definition) |>
  select(go_term:text,aspect, usage) |>
  write_tsv(here("results/pinniped/go_terms/go_term_info.tsv"),col_names = TRUE,
            append = TRUE)

go_and_busco <- go_details_filled |>
  select( go_term, ends_with("_rank"), ends_with("_p_val"), busco_id ) |>
  unnest( busco_id ) |>
  left_join(busco_data)

go_and_busco |>
  write_tsv(here("results/pinniped/go_terms/go_term_busco_stats.tsv"))

gerp_median <- all_thresholds |> filter(type == "gerp", by == "BUSCO", prob == "50%") |>  pluck("threshold")
fst_median <- all_thresholds |> filter(type == "fst", by == "BUSCO", prob == "50%") |>  pluck("threshold")

clr_lab <- \(x1,x2){
  c2 <- clr_darken(clrs[[2]])
  glue("<span style='color:{clrs[[1]]}'>{x1}</span><br>&<br><span style='color:{c2}'>{x2}</span>")
}

pbd_0 <- go_and_busco |>
  filter(fst_rank <= top_n) |>
  ggplot(aes(x = gerp_rs_mean, y = fst_mean)) +
  geom_vline(xintercept = gerp_median,
             color = "gray", linetype = 3, linewidth = .2) +
  geom_hline(yintercept = fst_median,
             color = "gray", linetype = 3, linewidth = .2) +
  geom_hdr(data = busco_data,
           aes(fill = after_stat(probs)),
           probs = c(0.99, 0.9, 0.66),
           alpha= 1,
           show.legend = FALSE) +
  geom_richtext(data = tibble(gerp_rs_mean = rep((gerp_median + c(0, .1))/2, each = 2),
                          fst_mean = rep((fst_median + c(0, 1))/2,2),
                          label = c(clr_lab("unconserved","undiverged"),
                                    clr_lab("unconserved","diverged"),
                                            clr_lab("conserved","undiverged"),
                                                    clr_lab("conserved","diverged"))),
            aes(label = label),
            fill = NA, label.color = NA,
            label.padding = grid::unit(rep(0, 4), "pt"),
            family = fnt_sel,
            size = 1.2 * fnt_sz / ggplot2::.pt)+
  theme_ms()

pbd_t <- go_and_busco |>
  filter(gerp_top_rank <= top_n) |>
  ggplot(aes(x = gerp_rs_mean, y = fst_mean)) +
  geom_vline(xintercept = gerp_median,
             color = "gray", linetype = 3, linewidth = .2) +
  geom_hline(yintercept = fst_median,
             color = "gray", linetype = 3, linewidth = .2) +
  geom_hdr(data = busco_data,
           aes(fill = after_stat(probs)),
           probs = c(0.99, 0.9, 0.66),
           alpha= 1,
           show.legend = FALSE) +
  geom_hdr_lines(linewidth = .3,
                 probs = c(0.99, 0.9, 0.66),
                 aes(color = go_term == "GO:0051965"),
                 # color = clrs[[1]],
                 xlim = c(0, .1), ylim = c(0, 1)) +
  facet_grid(. ~ glue::glue("{str_pad(gerp_top_rank, width = 2, pad = 0)}: {go_term}")  ) +
  labs(subtitle = "High GERP Score GO Terms") +
  scale_color_manual(values = c(`TRUE` = "black",
                                `FALSE` = clrs[[1]]),
                     guide = "none") +
  theme_ms()  +
  theme(axis.title.x = element_blank())

pbd_f <- go_and_busco |>
  filter(fst_rank <= top_n) |>
  ggplot(aes(x = gerp_rs_mean, y = fst_mean)) +
  geom_vline(xintercept = gerp_median,
             color = "gray", linetype = 3, linewidth = .2) +
  geom_hline(yintercept = fst_median,
             color = "gray", linetype = 3, linewidth = .2) +
  geom_hdr(data = busco_data,
           aes(fill = after_stat(probs)),
           probs = c(0.99, 0.9, 0.66),
           alpha= 1,
           show.legend = FALSE) +
  geom_hdr_lines(linewidth = .3,
                 probs = c(0.99, 0.9, 0.66),
                 aes(color = go_term == "GO:0051965"),
                 # color = clr_darken(clrs[[2]],.2),
                 xlim = c(0, .1), ylim = c(0, 1)) +
  facet_grid(. ~ glue::glue("{str_pad(fst_rank, width = 2, pad = 0)}: {go_term}")) +
  scale_color_manual(values = c(`TRUE` = "black",
                                `FALSE` = clr_darken(clrs[[2]],.2)),
                     guide = "none") +
  labs(subtitle = "High *F<sub>ST</sub>* GO Terms") +
  theme_ms()

pp_2d_dens <- pbd_0 + (pbd_t + pbd_f + plot_layout(heights = c(.95, 1))) +
  plot_layout(guides = "collect",
              widths = c(.4, 1)) +
  plot_annotation(tag_levels = "a") &
  # scale_color_manual(values = c("gray", rev(clrs))) &
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,"Greys")[2:4]) &
  scale_x_continuous(breaks = c(0, .04, .08)) &
  scale_y_continuous(breaks = c(0, .5, 1)) &
  labs(x = "Average GERP RS Score per BUSCO",
       y = "Average *F<sub>ST</sub>* per BUSCO") &
  guides(alpha = guide_legend(title = "Density Mass",
                              override.aes = list(color = rgb(.1,.1,.1),
                                                  linewidth = .75))) &
  coord_cartesian(xlim = c(0, .1),
                  ylim = c(0, 1)) &
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = "lightgray", linewidth = .2),
        axis.title.y = element_markdown(),
        plot.subtitle = element_markdown(),
        plot.tag = element_text(family = fnt_sel),
        legend.margin = margin())

ggsave(plot = pp_2d_dens,
       filename = here("results/img/busco_go_term_2d_dens.pdf"),
       width = 9, height = 3.5, device = cairo_pdf)

go_details_filled |>
  unnest(definition) |>
  select(go_term,
         gerp_rank = gerp_top_rank,
         gerp_p = gerp_top_p_val,
         fst_rank,
         fst_p = fst_p_val,
         term = name,
         description = text) |>
  mutate(across(ends_with("p"), \(x){sprintf("%.4f",x)})) |>
  knitr::kable(format = "latex") |>
  str_remove_all("\\\\hline\n") |>
  write_lines(here("results/tab/top_go_terms.tex"))
