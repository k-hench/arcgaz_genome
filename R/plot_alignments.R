#!/usr/bin/env Rscript
library(tidyverse)
library(circlize)
library(glue)
library(prismatic)

fnt_sel <- "Josefin sans"

get_psl <- \(idx){
  read_tsv(glue("results/psl/genome_{idx}.psl"),
           col_names = c("matches", "misMatches", "repMatches", "nCount",
                         "qNumInsert", "qBaseInsert", "tNumInsert",
                         "tBaseInsert", "strand", "qName", "qSize", "qStart",
                         "qEnd", "tName", "tSize", "tStart", "tEnd",
                         "blockCount", "blockSizes", "qStarts", "tStarts")) |> 
    select( tName, tStart, tEnd, qName, qStart, qEnd, strand ) |> 
    mutate( tName = str_c("g1_",tName),
            qName = glue("g{idx}_{qName}"))
}

data <- 2:3 |> map_dfr(get_psl) 

get_size <- \(idx){
  read_tsv(glue("results/genome/genome_{idx}.size"),
           col_names = c("chr", "end")) |> 
    mutate(chr = glue("g{idx}_{chr}"),
           start = 0) |> 
    select(chr, start, end)
  }

sizes <- 1:3 |> map_dfr(get_size)
clrs <- c(g1 = "gray30", 
          g2 = clr_darken("#f46d43",.3),
          g3 = clr_darken("#66c2a5",.3))

clrs_strand <- c(clrs, clr_darken(clrs, .5)) |> 
  set_names(nm = c(str_c(names(clrs),'+'),
                   str_c(names(clrs),'-')))

get_bed <- \(idx){read_tsv(glue("results/bed/inversions_{idx}.bed"))}
regions <- 1:3 |> map_dfr(get_bed)
clr_al <- clrs_strand[str_c(str_sub(data$qName,1,2), data$strand)]

svg(filename = "img/alignment.svg", width = 5, height = 5)
circos.clear()
col_text <- "white"
circos.par("track.height" = 0.8, gap.degree = 5, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors = sizes$chr, 
                  xlim = sizes[,2:3] |>
                    as.matrix())
  
# genomes
circos.track(ylim=c(0, 1), 
             panel.fun=function(x, y) {
               chr=CELL_META$sector.index
               xlim=CELL_META$xlim
               ylim=CELL_META$ylim
               circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
                           facing = "bending.inside",
                           niceFacing = TRUE,
                           family = fnt_sel)},
             bg.col= rep(clrs, each = 2),
             bg.border = FALSE, 
             track.height = 0.06)

# genomes x axis
x_scale <- 1e3
brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)*x_scale * 10
circos.track(track.index = get.current.track.index(), 
             panel.fun=function(x, y) {
               circos.axis(h="top", major.at=brk, labels = round(brk/x_scale, 1), 
                           labels.cex = 0.4,
                           col = "black", labels.col="black",
                           lwd = 0.7, labels.facing="clockwise")
             }, bg.border = FALSE)


circos.genomicTrack(regions,
                    stack = TRUE,
                    panel.fun = function(region_in, value, ...) {
                      circos.genomicRect(region_in, value,
                                         col = clr_alpha(clrs[str_sub(CELL_META$sector.index,1,2)], .7),
                                         border = NA, ...)
                    },
                    bg.border = FALSE, 
                    track.height = 0.03)


circos.genomicLink(data[1:3], data[4:6],
                   col = clr_alpha(clr_al,.4),
                   border = clr_darken(clr_al))
dev.off()
