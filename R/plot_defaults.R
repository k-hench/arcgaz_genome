mm_to_inch <- 25.4
fwidth <- 178 / mm_to_inch
fhalfwidth <- 87 / mm_to_inch

fnt_sel <- "CMU Sans Serif"
fnt_sz <- 12 / ggplot2::.pt

clrs <- c("#1f4971", "#8cce28")
clrs_n <- \(n){scales::colour_ramp(colors = clrs)((0:(n-1))/(n-1))}
plt_lwd <- 0.15
