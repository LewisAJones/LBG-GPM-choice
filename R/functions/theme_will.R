theme_will <- function(...) {
  theme(axis.ticks = element_line(color = "black", linewidth = 1),
        axis.line = element_blank(),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(colour = "black", size = 12),
        panel.border = element_rect(linetype = "solid", colour = "black",
                                    fill = NA, linewidth = 2),
        panel.background = element_blank(),
        ...)
}