library(ggplot2)
library(gridExtra)
library(grid)

my_fun <- function(id){
  ggplot() + ggtitle(paste(id, "hours-feed"))
}

pl <- lapply(seq_len(10), my_fun)

lg <- tableGrob(c("", "26ppm", "39ppm"), theme= ttheme_minimal())
rg <- arrangeGrob(grobs = pl, ncol=5,
                  top = textGrob("Fish 11 - feed",gp=gpar(fontsize=18)))

grid.newpage()
grid.draw(cbind(lg, rg, size = "last"))
