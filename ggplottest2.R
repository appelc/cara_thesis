data(iris)
a <- ggplot(iris, aes(x=Species, y=Petal.Width)) + geom_boxplot(color="black") + ylab(expression(Foo~Bar~(g~cm^{-3})))
b <- ggplot(iris, aes(x=Species, y=Petal.Length*100)) + geom_boxplot(color="black") + ylab("foobar (mm)")
c <- ggplot(iris, aes(x=Species, y=Sepal.Width)) + geom_boxplot(color="black") + ylab("foobar (%)")
d <- ggplot(iris, aes(x=Species, y=log10(Sepal.Length))) + geom_boxplot(color="black") + ylab("foobar (cm)")

plots <- list(a,b,c,d)
grobs = lapply(plots, ggplotGrob)
g = do.call(rbind, c(grobs, size="first"))

g$widths = do.call(unit.pmax, lapply(grobs, "[[", "widths"))
grid.newpage()
grid.draw(g)


plots1 <- list(a,b)
plots2 <- list(c,d)
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)

g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(rbind, c(list(g1,g2), size="first")) #this does not work

grid.newpage()
grid.draw(g3)



# Tweak the margins to use up empty space.  Margins: Top, Right, Bottom, Left
a1 <- a + theme(axis.title.x = element_blank(), 
                axis.text.x = element_blank(), 
                axis.ticks.x= element_blank(), 
                plot.margin= unit(c(1, 1, -0.5, 0.5), "lines") ) 
b1 <- b + theme(axis.title.x = element_blank(), 
                axis.text.x = element_blank(), 
                axis.ticks.x= element_blank(), 
                plot.margin= unit(c(1, 1, -0.5, 0.5), "lines") ) 
c1 <- c + theme(plot.margin= unit(c(0, 1, 0.5, 0.5), "lines") )  
d1 <- d + theme(plot.margin= unit(c(0, 1, 0.5, 0.5), "lines") )  

grobz <- lapply(list(a1, b1, c1, d1), ggplotGrob)
grobz.plot <- arrangeGrob( grobs = list(rbind(grobz[[1]], grobz[[3]], size = "last"),
                                        rbind(grobz[[2]], grobz[[4]], size = "last")),
                           ncol = 2)
grbzplot <- arrangeGrob(grobs = list(rbind(a1, b1, size = 'last'),
                                     rbind(c1, d1, size = 'last')),
                        ncol = 2)

grid.draw(grobz.plot)
