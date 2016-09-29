####################################
### Porcupine habitat selection analysis
### In ArcMap, used spatial join to calculate loc points within each veg polygon
### Following analysis directions from adehabitatHS package vignette
####################################

install.packages("data.table")
library(data.table)
library(adehabitatHS)

setwd("C:/Users/Cara/Documents/__RESEARCH/_DATA/Analysis/0304516")

####################################
#### DESIGN 1 HABITAT SELECTION
####################################

## SUMMER
porcpts.sum <- read.csv("join_output_sum.txt")
head(porcpts.sum)
# extract just the relevant columns:
porcpts.sum <- data.frame(porcpts.sum$Class_2, porcpts.sum$Area_1, porcpts.sum$Sum_Count)
colnames(porcpts.sum) <- c("class", "area", "count")
head(porcpts.sum)
porcpts.sum$area <- as.numeric(porcpts.sum$area)
porcpts.sum$count <- as.numeric(porcpts.sum$count)
## create a summary table with sum(area) and sum(count) for each class:
sum.pts <- aggregate(cbind(area, count) ~ class, data=porcpts.sum, FUN = sum) 
## add columns with calculations for % used and % avail:
sum.pts$used <- sum.pts$count / (sum(sum.pts$count)) 
sum.pts$avail <- sum.pts$area / (sum(sum.pts$area))
head(sum.pts)
## write to a CSV
write.csv(sum.pts, "porc_used_avail_sum.csv")

## WINTER
porcpts.win <- read.csv("join_output_win.txt")
head(porcpts.win)
porcpts.win <- data.frame(porcpts.win$Class_2, porcpts.win$Area_1, porcpts.win$Sum_Count)
colnames(porcpts.win) <- c("class", "area", "count")
head(porcpts.win)
porcpts.win$area <- as.numeric(porcpts.win$area)
porcpts.win$count <- as.numeric(porcpts.win$count)
win.pts <- aggregate(cbind(area, count) ~ class, data=porcpts.win, FUN = sum)
win.pts$used <- win.pts$count / (sum(win.pts$count))
win.pts$avail <- win.pts$area / (sum(win.pts$area))
head(win.pts)
write.csv(win.pts, "porc_used_avail_win.csv")

## write one CSV for ease of making figures, etc.
used_avail <- cbind(sum.pts, win.pts$used)
write.csv(used_avail, "porc_used_avail_fig.csv")
## don't forget to save it as .xlsx when making figure

## BAR PLOTS
## revisit this (fix width, order, etc.)
par(mfrow=c(1,1))
barplot(sum.pts$avail, sum.pts$used, names=sum.pts$class, beside=TRUE, col=c("red", "blue"), axes=TRUE)

####################################
## DESIGN I: Selection Ratios
####################################

## make vectors for USED sum/win points
us <- sum.pts$count
names(us) <- sum.pts$class

uw <- win.pts$count
names(uw) <- sum.pts$class

## make vector for AVAILABLE percentage
## should be same sum/win!
tav <- sum.pts$avail 
names(tav) <- sum.pts$class

## SUMMER SELECTION RATIOS
class(us) <- NULL
class(tav) <- NULL
wi_s <- widesI(us, tav)
wi_s
plot(wi_s)
names(wi_s)
s_output <- data.frame(wi_s$used.prop, wi_s$se.used, wi_s$avail.prop, wi_s$se.avail, wi_s$wi, wi_s$se.wi, wi_s$chisquwi, wi_s$Bi, wi_s$comparisons)
write.csv(s_output, "wi_s.csv")

## WINTER SELECION RATIOS
class(uw) <- NULL
class(tav) <- NULL
wi_w <- widesI(uw, tav)
wi_w
plot(wi_w)
w_output <- data.frame(wi_w$used.prop, wi_w$se.used, wi_w$avail.prop, wi_w$se.avail, wi_w$wi, wi_w$se.wi, wi_w$chisquwi, wi_w$Bi, wi_w$comparisons)
write.csv(w_output, "wi_w.csv")

####################################
#### DESIGN 2 HABITAT SELECTION
####################################

porcpts <- read.csv("join_output_sum.txt")
ids <- porcpts.sum
head(porcpts)

## Need to figure out how to get matrix with loc by Animal ID
## could separate shapefiles and then do spatial join but way too time-consuming

## spatial join in R? 
## try "point.in.polygon." in package sp:
library(sp)
