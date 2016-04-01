##########################################################
### Exploring body mass change
### (1) create figure with trends for each animal
### (2) linear model to find informed seasonal cutoff
### (3) t-test for significance based on seasonal cutoff
##########################################################

library(googlesheets)
library(ggplot2)

setwd("C:/Users/Cara/Documents/__RESEARCH/ANALYSIS/033016")

## load from Google Drive
gs_ls()
wts <- gs_title("Porc weights & captures")
wts
porc.wts <- data.frame(gs_read(ss=wts, ws="Weights", is.na(TRUE), range=cell_cols(1:5)))
colnames(porc.wts) <- c("id", "date", "wt", "sex", "month")

## format correctly
porc.wts$id <- as.factor(porc.wts$id)
porc.wts$date <- as.Date(porc.wts$date, "%m/%d/%Y")
porc.wts$sex <- as.factor(porc.wts$sex)
porc.wts$month <- as.factor(porc.wts$month)

## subset males/females
m.wts <- subset(porc.wts, sex == "M")
m.wts <- droplevels(m.wts)
f.wts <- subset(porc.wts, sex == "F")
f.wts <- droplevels(f.wts)

##########################################################
### (1) create figure with trends for each animal
##########################################################

ggplot(data=porc.wts, aes(x=date, y=wt, group=id)) +
          geom_point(colour = 'black', size = 3) + geom_line(colour = 'black', size = 1) +
          geom_point(data = m.wts, aes(y = m.wts$wt),
             colour = 'grey', size = 3) + geom_line(data=m.wts, colour = 'grey', size = 1)

###################################

lm1 <- lm(wt ~ date, data=porc.wts)
summary(lm1)

plot(lm1)
plot(wt ~ date, data=porc.wts)
abline(lm1)

ids <- unique(porc.wts$id)
frame()

for(i in ids){
  porc.i <- subset(porc.wts, id == i)
  droplevels(porc.i)
  lm.i <- lm(wt ~ date, data=porc.i)
  plot(wt ~ date, data=porc.i, main=i)
  abline(lm.i)
}

