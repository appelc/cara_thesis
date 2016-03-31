#############################
### Trying to find an informed seasonal cutoff
### Model weight against season
############################

library(googlesheets)

setwd("C:/Users/Cara/Documents/__RESEARCH/ANALYSIS/033016")

## load location points
gs_ls()
wts <- gs_title("Porc weights & captures")
wts
porc.wts <- data.frame(gs_read(ss=wts, ws="Weights", is.na(TRUE), range=cell_cols(1:5)))
colnames(porc.wts) <- c("id", "date", "wt", "sex", "month")

porc.wts$id <- as.factor(porc.wts$id)
porc.wts$date <- as.Date(porc.wts$date, "%m/%d/%Y")
porc.wts$sex <- as.factor(porc.wts$sex)
porc.wts$month <- as.factor(porc.wts$month)

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
