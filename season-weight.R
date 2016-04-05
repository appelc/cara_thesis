##########################################################
### Exploring body mass change
### (1) create figure with trends for each animal
### (2) linear model to find informed seasonal cutoff
### (3) t-test for significance based on seasonal cutoff
##########################################################

library(googlesheets)
library(ggplot2)
library(dplyr)

setwd("C:/Users/Cara/Documents/__RESEARCH/ANALYSIS/033016")

## load from Google Drive
gs_ls()
wts <- gs_title("Porc weights & captures")
wts
porc.wts <- data.frame(gs_read(ss=wts, ws="Weights", is.na(TRUE), range=cell_cols(1:5)))
colnames(porc.wts) <- c("id", "date", "wt", "sex", "month")
head(porc.wts)

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

##########################################################
### (2) linear model to find informed seasonal cutoff
##########################################################

## doesn't really make sense as a LM now that I think about it;
## maybe random effect or ANOVA.
## how could I look at weight change to find a seasonal switch?

lm1 <- lm(wt ~ date, data=porc.wts)
summary(lm1)

plot(wt ~ date, data=porc.wts)
abline(lm1)

##########################################################
### (3) t-test for significance based on seasonal cutoff
##########################################################

## keep only porcs 15.01 - 15.14 (b/c no summer data on others)
## (maybe I could actually keep them all until the end...)
porc_2015 <- grepl('^15', porc.wts$id) ## ask which ones start with "15" 
porc_2015 # (TRUE/FALSE)
porc.wts <- porc.wts[porc_2015,] ## then only keep the rows with TRUE
porc.wts <- droplevels(porc.wts)

## subset m/f again for animals captured in 2015
f.wts15 <- subset(porc.wts, sex="F")
f.wts15 <- droplevels(f.wts15)
m.wts15 <- subset(porc.wts, sex="M")
m.wts15 <- droplevels(f.wts15)

##### FROM HERE FORWARD, NEED TO SELECT A CUTOFF DATE

## separate summer/winter weights
## could do this all iteratively again, I suppose
s.wts <- subset(porc.wts, date < as.Date("2015-10-01"))
w.wts <- subset(porc.wts, date >= as.Date("2015-10-01"))

## if multiple summer weights, need to take the average
## could use "table" to ask which ones have >1 summer weight (it's only 15.12), 
## but it can still calculate mean based on 1 weight and will be simpler

ids <- unique(s.wts$id)
s.avgs <- NULL

for(i in ids){
  s.avg.i <- mean(s.wts$wt[(which(s.wts$id == i))])
  sex.i <- s.wts$sex[which(s.wts$id == i)]
  wt.id <- data.frame(i, s.avg.i, sex.i)
  s.avgs <- rbind(s.avgs, wt.id)
}

colnames(s.avgs) <- c("id", "s_wt", "sex")
dups <- duplicated(s.avgs[, ("id")])
s.avgs <- s.avgs[!dups,]

## if multiple winter weights, need to take the average
## several do have multiple winter weights

ids <- unique(w.wts$id)
w.avgs <- NULL

for(i in ids){
    w.avg.i <- mean(w.wts$wt[(which(w.wts$id == i))])
    wt.id <- data.frame(i, w.avg.i)
    w.avgs <- rbind(w.avgs, wt.id)
}

colnames(w.avgs) <- c("id", "w_wt")
w.avgs$w_wt <- round(w.avgs$w_wt, 2)

## now combine them 
## hmm, maybe I could even leave in ones that only have winter until this point
avg.wts <- left_join(s.avgs, w.avgs) ## will get a warning; see next step
avg.wts$id <- as.factor(avg.wts$id) 

## now get rid of rows with NA
avg.wts <- avg.wts[!is.na(avg.wts$w_wt),]
avg.wts <- avg.wts[,c("id", "sex", "s_wt", "w_wt")]
avg.wts

## now do statistics!
t.test(avg.wts$s_wt, avg.wts$w_wt, paired=TRUE)

mean_s <- mean(avg.wts$s_wt)
mean_s # 8.497 kg
mean_w <- mean(avg.wts$w_wt)
mean_w # 7.79 kg
means <- cbind(mean_s, mean_w)

par(mfrow=c(1,1))
barplot(means, ylab="Body mass (kg)")
write.csv(avg.wts, "avg.wts.csv")


