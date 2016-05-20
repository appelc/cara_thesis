##########################################################
### Exploring body mass change
### (1) create figure with trends for each animal
### (2) linear model to find informed seasonal cutoff
### (3) t-test for significance based on seasonal cutoff
##########################################################

library(googlesheets)
library(ggplot2)
library(dplyr)

## load from Google Drive
gs_ls()
wts <- gs_title("Porc weights & captures")
porc.wts <- data.frame(gs_read(ss=wts, ws="Weights", is.na(TRUE), range=cell_cols(1:5)))
colnames(porc.wts) <- c("id", "date", "kg", "sex", "month")
head(porc.wts)

## format correctly
porc.wts$id <- as.factor(porc.wts$id)
porc.wts$date <- as.Date(porc.wts$date, "%m/%d/%Y")
porc.wts$sex <- as.factor(porc.wts$sex)
porc.wts$month <- as.factor(porc.wts$month)
porc.wts$kg <- as.numeric(porc.wts$kg)

## subset males/females
m.wts <- porc.wts[porc.wts$sex == "M",]
m.wts <- droplevels(m.wts)
f.wts <- porc.wts[porc.wts$sex == "F",]
f.wts <- droplevels(f.wts)

## get rid of singletons (Bowie and Badger)?
f.wts <- f.wts[f.wts$id != "15.05" & f.wts$id != "16.16",]

##########################################################
### (1) create figure with trends for each animal
##########################################################

ggplot(data=porc.wts, aes(x = date, y = wt, group = id)) +
          geom_point(data = f.wts, aes(y = f.wts$wt, colour = 'Female', shape = 'Female'), size = 5) + 
          geom_line(data = f.wts, colour = 'black', size = 2) +
          geom_point(data = m.wts, aes(y = m.wts$wt, colour = 'Male', shape = 'Male'), size = 5) + 
          geom_line(data = m.wts, colour = 'darkgray', size = 2, lty = 2) +
          ylab('Body mass (kg)') + xlab('Date') + 
          scale_y_continuous(limits = c(4,11)) +
          theme(axis.text=element_text(size=20, colour = 'black'),
                axis.title=element_text(size=24, colour = 'black'),
                axis.line.x = element_line(size = 1, colour = 'black'),
                axis.line.y = element_line(size = 1, colour = 'black'),
                axis.ticks = element_line(size = 1, colour = 'black'),
                axis.ticks.length = unit(.3, 'cm'),
                panel.background = element_rect(fill = 'white'),
                legend.background = element_rect(colour = 'white'),
                legend.key = element_rect(fill = 'white'),
                legend.text = element_text(size = 20),
                legend.key.size = unit(1.5, "cm")) +     
          scale_colour_manual("", 
                      breaks = c('Female', 'Male'),
                      values = c('black', 'darkgray')) +
          scale_shape_manual("",
                      breaks = c('Female', 'Male'),
                      values = c(15, 16))


##########################################################
### (1b) descriptive stats
##########################################################

## weight at capture (summer): males vs. females
sum.f <- porc.wts[c(1, 2, 5, 7, 8, 9, 10, 12, 13),]
sum.m <- porc.wts[c(3, 4, 6, 11, 15),]

mean(sum.f$kg)
(sd(sum.f$kg)) / (sqrt(nrow(sum.f)))

mean(sum.m$kg)
(sd(sum.m$kg)) / (sqrt(nrow(sum.m)))  

t.test(sum.f$kg, sum.m$kg, paired = FALSE) # no sig diff between weights at capture

## weight at capture (winter): males vs. females
win.f <- porc.wts[c(43, 44),]
win.m <- porc.wts[c(41, 45),]

mean(win.f$kg)
(sd(win.f$kg)) / (sqrt(2))

mean(win.m$kg)
(sd(win.m$kg)) / (sqrt(2))

t.test(win.f$kg, win.m$kg, paired = FALSE) # no sig diff between weights at capture

## weight at capture (all together)
cap.wts <- rbind(sum.f, sum.m, win.f, win.m)
mean(cap.wts$kg)
(sd(cap.wts$kg)) / (sqrt(nrow(cap.wts)))



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
porc.wts15 <- porc.wts[porc_2015,] ## then only keep the rows with TRUE
porc.wts15 <- droplevels(porc.wts15)

## subset m/f again for animals captured in 2015
m.wts15 <- porc.wts15[porc.wts15$sex == "M",]
m.wts15 <- droplevels(m.wts15)
f.wts15 <- porc.wts15[porc.wts15$sex == "F",]
f.wts15 <- droplevels(f.wts15)

##### FROM HERE FORWARD, NEED TO SELECT A CUTOFF DATE

## separate summer/winter weights
## could do this all iteratively again, I suppose
s.wts <- porc.wts15[porc.wts15$date < as.Date("2015-11-01"),]
w.wts <- porc.wts15[porc.wts15$date >= as.Date("2015-11-01"),]

## if multiple summer weights, need to take the average
## could use "table" to ask which ones have >1 summer weight (it's only 15.12), 
## but it can still calculate mean based on 1 weight and will be simpler

ids <- unique(s.wts$id)
s.avgs <- NULL

for(i in ids){
  s.avg.i <- mean(s.wts$kg[s.wts$id == i])
  sex.i <- s.wts$sex[s.wts$id == i]
  wt.id <- data.frame(i, s.avg.i, sex.i)
  colnames(wt.id) <- c("id", "s_avg", "sex")
  s.avgs <- rbind(s.avgs, wt.id)
}

dups <- duplicated(s.avgs[, ("id")])
s.avgs <- s.avgs[!dups,]

## if multiple winter weights, need to take the average
## several do have multiple winter weights

ids <- unique(w.wts$id)
w.avgs <- NULL

for(i in ids){
    w.avg.i <- mean(w.wts$kg[w.wts$id == i])
    wt.id <- data.frame(i, w.avg.i)
    colnames(wt.id) <- c("id", "w_avg")
    w.avgs <- rbind(w.avgs, wt.id)
}

w.avgs$w_avg <- round(w.avgs$w_avg, 3)

## now combine them 
## hmm, maybe I could even leave in ones that only have winter until this point
## need dplyr for 'left_join'
avg.wts <- left_join(s.avgs, w.avgs) ## will get a warning; see next step
avg.wts$id <- as.factor(avg.wts$id) 

## now get rid of rows with NA
avg.wts <- avg.wts[!is.na(avg.wts$w_avg),]
avg.wts <- avg.wts[,c("id", "sex", "s_avg", "w_avg")]
avg.wts

## now do statistics!
t.test(avg.wts$s_avg, avg.wts$w_avg, paired=TRUE)

## females
avg.wts.f <- avg.wts[avg.wts$sex == 'F',]
t.test(avg.wts.f$s_avg, avg.wts.f$w_avg, paired = TRUE)

## males
avg.wts.m <- avg.wts[avg.wts$sex == 'M',]
t.test(avg.wts.m$s_avg, avg.wts.m$w_avg, paired = TRUE)

## 
mean_s <- mean(avg.wts$s_avg)
mean_s # 8.497 kg ## why do I get a different answer now? 8.192 kg (05/05/16)
mean_w <- mean(avg.wts$w_avg)
mean_w # 7.79 kg ## 7.403 kg (05/05/16)
means <- cbind(mean_s, mean_w)

par(mfrow=c(1,1))
barplot(means, ylab="Body mass (kg)")
write.csv(avg.wts, "avg.wts.csv")


