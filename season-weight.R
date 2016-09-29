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

## get rid of singletons for figure (04, 05, 16)?
f.wts <- f.wts[f.wts$id != "15.05" & f.wts$id != "16.16",]
f.wts <- droplevels(f.wts)
m.wts <- m.wts[m.wts$id != '15.04',]
m.wts <- droplevels(m.wts)

##########################################################
### (1) create figure with trends for each animal
##########################################################
#par(mar = c())
ggplot(data=porc.wts, aes(x = date, y = kg, group = id)) +
  scale_y_continuous(limits = c(4,11)) +   
  scale_x_date() +
  geom_vline(xintercept = as.numeric(porc.wts$date[26]), linetype='dashed', lwd = 1.5) + # the 26th entry corresponds to the first Nov entry
          geom_vline(xintercept = as.numeric(porc.wts$date[46]), linetype='dashed', lwd = 1.5) + # the 46th entry corresponds to the first March entry        
          geom_point(data = f.wts, aes(y = f.wts$kg, colour = 'Female', shape = 'Female'), size = 5) + 
          geom_line(data = f.wts, colour = 'black', size = 2) +
          geom_point(data = m.wts, aes(y = m.wts$kg, colour = 'Male', shape = 'Male'), size = 5) + 
          geom_line(data = m.wts, colour = 'darkgray', size = 2, lty = 2) +
          ylab('Body mass (kg)') + xlab('Date') + 
        theme(axis.text=element_text(size=18, colour = 'black'),
                axis.title.x=element_text(size=22, colour = 'black', margin(20,0,0,0)), ## fix margins
                axis.title.y=element_text(size=22, color = 'black', margin(0,20,0,0)),
                axis.line.x = element_line(size = 1, colour = 'black'),
                axis.line.y = element_line(size = 1, colour = 'black'),
                axis.ticks = element_line(size = 1, colour = 'black'),
                axis.ticks.length = unit(.3, 'cm'),
                panel.background = element_rect(fill = 'white'),
                legend.background = element_rect(colour = 'white'),
                legend.key = element_rect(fill = 'white'),
                legend.text = element_text(size = 18),
                legend.key.size = unit(1.5, "cm")) +     
          scale_colour_manual("", 
                      breaks = c('Female', 'Male'),
                      values = c('black', 'darkgray')) +
          scale_shape_manual("",
                      breaks = c('Female', 'Male'),
                      values = c(15, 16))

## I tried 'geom_rect' to shade the area between vertical lines (winter) but it didn't work

##########################################################
### (1b) descriptive stats
##########################################################

## weight at capture (summer): males vs. females
cap <- porc.wts[match(unique(porc.wts$id), porc.wts$id),] ## these are the 1st wts for each porc
cap.f <- cap[cap$sex == 'F',]
cap.f.s <- cap.f[-c(10:11),] ## females captured in summer
cap.f.w <- cap.f[c(10:11),] ## females captured in winter
cap.m <- cap[cap$sex == 'M',]
cap.m.s <- cap.m[-c(6:7),] ## males captured in summer
cap.m.s15 <- cap.m.s[c(1:5),] ## 2015 only
cap.m.s16 <- cap.m.s[-c(1:5),] ## 2016 only
cap.m.w <- cap.m[c(6:7),] ## males captured in winter
  
mean(cap.f.s$kg)
  (sd(cap.f.s$kg)) / (sqrt(nrow(cap.f.s)))
mean(cap.m.s$kg)
  (sd(cap.m.s$kg)) / (sqrt(nrow(cap.m.s)))  
mean(cap.m.s15$kg)
  (sd(cap.m.s15$kg)) / sqrt(nrow(cap.m.s15))
mean(cap.m.s16$kg)
  (sd(cap.m.s16$kg)) / sqrt(nrow(cap.m.s16))
t.test(cap.f.s$kg, cap.m.s15$kg, paired = FALSE) # no sig diff between weights at capture (summer)

mean(cap.f.w$kg)
  (sd(cap.f.w$kg)) / (sqrt(nrow(cap.f.w)))
mean(cap.m.w$kg)
  (sd(cap.m.w$kg)) / (sqrt(nrow(cap.m.w)))
t.test(cap.f.w$kg, cap.m.w$kg, paired = FALSE) # no sig diff between weights at capture (winter)

## weight at capture (all together)
mean(cap$kg)
(sd(cap$kg)) / (sqrt(nrow(cap)))

##########################################################
### (2) linear model to find informed seasonal cutoff
##########################################################

## doesn't really make sense as a LM now that I think about it;
## maybe random effect or ANOVA.
## how could I look at weight change to find a seasonal switch?

lm1 <- lm(kg ~ date, data=porc.wts)
summary(lm1)

plot(kg ~ date, data=porc.wts)
abline(lm1)

##########################################################
### (3) t-test for significance between seasonal weights
##########################################################

##### FROM HERE FORWARD, NEED TO SELECT A CUTOFF DATE

## separate summer/winter weights
s.wts <- porc.wts[porc.wts$date < as.Date('2015-11-01') | porc.wts$date > as.Date('2016-02-29'),]
w.wts <- porc.wts[porc.wts$date > as.Date('2015-10-31') & porc.wts$date < as.Date('2016-03-01'),]

## if multiple summer weights, need to take the average
ids <- unique(s.wts$id)
s.avgs <- NULL

for(i in ids){
  s.avg.i <- mean(s.wts$kg[s.wts$id == i])
  sex.i <- s.wts$sex[s.wts$id == i]
  wt.id <- data.frame(i, s.avg.i, sex.i[1]) ## it will make rows for each; only need 1
  colnames(wt.id) <- c("id", "s_avg", "sex")
  s.avgs <- rbind(s.avgs, wt.id)
}

## if multiple winter weights, need to take the average
ids <- unique(w.wts$id)
w.avgs <- NULL

for(i in ids){
    w.avg.i <- mean(w.wts$kg[w.wts$id == i])
    wt.id <- data.frame(i, w.avg.i)
    colnames(wt.id) <- c("id", "w_avg")
    w.avgs <- rbind(w.avgs, wt.id)
}

## now combine them (need dplyr for 'left_join')
avg.wts <- left_join(s.avgs, w.avgs) ## will get a warning; see next step
avg.wts$id <- as.factor(avg.wts$id) 

## now get rid of rows with NA (it means the animal didn't have both summer & winter weights)
avg.wts <- avg.wts[!is.na(avg.wts$w_avg),]
avg.wts <- avg.wts[,c("id", "sex", "s_avg", "w_avg")] # reorder
avg.wts

## now do statistics!
t.test(avg.wts$s_avg, avg.wts$w_avg, paired=TRUE)

## females
avg.wts.f <- avg.wts[avg.wts$sex == 'F',]
t.test(avg.wts.f$s_avg, avg.wts.f$w_avg, paired = TRUE)

## males
avg.wts.m <- avg.wts[avg.wts$sex == 'M',]
t.test(avg.wts.m$s_avg, avg.wts.m$w_avg, paired = TRUE)

## means:
(mean_s <- mean(avg.wts$s_avg)) # 8.06 kg
(mean_w <- mean(avg.wts$w_avg)) # 7.41 kg
means <- cbind(mean_s, mean_w)

## barplot:
par(mfrow=c(1,1))
barplot(means, ylab="Body mass (kg)")
write.csv(avg.wts, 'csvs/avg.wts.091016.csv')

## how much did they lose?
(1 - (mean_w / mean_s)) * 100 # overall: 8.22%
(1 - (mean(avg.wts.f$w_avg)) / mean(avg.wts.f$s_avg)) * 100 # females: 6.39%

(1 - (mean(avg.wts.m$w_avg)) / mean(avg.wts.m$s_avg)) * 100 # males: 9.99%

avg.wts$pct_dec <- (1-(avg.wts$w_avg / avg.wts$s_avg))*100

## what about max wt difference for each porcupine (not just between avg summer/winter)
ids <- unique(porc.wts$id)
wt.diff <- NULL

for (i in ids){
      porc.wts.i <- porc.wts[porc.wts$id == i,]    
      max_diff <- max(porc.wts.i$kg) - min(porc.wts.i$kg)  
      diff_pct <- (1-(min(porc.wts.i$kg)/max(porc.wts.i$kg)))*100
      month_max <- porc.wts.i$month[which.max(porc.wts.i$kg)]
      month_min <- porc.wts.i$month[which.min(porc.wts.i$kg)]
      wt.diff.i <- data.frame(i, month_max, month_min, max_diff, diff_pct)
      wt.diff <- rbind(wt.diff, wt.diff.i)
}

wt.diff <- wt.diff[wt.diff$max_diff != 0,] ## get rid of animals with only 1 weight
wt.diff <- wt.diff[-c(16:17),] ## and those whose weights were in the same season
wt.diff$sex <- c('f', 'f', 'm', 'm', 'f', 'f', 'f', 'f', 'm', 'f', 'f', 'm', 'm', 'f', 'm')

## the greatest percent loss was 37.98% for 15.11, between June and January (when he died)
## then 28.70 for 15.14 (also June-January), 19.46% for 15.03 (June-March), 16.50% for 15.12
## (July-December), 12.33 for 15.01 (May-January), etc.

mean(wt.diff$diff_pct[wt.diff$sex == 'f']) ## mean weight diff % for females
  sd(wt.diff$diff_pct[wt.diff$sex == 'f']) / length(wt.diff$diff_pct[wt.diff$sex == 'f'])
mean(wt.diff$diff_pct[wt.diff$sex == 'm']) ## mean weight diff % for males
  sd(wt.diff$diff_pct[wt.diff$sex == 'm']) / length(wt.diff$diff_pct[wt.diff$sex == 'm'])

## females reached their max weights in May, Oct, June, June, June, October, July, July, June
## females reached their min weights in Jan, May, Jan, Oct, Oct, June, Dec, Jan, Feb

## males reached their max weights in June, June, June, June, April, April
## males reached their min weights in March, Nov, Jan, Jan, Feb