##########################################################
### Exploring body mass change
### (1) create figure with trends for each animal
### (2) linear model to find informed seasonal cutoff
### (3) t-test for significance based on seasonal cutoff
### (4) report intervals between captures for each porcupine (exploratory)
### (5) split into capture periods instead of just summer/winter (exploratory)
##########################################################

library(googlesheets)
library(ggplot2)
library(dplyr)

## load from csv (only data used in thesis/paper)
porc.wts <- read.csv('c:/Users/Cara/Documents/__RESEARCH/porc_ms/data/porc_weights.csv')
porc.wts$Animal.ID <- as.factor(porc.wts$Animal.ID)
porc.wts$Date <- as.Date(porc.wts$Date, '%m/%d/%Y')
colnames(porc.wts) <- c('id', 'date', 'kg', 'sex', 'month')

## for manuscript figure, remove last 2 measurements (16.17 in October)
porc.wts <- porc.wts[-c(53:55),]


## OR:
## load from Google Drive (for current wts)
#gs_ls()
#wts <- gs_title("Porc weights & captures")
#porc.wts <- data.frame(gs_read(ss=wts, ws="Weights", is.na(TRUE), range=cell_cols(1:5)))
#colnames(porc.wts) <- c("id", "date", "kg", "sex", "month")
#head(porc.wts)
    ## format correctly
  #  porc.wts$id <- as.factor(porc.wts$id)
  #  porc.wts$date <- as.Date(porc.wts$date, "%m/%d/%Y")
  #  porc.wts$sex <- as.factor(porc.wts$sex)
  #  porc.wts$month <- as.factor(porc.wts$month)
  #  porc.wts$kg <- as.numeric(porc.wts$kg)

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
###       (a) female/male on same figure
##########################################################

#par(mar = c())
ggplot(data=porc.wts, aes(x = date, y = kg, group = id)) +
  scale_y_continuous(limits = c(4,11)) +   
  scale_x_date(limits = as.Date(c('2015-05-27','2016-09-10'))) + # for new data: change to '2016-07-01' and '2018-02-10'
  geom_vline(xintercept = as.numeric(as.Date('2015-11-01')), linetype='dotted', lwd = 1.5) + 
          geom_vline(xintercept = as.numeric(as.Date('2016-03-01')), linetype='dotted', lwd = 1.5) +      
          geom_point(data = f.wts, aes(y = f.wts$kg, colour = 'Female', shape = 'Female'), size = 4) + 
          geom_line(data = f.wts, colour = 'black', size = 1.5) +
          geom_point(data = m.wts, aes(y = m.wts$kg, colour = 'Male', shape = 'Male'), size = 4) + 
          geom_line(data = m.wts, colour = 'darkgray', size = 1.5, lty = 2) +
          ylab('Body mass (kg)') + #xlab('Date') + 
        theme(axis.text=element_text(size=20, colour = 'black'),
                axis.text.x=element_text(angle = 35, hjust = 1),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=22, color = 'black', margin = margin(0,15,0,0)),
                axis.line.x = element_line(size = 1, colour = 'black'),
                axis.line.y = element_line(size = 1, colour = 'black'),
                axis.ticks = element_line(size = 1, colour = 'black'),
                axis.ticks.length = unit(.3, 'cm'),
                panel.background = element_rect(fill = 'white')) +     
          scale_colour_manual("", 
                      breaks = c('Female', 'Male'),
                      values = c('black', 'darkgray'), guide = FALSE) +
          scale_shape_manual("",
                      breaks = c('Female', 'Male'),
                      values = c(16, 17), guide = FALSE)

## EXPORT AT 662 X 416

## I tried 'geom_rect' to shade the area between vertical lines (winter) but it didn't work

## For legend (add to 'theme': #legend.background = element_rect(colour = 'white'), legend.key = element_rect(fill = 'white'),
# legend.text = element_text(size = 18), legend.key.size = unit(1.5, "cm"), plot.margin=unit(c(0.8,0.25,0.75,0.5), 'cm')))

##########################################################  
###       (b) female/male on separate panels
##########################################################

porc.wts$sex <- as.character(porc.wts$sex)  
porc.wts$sex[porc.wts$sex == 'F'] <- 'Female'  #spell out for facet wrap labels
porc.wts$sex[porc.wts$sex == 'M'] <- 'Male'

## get rid of singletons for figure (04, 05, 16)?
porc.wts <- porc.wts[porc.wts$id != '15.04' & porc.wts$id != '15.05' & porc.wts$id != '16.16',]
  porc.wts <- droplevels(porc.wts)

ggplot(data = porc.wts, aes(x = date, y = kg, group = id)) +
    scale_y_continuous(limits = c(4,12)) +    # either (4,11) or (0,12.5) 
    scale_x_date(limits = as.Date(c('2015-05-27','2016-09-10'))) + 
    facet_wrap(~ sex, ncol = 1, scales = 'free_y') +
    geom_point(size = 3) + 
    geom_line(size = 1) +
      geom_vline(xintercept = as.numeric(as.Date('2015-11-01')), linetype='dashed', lwd = 1) + 
      geom_vline(xintercept = as.numeric(as.Date('2016-03-01')), linetype='dashed', lwd = 1) +  
    ylab('Body mass (kg)') + 
    annotate('segment', x = as.Date(-Inf, origin="1970-01-01"), xend = as.Date(Inf, origin="1970-01-01"), y = -Inf, yend = -Inf, size = 2, colour = 'black') + ## add x-axis to top panel
    theme(axis.text=element_text(size=20, colour = 'black'),
          axis.text.x=element_text(angle = 35, hjust = 1),
          axis.title.x=element_blank(), 
          axis.title.y=element_text(size=22, color = 'black', margin = margin(0,15,0,0)),
          axis.line.x = element_blank(),  #element_line(size = 1, colour = 'black'), ## otherwise 'annotate' will add an extra x-axis to bottom panel, making it look thicker
          axis.line.y = element_line(size = 1, colour = 'black'),
       #   axis.ticks = element_line(size = 1, colour = 'black'), # use default; 1 is too thick
          axis.ticks.length = unit(.3, 'cm'),
          panel.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 20))

## export at 550 x 700 TIFF (or 5.5 x 7 portrait for PDF)
## or 5.5 x 5.5 (old)
  
##########################################################
### (1b) descriptive stats
##########################################################

## weight at capture (summer): males vs. females
cap <- porc.wts[match(unique(porc.wts$id), porc.wts$id),] ## these are the 1st wts for each porc
cap.f <- cap[cap$sex == 'F',]
cap.f.s <- cap.f[-c(10:11),] ## females captured in summer (all 2015)
cap.f.w <- cap.f[c(10:11),] ## females captured in winter
cap.m <- cap[cap$sex == 'M',]
cap.m.s <- cap.m[-c(6:7),] ## males captured in summer
cap.m.s15 <- cap.m.s[c(1:5),] ## 2015 only
cap.m.s16 <- cap.m.s[-c(1:5),] ## 2016 only
cap.m.w <- cap.m[c(6:7),] ## males captured in winter
  
## means in summer
mean(cap.f.s$kg)
  (sd(cap.f.s$kg)) / (sqrt(nrow(cap.f.s))) #females 2015
mean(cap.m.s$kg)
  (sd(cap.m.s$kg)) / (sqrt(nrow(cap.m.s)))  #males (both)
mean(cap.m.s15$kg)
  (sd(cap.m.s15$kg)) / sqrt(nrow(cap.m.s15))  #males 2015
mean(cap.m.s16$kg)
  (sd(cap.m.s16$kg)) / sqrt(nrow(cap.m.s16))  #males 2016

## test for equal variance & t-test
var.test(cap.f.s$kg, cap.m.s15$kg)
  qf(0.95, 8, 4) ## find tabulated F using p=0.95, num df, denom df (if > computed F, assume homogenous variances)
t.test(cap.f.s$kg, cap.m.s15$kg, var.equal = TRUE, paired = FALSE) # no sig diff between weights at capture (summer), p=0.176
  qt(0.975, 12) ## p=0.975 b/c t-distribution is 2-tailed (greater than computed t, so confirm no difference between means)

## means in winter  
mean(cap.f.w$kg)
  (sd(cap.f.w$kg)) / (sqrt(nrow(cap.f.w)))  #females 
mean(cap.m.w$kg)
  (sd(cap.m.w$kg)) / (sqrt(nrow(cap.m.w)))  #males

## test for equal variance & t-test
var.test(cap.f.w$kg, cap.m.w$kg)
  qf(0.95, 1, 1) ## assume homogenous variance (p>0.05; tabulated F > computed F)
t.test(cap.f.w$kg, cap.m.w$kg, var.equal = TRUE, paired = FALSE) # no sig diff between weights at capture (winter)
  qt(0.975, 2)

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

##### FROM HERE FORWARD, NEED TO SELECT CUTOFF DATEs
## currently 2015-11-01 and 2016-02-29

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
### *** wrong method here; see bottom of script... ***
(1 - (mean_w / mean_s)) * 100 # overall: 8.22%
(1 - (mean(avg.wts.f$w_avg)) / mean(avg.wts.f$s_avg)) * 100 # females: 6.39%

(1 - (mean(avg.wts.m$w_avg)) / mean(avg.wts.m$s_avg)) * 100 # males: 9.99%

avg.wts$pct_dec <- (1-(avg.wts$w_avg / avg.wts$s_avg))*100

##############

## what about max wt difference for each porcupine (not just between avg summer/winter)

ids <- unique(porc.wts$id)
wt.diff <- NULL

for (i in ids){
      porc.wts.i <- porc.wts[porc.wts$id == i,]    
      max_diff <- max(porc.wts.i$kg) - min(porc.wts.i$kg)  
      diff_pct <- (1-(min(porc.wts.i$kg)/max(porc.wts.i$kg)))*100
      month_max <- porc.wts.i$month[which.max(porc.wts.i$kg)]
      month_min <- porc.wts.i$month[which.min(porc.wts.i$kg)]
      wt.diff.i <- data.frame(i, month_max, month_min, 'min_wt' = min(porc.wts.i$kg), 'max_wt' = max(porc.wts.i$kg), max_diff, diff_pct)
      wt.diff <- rbind(wt.diff, wt.diff.i)
}

wt.diff$sex <- c('f', 'f', 'm', 'm', 'f', 'm', 'f', 'f', 'f', 'f','m', 'f', 'f', 'm', 'm', 'f', 'f', 'm', 'm', 'm')
wt.diff <- wt.diff[wt.diff$max_diff != 0,] ## get rid of animals with only 1 weight
wt.diff <- wt.diff[-c(16:17),] ## and those whose weights were in the same season

## the greatest percent loss was 37.98% for m15.11, between June and January (when he died)
## then 28.70% for m15.14 (also June-January), 19.46% for m15.03 (June-March), 16.50% for f15.12
## (July-December), 12.33% for f15.01 (May-January), etc.

mean(wt.diff$diff_pct[wt.diff$sex == 'f']) ## mean weight diff % for females
  sd(wt.diff$diff_pct[wt.diff$sex == 'f']) / length(wt.diff$diff_pct[wt.diff$sex == 'f'])
mean(wt.diff$diff_pct[wt.diff$sex == 'm']) ## mean weight diff % for males
  sd(wt.diff$diff_pct[wt.diff$sex == 'm']) / length(wt.diff$diff_pct[wt.diff$sex == 'm'])

## females reached their max weights in May, June, June, June, June, July, July, Oct, Oct 
## females reached their min weights in Dec, Jan, Jan, Jan, Feb, May*, June^, Oct, Oct [*15.02; ^15.10 (not weighed after Oct)]

## males reached their max weights in April, April, June, June, June, June, Sept*, Sept* (*both 16.19 and 16.20 not weighed after early Sept; also potentially young and not parcitipating in the breeding season?)
## males reached their min weights in Nov, Jan, Jan, Jan, Feb, March, July*, July* (but both only compared to Sept; 16.19 and 16.20)
  
## PERCENT LOSS:

  # both sexes
  mean(wt.diff$max_diff / wt.diff$max_wt) * 100
    sd(wt.diff$max_diff / wt.diff$max_wt) / sqrt(nrow(wt.diff)) # SE

  # females  
  mean(wt.diff$max_diff[wt.diff$sex == 'f'] / wt.diff$max_wt[wt.diff$sex == 'f']) * 100
    sd(wt.diff$max_diff[wt.diff$sex == 'f'] / wt.diff$max_wt[wt.diff$sex == 'f']) / sqrt(nrow(wt.diff[wt.diff$sex == 'f',])) # SE
  
  # males
  mean(wt.diff$max_diff[wt.diff$sex == 'm'] / wt.diff$max_wt[wt.diff$sex == 'm']) * 100
    sd(wt.diff$max_diff[wt.diff$sex == 'm'] / wt.diff$max_wt[wt.diff$sex == 'm']) / sqrt(nrow(wt.diff[wt.diff$sex == 'm',]))

## t-tests on max/min instead of summer/winter
  t.test(wt.diff$min_wt, wt.diff$max_wt, paired = TRUE)
    
  # females
  t.test(wt.diff$min_wt[wt.diff$sex == 'f'], wt.diff$max_wt[wt.diff$sex == 'f'], paired = TRUE)
    
  # males
  t.test(wt.diff$min_wt[wt.diff$sex == 'm'], wt.diff$max_wt[wt.diff$sex == 'm'], paired = TRUE)

###
## should actually be max summer - min winter (yes kind of misleading bc of seasonal cutoffs, like 15.03 min is in March)

head(porc.wts)
  
  porc.wts$season[porc.wts$date < as.Date('2015-11-01') | porc.wts$date > as.Date('2016-02-29')] <- 's'
  porc.wts$season[porc.wts$date > as.Date('2015-10-31') & porc.wts$date < as.Date('2016-03-01')] <- 'w'
  porc.wts$season[porc.wts$date > as.Date('2016-10-31') & porc.wts$date < as.Date('2017-03-01')] <- 'w'

# remove animals whose min/max wts are in the same season (16.19, 16.20)
  
ids <- levels(porc.wts$id)
pct.decline <- NULL
  
  for (i in ids){
    porc.wts.i <- porc.wts[porc.wts$id == i,]    
    max_sum <- max(porc.wts.i$kg[porc.wts.i$season == 's'])
    min_win <- min(porc.wts.i$kg[porc.wts.i$season == 'w'])
    decline <- ((max_sum - min_win) / max_sum) * 100
    decline_df <- data.frame(i, max_sum, min_win, decline, 'sex' = porc.wts.i$sex[1])
    pct.decline <- rbind(pct.decline, decline_df)
  }

# remove porcs without wts in both seasons
pct.decline <- pct.decline[pct.decline$decline != '-Inf' & pct.decline$decline != 'NaN',] 

  mean(pct.decline$decline, na.rm = TRUE)
    sd(pct.decline$decline, na.rm = TRUE) / sqrt(nrow(pct.decline))
    # overall: mean 12.77%, SE 3.13
  
  mean(pct.decline$decline[pct.decline$sex == 'F'], na.rm = TRUE)
    sd(pct.decline$decline[pct.decline$sex == 'F'], na.rm = TRUE) / sqrt(nrow(pct.decline[pct.decline$sex == 'F',]))
    # overall: mean 8.32%, SE 2.19

  mean(pct.decline$decline[pct.decline$sex == 'M'], na.rm = TRUE)  
    sd(pct.decline$decline[pct.decline$sex == 'M'], na.rm = TRUE) / sqrt(nrow(pct.decline[pct.decline$sex == 'M',]))
    # overall: mean 17.22%, SE 5.52
    

##########################################################
### (4) report intervals between captures for each porcupine
##########################################################
    
head(porc.wts)
  
porc.wts$interval <- NA
porc.wts$mass_change <- NA

  for (y in 1:(nrow(porc.wts) - 1)){
    z <- y + 1
      if (porc.wts$id[y] == porc.wts$id[z]) {   ## so we're only comparing with same porc
        interval <- difftime(porc.wts$date[z], porc.wts$date[y], units = 'days')
          porc.wts$interval[z] <- interval
        mass_change <- porc.wts$kg[z] - porc.wts$kg[y]
          porc.wts$mass_change[z] <- mass_change
    }
  }

mean(porc.wts$interval, na.rm = TRUE)  ## mean interval between all captures (69 days)
min(porc.wts$interval, na.rm = TRUE)   ## min inerval (8 days)
max(porc.wts$interval, na.rm = TRUE)   ## max interval (210 days)

## but what about just between initial capture and first re-capture? (for summer captures)
## there's definitely a better way to do this...

## all summer captures (May 28 – July 23) were weighed again in October EXCEPT for: 
    ## 15.06 (weighed 11/06) and 15.07 (weighed 01/12)
## what about 15.12? captured 7/2, weighed 7/20 (18 days) and again 10/30 (102 days, or 120 days since capture)

first.wts <- porc.wts[c(2,6,9,16,20,22,24,26,31,35,38),]  ## 30 or 31 for 15.12? remove 18 for 15.07?

## of the 14 porcupines captured in summer 2015 (between 05-28 and 07-23), we re-weighed
## 11 of them between 10-09 and 11-06 (intervals ranging 99-149 days; mean interval 121)
## + 1 we didn't weigh until 01/12 (210 days)
## + 1 died in summer + 1 dispersed

## furthermore, of the 11 re-weighed in Oct-Nov, 9 had LOST mass and only 2 had gained mass (both F)
## avg mass difference was -0.51 kg

mean(first.wts$mass_change)

##########################################################
### (5) split into capture periods instead of just summer/winter
##########################################################

porc.wts$period[porc.wts$date <= as.Date('2015-07-31')] <- 'may_july_15'
porc.wts$period[porc.wts$date >= as.Date('2015-10-01') & porc.wts$date <= as.Date('2015-11-06')] <- 'oct_nov_15'
porc.wts$period[porc.wts$date >= as.Date('2016-01-01') & porc.wts$date <= as.Date('2016-02-28')] <- 'jan_feb_16'
porc.wts$period[porc.wts$date >= as.Date('2016-03-01') & porc.wts$date <= as.Date('2016-04-30')] <- 'mar_apr_16'
porc.wts$period[porc.wts$date >= as.Date('2016-06-01') & porc.wts$date <= as.Date('2016-07-31')] <- 'jun_jul_16'

porc.wts$season[porc.wts$date < as.Date('2015-11-01')] <- 'summer_15'
porc.wts$season[porc.wts$date >= as.Date('2015-11-01') & porc.wts$date < as.Date('2016-03-01')] <- 'winter_15_16'
porc.wts$season[porc.wts$date >= as.Date('2016-03-01') & porc.wts$date < as.Date('2016-11-01')] <- 'summer_16'

## if a porc has multiple weights in a single period, need to take the average during that period
wts_by_period <- aggregate(porc.wts$kg, list(porc.wts$id, porc.wts$period, porc.wts$sex), FUN = mean)
  colnames(wts_by_period) <- c('id', 'period', 'sex', 'kg')
avg_wts_period <- aggregate(wts_by_period$kg, list(wts_by_period$sex, wts_by_period$period), FUN = mean)
  avg_wts_period_length <- aggregate(wts_by_period$kg, list(wts_by_period$sex, wts_by_period$period), FUN = length)
  avg_wts_period$n <- as.numeric(avg_wts_period_length$x)
  colnames(avg_wts_period) <- c('sex', 'period', 'avg_kg', 'n')
  
## same by season (already did this above -- this is cleaner)
wts_by_season <- aggregate(porc.wts$kg, list(porc.wts$id, porc.wts$season, porc.wts$sex), FUN = mean)
  wts_by_season_length <- aggregate(porc.wts$kg, list(porc.wts$id, porc.wts$season, porc.wts$sex), FUN = length)
  wts_by_season$n <- as.numeric(wts_by_season_length$x)
  colnames(wts_by_season) <- c('id', 'season', 'sex', 'kg', 'n')  
  
avg_wts_season <- aggregate(wts_by_season$kg, list(wts_by_season$sex, wts_by_season$season), FUN = mean)
  avg_wts_season_length <- aggregate(wts_by_season$kg, list(wts_by_season$sex, wts_by_season$season), FUN = length)
  avg_wts_season$n <- as.numeric(avg_wts_season_length$x)
  colnames(avg_wts_season) <- c('sex', 'season', 'avg_kg', 'n')  

  
## t-tests (capture periods):
  
a <- wts_by_period[wts_by_period$period == 'may_july_15' & wts_by_period$sex == 'F',]
b <- wts_by_period[wts_by_period$period == 'may_july_15' & wts_by_period$sex == 'M',]
  t.test(a$kg, b$kg, paired = FALSE)
  # no difference between m & f in summer capture period

c <- wts_by_period[wts_by_period$period == 'may_july_15' & wts_by_period$sex == 'F',]
d <- wts_by_period[wts_by_period$period == 'oct_nov_15' & wts_by_period$sex == 'F',] 
  t.test(c$kg, d$kg, paired = FALSE)
    c <- c[-c(3,4),]
      t.test(c$kg, d$kg, paired = TRUE)
  # no difference between summer & fall for females (even when paired)
  
e <- wts_by_period[wts_by_period$period == 'may_july_15' & wts_by_period$sex == 'M',] 
f <- wts_by_period[wts_by_period$period == 'oct_nov_15' & wts_by_period$sex == 'M',] 
  t.test(e$kg, f$kg, paired = FALSE)       
    e <- e[-2,]
      t.test(e$kg, f$kg, paired = TRUE)
  # significant difference between summer & fall for males (when paired) but n = 4

g <- wts_by_period[wts_by_period$period == 'oct_nov_15' & wts_by_period$sex == 'F',]
h <- wts_by_period[wts_by_period$period == 'oct_nov_15' & wts_by_period$sex == 'M',]
  t.test(g$kg, h$kg, paired = FALSE)  
  # no difference between m & f in fall capture period
  
i <- wts_by_period[wts_by_period$period == 'may_july_15' & wts_by_period$sex == 'F',]
j <- wts_by_period[wts_by_period$period == 'jan_feb_16' & wts_by_period$sex == 'F',]
  t.test(i$kg, j$kg, paired = FALSE)
    i <- i[-c(3, 5:8),] 
    j <- j[-c(5,6),]
      t.test(i$kg, j$kg, paired = TRUE)
  # no difference between summer & winter for females (paired or unpaired)
      
k <- wts_by_period[wts_by_period$period == 'may_july_15' & wts_by_period$sex == 'M',]
l <- wts_by_period[wts_by_period$period == 'jan_feb_16' & wts_by_period$sex == 'M',]
  t.test(k$kg, l$kg, paired = FALSE)
    k <- k[-c(2:3),] 
    l <- l[-c(4,5),]
      t.test(k$kg, l$kg, paired = TRUE)
  # significanta difference between summer & winter for males (when paired) but n = 3
      
      
