## Try an RSF
## with simulated random "available" points
## and change summer/winter cutoff to model likelihood over time

library(spatstat)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(sp)
library(googlesheets)
library(ResourceSelection)
library(dplyr)

setwd("C:/Users/Cara/Documents/__RESEARCH/ANALYSIS/032816")

## load location points
gs_ls()
locs <- gs_title("Porc relocation data")
porc.locs.all <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE), range=cell_cols(1:16)))
colnames(porc.locs.all)
colnames(porc.locs.all) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n", 
                             "obs", "loc", "pos", "notes", "xvar", "yvar", "cov", "error")
porc.locs <- subset(porc.locs.all, type %in% c("V","V*","P","P*","L"))
porc.locs$date <- as.Date(porc.locs$date, "%m/%d/%Y")

## now make it spatial points data frame
locs <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                               data=data.frame(porc.locs),
                               proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

## load veg polygons shapefile
veg <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(locs) ## it's almost the same but not exactly; needs to be the same for "over" function
is(veg) # it's a spatial polygons data frame

## I collapsed veg classes into 5 categories 
## QUESTION: how to do "ifelse" to create a new column based on veg type, for multiple tests
head(veg)
unique(veg$Class_3)

# Design 1: all points for all animals are pooled (used), and 
# compared to the study area (available)

# Use my study area extent shapefile for now
# but could try the buffering again, based on greatest distance between points

study.area <- readOGR(dsn="D:/GIS DATA/Relocations", layer="Study_site_N", verbose=TRUE)
proj4string(study.area) <- proj4string(locs) ## it's almost the same but not exactly; needs to be the same for "over" function

# Plot all the data
plot(study.area) 
plot(veg, add=TRUE)
plot(locs, add=TRUE, col="red")
plot(locs[locs$id=="15.01",], col="blue", add=TRUE) ## can see diff individuals

# Create random points within the available study area
# use "length(locs)" to make sure we have the same no. of avail. pts as used pts
available.pts <- spsample(study.area, n=(length(locs)), type="random")
plot(available.pts, add=TRUE, col="green")

# Now add 0/1 for pres/avail
# and add columns for "season"

porc.pts <- as.data.frame(locs) # convert to data.frame
porc.pts$pres <- rep(1, nrow(porc.pts)) # add a 1 for all presence (used) pts
porc.pts$name_p <- rep("present", nrow(porc.pts)) # add a column for description
porc.pts$season <- rep(NA, nrow(porc.pts)) # add a column for season (will fill in later)
porc.pts$name_s <- rep(NA, nrow(porc.pts)) # add a column for season description
colnames(porc.pts)
porc.pts <- porc.pts[,c(1, 2, 17:22)] # remove unnecessary columns
colnames(porc.pts) <- c("date", "id", "x", "y", "pres", "name_p", "season", "name_s")
head(porc.pts)

# Do the same for the "available" points we created

available.pts <- as.data.frame(available.pts)
available.pts$pres <- rep(0, nrow(available.pts)) # add a 0 for all available pts
available.pts$name_p <- rep("available", nrow(available.pts))
available.pts$season <- rep(0, nrow(available.pts))
available.pts$name_s <- rep("summer", nrow(available.pts))
available.pts$id <- rep(NA, nrow(available.pts))
# set these all to June 1 so the for-loop will classify them as "summer" (0)
# available.pts$date <- rep(as.Date("2015-06-01"), nrow(available.pts))
# actually, set date as NA
available.pts$date <- rep(NA, nrow(available.pts))
# reorder columns b/c OCD
available.pts <- available.pts[,c("date","id", "x", "y", "pres", "name_p", "season", "name_s")]
# duplicate the available points to use as the winter set
avail2 <- available.pts 
# set these all to Dec 1 so the for-loop will classify them as "winter" (1)
# avail2$date <- rep(as.Date("2015-12-31"), nrow(avail2))
avail2$season <- rep(1, nrow(avail2))
avail2$name_s <- rep("winter", nrow(avail2))
# combine summer/winter available points
avail.pts <- rbind(available.pts, avail2)

# combine used and available pts
all.pts <- rbind(porc.pts, avail.pts)
nrow(all.pts) # good, it's the no. of locations * 3
head(all.pts)

# Convert the data back to a SpatialPointsDataFrame so we
# can extract info from our predictor layers
sp.pres.pts <- SpatialPointsDataFrame(data.frame(porc.pts$x, porc.pts$y), 
                                     porc.pts, proj4string=locs@proj4string)
sp.avail.pts <- SpatialPointsDataFrame(data.frame(avail.pts$x, avail.pts$y),
                                       avail.pts, proj4string=locs@proj4string)
## OR:
sp.all.pts <- SpatialPointsDataFrame(data.frame(all.pts$x, all.pts$y),
                                     all.pts, proj4string=locs@proj4string)

# Use the "over" function to get data from a polygon or 
# other shapefile / SpatialPolygonsDataFrame
sp.all.pts$veg <- over(sp.all.pts, veg)$Class_3
head(sp.all.pts) #yay!

sp.pres.pts$veg <- over(sp.pres.pts, veg)$Class_3
sp.avail.pts$veg <- over(sp.avail.pts, veg)$Class_3

# Now use a for-loop to iteratively change the season cutoff
# and incorporate the exponential RSF using package "Resource Selection"

# try using glm instead of rsf to see if it's faster
# and cutoff sequence is every 5 days

# first, do it for the whole study NOT subsampling for equal no. pts in sum/winter 
cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  sp.all.pts$season <- ifelse(sp.all.pts$date >= i, 1, 0)
  sp.all.pts$name_s <- ifelse(sp.all.pts$date >= i, "winter", "summer")
  sp.all.pts$season[(length(locs)+1):(2*length(locs))] <- 0
  sp.all.pts$season[(2*length(locs)+1):nrow(sp.all.pts)] <- 1
  sp.all.pts$season <- as.factor(sp.all.pts$season)
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=sp.all.pts)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

exp(cbind(OR = coef(glm.i), confint(glm.i)))
plot(ll ~ cutoff, main="all animals, no subsampling")
plot(aic ~ cutoff, main="all animals, no subsampling")

###########################################

## Now try a few individual porcupines

###########################################

# FOR-LOOP TO DO ALL PORCS
cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(j in sp.all.pts$id){
  for (i in cutoff){
    locs.j <- subset(sp.pres.pts, id == j)
    locs.j$season <- ifelse(locs.j$date >= i, 1, 0)
    locs.j$name_s <- ifelse(locs.j$date >= i, "winter", "summer")
    all.locs <- rbind(locs.j, sp.avail.pts)
    all.locs <- all.locs[!is.na(all.locs$veg),]
    glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=all.locs, na.exclude)
    aic <- c(aic, glm.i$aic)
    date <- as.Date.character(i, "%Y-%m-%d")
    ll <- c(ll, (logLik(glm.i)))
    odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
    odds <- bind_rows(odds, odds.rat)
  }
    mypath1 <- file.path("C:","Users","Cara","Documents","__RESEARCH", "_DATA", "Analysis", "AIC plots",
                        paste(j, "_AIC", ".png", sep = ""))
    png(file=mypath1)
        mytitle1 = paste("AIC_", j)
        plot(aic ~ cutoff, main = mytitle1)
        dev.off()
    mypath2 <- file.path("c:", "Users","Cara","Documents","__RESEARCH","_DATA","Analysis", "logLik plots",
                         paste(j, "_logLik", ".png", sep = ""))
    png(file=mypath2)
        mytitle2 = paste("logLik_", j)
        plot(ll ~ cutoff, main = mytitle2)
        dev.off()
}
#############################################
# for-loop take 2

ids <- unique(sp.pres.pts$id)

for(j in ids){
  cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)  
  aic <- NULL
  odds <- NULL
  ll <- NULL
  locs.j <- subset(sp.pres.pts, id == j)
      for (i in cutoff){
      locs.j$season <- ifelse(locs.j$date >= i, 1, 0)
      locs.j$name_s <- ifelse(locs.j$date >= i, "winter", "summer")
      all.locs <- rbind(locs.j, sp.avail.pts)
      all.locs <- all.locs[!is.na(all.locs$veg),]
      glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=all.locs)
      aic <- c(aic, glm.i$aic)
      date <- as.Date.character(i, "%Y-%m-%d")
      ll <- c(ll, (logLik(glm.i)))
      odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
      odds <- bind_rows(odds, odds.rat)
  }
  mypath1 <- file.path("C:","Users","Cara","Documents","__RESEARCH", "_DATA", "Analysis", "AIC plots",
                       paste(j, "_AIC", ".png", sep = ""))
  png(file=mypath1)
  mytitle1 = paste("AIC", j)
  plot(aic ~ cutoff, main = mytitle1)
  dev.off()
  mypath2 <- file.path("c:", "Users","Cara","Documents","__RESEARCH","_DATA","Analysis", "logLik plots",
                       paste(j, "_logLik", ".png", sep = ""))
  png(file=mypath2)
  mytitle2 = paste("logLik", j)
  plot(ll ~ cutoff, main = mytitle2)
  dev.off()
}

#############################################

# Henrietta
locs.01 <- subset(sp.all.pts, id %in% c("15.01", NA))
unique(locs.01$id)
head(locs.01)
View(locs.01)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.01$season <- ifelse(locs.01$date >= i, 1, 0)
  locs.01$name_s <- ifelse(locs.01$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.01)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

plot(ll ~ cutoff, main="15.01")
plot(aic ~ cutoff, main="15.01")
?plot

########
## anna
locs.02 <- subset(sp.all.pts, id %in% c("15.02", NA))
unique(locs.02$id)
head(locs.02)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.02$season <- ifelse(locs.02$date >= i, 1, 0)
  locs.02$name_s <- ifelse(locs.02$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.02)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

plot(ll ~ cutoff, main="15.02")
plot(aic ~ cutoff, main="15.02")
?plot

########
## rip
locs.03 <- subset(sp.all.pts, id %in% c("15.03", NA))
unique(locs.03$id)
head(locs.03)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.03$season <- ifelse(locs.03$date >= i, 1, 0)
  locs.03$name_s <- ifelse(locs.03$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.03)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

par(mfrow=c(1,1))
plot(ll ~ cutoff, main="15.03")
plot(aic ~ cutoff, main="15.03")

########
## stevie
locs.07 <- subset(sp.all.pts, id %in% c("15.07", NA))
unique(locs.07$id)
head(locs.07)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.07$season <- ifelse(locs.07$date >= i, 1, 0)
  locs.07$name_s <- ifelse(locs.07$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.07)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

plot(ll ~ cutoff, main="15.07")
plot(aic ~ cutoff, main="15.07")
?plot

########
## goldberry
locs.12 <- subset(sp.all.pts, id %in% c("15.12", NA))
unique(locs.12$id)
head(locs.12)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.12$season <- ifelse(locs.12$date >= i, 1, 0)
  locs.12$name_s <- ifelse(locs.12$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.12)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

plot(ll ~ cutoff, main="15.12")
plot(aic ~ cutoff, main="15.12")
?plot

########
## scully
locs.13 <- subset(sp.all.pts, id %in% c("15.13", NA))
unique(locs.13$id)
head(locs.13)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.13$season <- ifelse(locs.13$date >= i, 1, 0)
  locs.13$name_s <- ifelse(locs.13$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.13)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

plot(ll ~ cutoff, main="15.13")
plot(aic ~ cutoff, main="15.13")

########
## stormy
locs.14 <- subset(sp.all.pts, id %in% c("15.14", NA))
unique(locs.14$id)
head(locs.14)

cutoff <- seq(as.Date("2015-06-10"), as.Date("2016-03-28"), 5)
aic <- NULL
odds <- NULL
ll <- NULL

for(i in cutoff){
  locs.14$season <- ifelse(locs.14$date >= i, 1, 0)
  locs.14$name_s <- ifelse(locs.14$date >= i, "winter", "summer")
  glm.i <- glm(pres ~ veg*season, family="binomial"(link="logit"), data=locs.14)
  aic <- c(aic, glm.i$aic)
  date <- as.Date.character(i, "%Y-%m-%d")
  ll <- c(ll, (logLik(glm.i)))
  odds.rat <- cbind(i, data.frame(rbind(exp(coef(glm.i)))))
  odds <- bind_rows(odds, odds.rat) #instead of "rbind" b/c not all have "fruit" column
}

plot(ll ~ cutoff, main="15.14")
plot(aic ~ cutoff, main="15.14")

hist(porc.locs.all$date)

