################
## Home ranges
################

library(adehabitatHR)
library(googlesheets)
library(raster)

######################
## 1. First, load porcupine location data
######################
gs_ls()
locs <- gs_title("Porc relocation data")
porc.vhf <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE), range=cell_cols(1:8)))
colnames(porc.vhf) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n")
porc.vhf <- subset(porc.vhf, type %in% c("V","V*","P","P*","L"))
porc.vhf$utm_e <- as.numeric(porc.vhf$utm_e)
porc.vhf$utm_n <- as.numeric(porc.vhf$utm_n)
## check date class before running one of the next lines
porc.vhf$date <- as.Date(porc.vhf$date, "%m/%d/%Y") 
#porc.vhf$date <- as.Date(porc.vhf$date, origin = as.Date("1899-12-30"))

## Incorporate GPS data (1 random point per day)
gps.pts <- read.csv("daily.gps.csv")
gps.pts$type <- rep("gps", nrow(gps.pts))
gps.pts$az <- rep(NA, nrow(gps.pts))
gps.pts.df <- data.frame(gps.pts$Date, gps.pts$Animal.ID, gps.pts$Session, gps.pts$type,
                         gps.pts$Time, gps.pts$az, gps.pts$UTM.E, gps.pts$UTM.N)
colnames(gps.pts.df) <- colnames(porc.vhf)
gps.pts.df$date <- as.Date(gps.pts.df$date, "%Y-%m-%d")
gps.pts.df$id <- as.factor(gps.pts.df$id)

## combine VHF with GPS points
porc.locs <- rbind(porc.vhf, gps.pts.df)

## exclude 15.04
porc.locs <- porc.locs[porc.locs$id != '15.04',]

## subset summer locations (before Nov 1 or after March 1) and winter locations
sum.locs <- porc.locs[(porc.locs$date < "2015-11-01") | (porc.locs$date > "2016-02-29"), ]
win.locs <- porc.locs[(porc.locs$date > "2015-10-31") & (porc.locs$date < "2016-03-01"), ]

## Keep only animals with >= 5 locations
n <- table(porc.locs$id)
porc.locs <- subset(porc.locs, id %in% names(n[n >= 5]), drop=TRUE)
porc.locs <- droplevels(porc.locs)

n <- table(sum.locs$id)
sum.locs <- subset(sum.locs, id %in% names(n[n >= 5]), drop=TRUE)
sum.locs <- droplevels(sum.locs)

n <- table(win.locs$id)
win.locs <- subset(win.locs, id %in% names(n[n >= 5]), drop=TRUE)
win.locs <- droplevels(win.locs)

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

## Calculate grid & extent based on desired cell size (# meters on each side)
## For for each animal separately 

ids <- unique(porc.locs$id)
kde.areas <- list()
kud.all <- list()
contours <- list()

for (i in ids){
        locs.i <- porc.locs[porc.locs$id == i,]
        locs.i$id_season <- rep(paste(i, "_all", sep = ""), nrow(locs.i))
        locs.sum.i <- sum.locs[sum.locs$id == i,]
        locs.sum.i$id_season <- rep(paste(i, "_sum", sep = ""), nrow(locs.sum.i))
        locs.win.i <- win.locs[win.locs$id == i,]
        locs.win.i$id_season <- rep(paste(i, "_win", sep = ""), nrow(locs.win.i))
        locs.all.i <- rbind(locs.i, locs.sum.i, locs.win.i)
        sp.i <- SpatialPointsDataFrame(data.frame(locs.all.i$utm_e, locs.all.i$utm_n),
                                       data=data.frame(locs.all.i$id_season),
                                       proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
        c = 10   ## desired cell size (meters)
        fake.kern <- kernelUD(xy = sp.i, extent = 1)
        spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
        eas <- diff(range(spdf@extent[1:2]))
        nor <- diff(range(spdf@extent[3:4]))
        if(eas > nor){
          g <- (eas/c)
        } else {
          g <- (nor/c)
        }
        # calculate UD on all IDs ("all", "winter", and "summer") with same4all = TRUE
        kern.i <- kernelUD(xy = sp.i, h = 60, grid = g, extent = 1, same4all = TRUE)
        kde.i <- kernel.area(kern.i, percent = c(50, 90, 95, 99), unin = "m", unout = "km2", standardize = FALSE)
        data.frame(kde.i, row.names = c("50", "90", "95", "99"))
        kde.areas[[i]] <- kde.i
        kud.all[[i]] <- kern.i
        
        # make contours
        # make 95% contours (full and summer)
        cont95 <- list()
        for (j in 1:length(kern.i)){
            cont95.i <- getverticeshr.estUD(kern.i[[j]], percent = 95, unin = "m", unout = "km2", standardize = FALSE)
            cont95[[j]] <- cont95.i
        }
        contours[[i]] <- cont95
}

## output KDE areas
write.csv(kde.areas, "csvs/kde_areas_w-gps_061416.csv")

## manipulate list
ids <- names(kde.areas)
kde.all <- NULL
for(i in ids){
    kdes.i <- NULL
  for (j in 1:length(kde.areas[[i]])){
        kde <- (kde.areas[[i]])[j]
        kde$season <- rep(j, nrow(kde))
        kde$method <- rep('kde', nrow(kde))
        colnames(kde) <- c('km2', 'season', 'method')
        per <- c('50', '90', '95', '99')
        kde.df <- data.frame(i, per, kde$km2, kde$season, kde$method)
        colnames(kde.df) <- c('id', 'percent', 'area', 'season', 'method')
        kde.df$season[kde.df$season == 1] <- 'all'
        kde.df$season[kde.df$season == 2] <- 'summer'
        kde.df$season[kde.df$season == 3] <- 'winter'
        kdes.i <- rbind(kdes.i, kde.df)
        }
      kde.all <- rbind(kde.all, kdes.i)    
}

# quick plot:
f <- c('15.01', '15.02', '15.05', '15.07', '15.08', '15.09', '15.10', '15.12', '15.13', '16.16', '16.17')
m <- c('15.03', '15.04', '15.06', '15.11', '15.14', '16.15', '16.18')
plot(veg.ext)
for (i in f){
  plot(contours[[1:5]][[2]], add = TRUE, lty = 1, lwd = 2)
}

##################################
## Calculate MCP and export as CSV
##################################

## all points
all.df <- data.frame(porc.locs$utm_e, porc.locs$utm_n, porc.locs$id, porc.locs$date)
colnames(all.df) <- c("x", "y", "id", "date")
all.spdf <- SpatialPointsDataFrame(data.frame(all.df$x, all.df$y),
                                   data=data.frame(all.df$id),
                                   proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

mcp.all95 <- mcp(all.spdf, percent = 95, unin="m", unout="km2")
plot(mcp.all95, main="95% MCP, All")
mcp.all95.df <- data.frame(mcp.all95$id, rep('95', nrow(mcp.all95)), mcp.all95$area,
                           rep('all', nrow(mcp.all95)), rep('mcp', nrow(mcp.all95)))
colnames(mcp.all95.df) <- c('id', 'percent', 'area', 'season', 'method')
write.csv(mcp.all95.df, "csvs/mcp_all_95.csv")

## summer
sum.df <- data.frame(sum.locs$utm_e, sum.locs$utm_n, sum.locs$id, sum.locs$date)
colnames(sum.df) <- c("x", "y", "id", "date")
sum.spdf <- SpatialPointsDataFrame(data.frame(sum.df$x, sum.df$y),
                                   data = data.frame(sum.df$id),
                                   proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
mcp.sum95 <- mcp(sum.spdf, percent = 95, unin = 'm', unout = 'km2')
plot(mcp.sum95, main = '95% MCP, Summer')
mcp.sum95.df <- data.frame(mcp.sum95$id, rep('95', nrow(mcp.sum95)), mcp.sum95$area,
                           rep('summer', nrow(mcp.sum95)), rep('mcp', nrow(mcp.sum95)))
colnames(mcp.sum95.df) <- c('id', 'percent', 'area', 'season', 'method')
write.csv(mcp.sum95.df, "csvs/mcp_sum_95.csv")

## winter
win.df <- data.frame(win.locs$utm_e, win.locs$utm_n, win.locs$id, win.locs$date)
colnames(win.df) <- c("x", "y", "id", "date")
win.spdf <- SpatialPointsDataFrame(data.frame(win.df$x, win.df$y),
                                   data = data.frame(win.df$id),
                                   proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
mcp.win95 <- mcp(win.spdf, percent = 95, unin = 'm', unout = 'km2')
plot(mcp.win95, main = '95% MCP, Winter')
mcp.win95.df <- data.frame(mcp.win95$id, rep('95', nrow(mcp.win95)), mcp.win95$area,
                           rep('winter', nrow(mcp.win95)), rep('mcp', nrow(mcp.win95)))
colnames(mcp.win95.df) <- c('id', 'percent', 'area', 'season', 'method')
write.csv(mcp.win95.df, "csvs/mcp_win_95.csv")

mcp.all <- rbind(mcp.all95.df, mcp.sum95.df, mcp.win95.df)

##############################
### Now do analysis & try overlap metrics
##############################

hr.all <- rbind(kde.all, mcp.all) # combine MCP and KDE
hr.all$sex[hr.all$id %in% c('15.01', '15.02', '15.05', '15.07', '15.08', '15.09', '15.10', '15.12', '15.13', '16.16', '16.17')] <- 'f'
hr.all$sex[hr.all$id %in% c('15.03', '15.04', '15.06', '15.11', '15.14', '16.15', '16.18')] <- 'm'

## aggregate and compute SE
hr.summary <- aggregate(hr.all$area, by = list(sex = hr.all$sex, season = hr.all$season, method = hr.all$method, percent = hr.all$percent),
                        FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
hr.summary <- do.call(data.frame, hr.summary) # make data frames instead of matrices
colnames(hr.summary) <- c('sex', 'season', 'method', 'percent', 'mean', 'sd', 'n')
hr.summary$se <- hr.summary$sd / sqrt(hr.summary$n)
hr.summary$names <- c(paste(hr.summary$sex, hr.summary$season, sep = '_'))
write.csv(hr.summary, 'homeranges061416.csv')

## now do this all over again using ONLY animals with both summer and winter home ranges
## b/c these are the ones used in paired t-tests below

hr.all.mod <- subset(hr.all, id %in% c('15.01', '15.02', '15.03', '15.07', '15.11', '15.12', '15.13', '15.14', '16.15', '16.17', '16.18'), drop = TRUE)
hr.all.mod$id <- droplevels(hr.all.mod$id)
hr.mod.summary <- aggregate(hr.all.mod$area, by = list(sex = hr.all.mod$sex, season = hr.all.mod$season, method = hr.all.mod$method, percent = hr.all.mod$percent),
                        FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
hr.mod.summary <- do.call(data.frame, hr.mod.summary) # make data frames instead of matrices
colnames(hr.mod.summary) <- c('sex', 'season', 'method', 'percent', 'mean', 'sd', 'n')
hr.mod.summary$se <- hr.mod.summary$sd / sqrt(hr.mod.summary$n)
hr.mod.summary$names <- c(paste(hr.mod.summary$sex, hr.mod.summary$season, sep = '_'))
write.csv(hr.mod.summary, 'homeranges_mod061416.csv')

#################################################
## now make plots
#################################################

## first, subset only 95% KDEs
hr.summary.kde <- hr.summary[hr.summary$method == 'kde' & hr.summary$percent == '95',]
hr.summary.kde$sex <- as.character(hr.summary.kde$sex)
hr.summary.kde$sex[hr.summary.kde$sex == 'm'] <- 'males'
hr.summary.kde$sex[hr.summary.kde$sex == 'f'] <- 'females'

hr.mod.summary.kde <- hr.mod.summary[hr.mod.summary$method == 'kde' & hr.mod.summary$percent == '95',]
hr.mod.summary.kde$sex <- as.character(hr.mod.summary.kde$sex)
hr.mod.summary.kde$sex[hr.mod.summary.kde$sex == 'm'] <- 'males'
hr.mod.summary.kde$sex[hr.mod.summary.kde$sex == 'f'] <- 'females'

## now plot with all animals
#par(mar = c(5, 6, 4, 5) + 0.1)
tabbedMeans <- tapply(hr.summary.kde$mean, list(hr.summary.kde$season, hr.summary.kde$sex), function(x) c(x = x))
tabbedSE <- tapply(hr.summary.kde$se, list(hr.summary.kde$season, hr.summary.kde$sex), function(x) c(x = x))
barCenters <- barplot(height = tabbedMeans, beside = TRUE, las = 1, ylim = c(0, 0.6), cex.names = 1,
                      main = "Porcupine home ranges", ylab = "95% KDE area (sq. km)",
                      border = "black", axes = TRUE, legend.text = TRUE,
                      args.legend = list(title = "Season", x = "topright", cex = .9))
segments(barCenters, tabbedMeans - tabbedSE, barCenters, tabbedMeans + tabbedSE, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE, barCenters, tabbedMeans + tabbedSE, lwd = 1.5, angle = 90, code = 3, length = 0.05)

## now plot with only animals used in paired t-tests
tabbedMeans <- tapply(hr.mod.summary.kde$mean, list(hr.mod.summary.kde$season, hr.mod.summary.kde$sex), function(x) c(x = x))
tabbedSE <- tapply(hr.mod.summary.kde$se, list(hr.mod.summary.kde$season, hr.mod.summary.kde$sex), function(x) c(x = x))
barCenters <- barplot(height = tabbedMeans, beside = TRUE, las = 1, ylim = c(0, 0.6), cex.names = 1,
                      main = "Porcupine home ranges", ylab = "95% KDE area (sq. km)",
                      border = "black", axes = TRUE, legend.text = TRUE,
                      args.legend = list(title = "Season", x = "topright", cex = .9))
segments(barCenters, tabbedMeans - tabbedSE, barCenters, tabbedMeans + tabbedSE, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE, barCenters, tabbedMeans + tabbedSE, lwd = 1.5, angle = 90, code = 3, length = 0.05)

#################################################
## t-tests: use 95% KDE for now
#################################################

##### female summer vs. winter
s.f <- hr.all[hr.all$sex == 'f' & hr.all$season == 'summer' & hr.all$method == 'kde' & hr.all$percent == 95,]
w.f <- hr.all[hr.all$sex == 'f' & hr.all$season == 'winter' & hr.all$method == 'kde' & hr.all$percent == 95,]
s.f2 <- s.f[-c(3, 5, 6, 7),] # only keep those with both summer & winter for paired t-test

t.test(s.f2$area, w.f$area, paired = TRUE) # n = 6
t.test(s.f$area, w.f$area, paired = FALSE) # n = 10 / 6

## aggregate, compute standard error, plot, and add error bars
females <- rbind(w.f, s.f2)
f.summary <- aggregate(females$area, by = list(season = females$season),
                       FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
f.summary <- do.call(data.frame, f.summary) # make data frames instead of matrices
f.summary$se <- f.summary$x.sd / sqrt(f.summary$x.n)
colnames(f.summary) <- c("season", "mean", "sd", "n", "se")
barCenters <- barplot(height = f.summary$mean, names.arg = f.summary$names, beside = true, las = 2,
                      ylim = c(0, 0.6), cex.names = 0.75, xaxt = "n",
                      main = "Female home ranges",
                      ylab = "95% KDE area (sq. km)", border = "black", axes = TRUE)
text(x = barCenters, y = par("usr")[3] - 0.03, adj = 0.5, labels = f.summary$season, xpd = TRUE)
segments(barCenters, f.summary$mean - f.summary$se, barCenters, f.summary$mean + f.summary$se, lwd = 1.5)
arrows(barCenters, f.summary$mean - f.summary$se * 2, barCenters, f.summary$mean + f.summary$se * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)

##### male summer vs. winter
s.m <- hr.all[hr.all$sex == 'm' & hr.all$season == 'summer' & hr.all$method == 'kde' & hr.all$percent == 95,]
w.m <- hr.all[hr.all$sex == 'm' & hr.all$season == 'winter' & hr.all$method == 'kde' & hr.all$percent == 95,]
s.m2 <- s.m[-c(2),]

t.test(s.m2$area, w.m$area, paired = TRUE) # n = 5
t.test(s.m$area, w.m$area, paired = FALSE) # n = 5 / 6

## summer female vs. male
t.test(s.f$area, s.m$area, paired = FALSE) # n = 10 / 6

## winter female vs. male
t.test(w.f$area, w.m$area, paired = FALSE) # n = 6 / 5

## both summer vs. winter
s.both <- hr.all[hr.all$season == 'summer' & hr.all$method == 'kde' & hr.all$percent == 95,]
w.both <- hr.all[hr.all$season == 'winter' & hr.all$method == 'kde' & hr.all$percent == 95,]
s.both2 <- s.both[-c(4, 5, 7, 8, 9),]

t.test(s.both2$area, w.both$area, paired = TRUE) # n = 11
t.test(s.both$area, w.both$area, paired = FALSE) # n = 16 / 11


######### try a different way
head(hr.all)
hr.summary <- aggregate(hr.all$area, by = list(sex = hr.all$sex, season = hr.all$season),
                         FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
hr.summary <- do.call(data.frame, hr.summary) # make data frames instead of matrices

## compute standard error
hr.summary$se <- hr.summary$x.sd / sqrt(hr.summary$x.n)
colnames(hr.summary) <- c("sex", "season", "mean", "sd", "n", "se")
hr.summary$names <- c(paste(hr.summary$sex, hr.summary$season, sep = '_'))

## plot
#par(mar = c(5, 6, 4, 5) + 0.1)
plotTop <- max(hr.summary$mean) + hr.summary[hr.summary$mean == max(hr.summary$mean), 6] * 3
barCenters <- barplot(height = hr.summary$mean, names.arg = hr.summary$names, beside = true, las = 2,
                      ylim = c(0, plotTop), cex.names = 0.75, xaxt = "n",
                      main = "Porcupine home ranges",
                      ylab = "95% KDE area", border = "black", axes = TRUE)

# Specify the groupings. use srt = 45 for a 45 degree string rotation
text(x = barCenters, y = par("usr")[3] - 0.03, srt = 45, adj = 1, labels = hr.summary$names, xpd = TRUE)
segments(barCenters, hr.summary$mean - hr.summary$se, barCenters, hr.summary$mean + hr.summary$se, lwd = 1.5)
arrows(barCenters, hr.summary$mean - hr.summary$se * 2, barCenters, hr.summary$mean + hr.summary$se * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)

#############################

#### home range size vs. capture weight
#### using "porc.wts" from "season-weight.R" script

head(porc.wts)
sum.cap.wts <- porc.wts[c(1:3, 5:13, 15),] ## keep all summer captures except 15.04

head(hr.all)
sum.hr <- hr.all[hr.all$percent == 95 & hr.all$season == 'summer' & hr.all$method == 'kde',]
sum.hr <- sum.hr[1:13,] ## keep only summer captures

sum.hr$cap.wt <- sum.cap.wts$kg
sum.hr

m1 <- lm(area ~ cap.wt, data = sum.hr)
summary(m1)
plot(area ~ cap.wt, data = sum.hr)
abline(m1)

## overlap: **FINISH** (did I store all UDs?)
for (i in ids){
  porcs.overlap <- kerneloverlaphr(kud.all[[i]], method="HR", conditional=FALSE)
  View(porcs.overlap) ## Cool! try other methods
}

