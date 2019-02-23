## This is code I removed from porc_homerange.R for simplicity. It includes trying to outline
## the study area by combining all animal home ranges (or by maximum distance moved) as well as
## overlap metrices.

# ---------------------------------------------------------------
## 3. Try overall home range for all animals combined to outline study area
##    - should remove outliers: 15.04 (dispersal), 15.06 (exploratory)
# ---------------------------------------------------------------

locs.all <- porc.locs ## to create SPDF
locs.all$id2 <- rep('all', nrow(locs.all))
sp.all <- SpatialPointsDataFrame(data.frame(locs.all$utm_e, locs.all$utm_n),
                                 data=data.frame(locs.all$id2),
                                 proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

veg_buf <- gBuffer(veg, byid = FALSE, width = 200) ## buffer 200 m around study area
plot(sp.all)
plot(veg_buf, add = TRUE, border = 'red')

sp.all.clip <- sp.all[veg_buf,] ## clip points to buffered study area to remove outliers

c = 10   ## desired cell size (meters)
fake.kern <- kernelUD(xy = sp.all.clip, extent = 1)
spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
eas <- diff(range(spdf@extent[1:2]))
nor <- diff(range(spdf@extent[3:4]))
if(eas > nor){
  g <- (eas/c)
} else {
  g <- (nor/c)
}
# calculate UD on all IDs ("all", "winter", and "summer") with same4all = TRUE
kern.i <- kernelUD(xy = sp.all.clip, h = 60, grid = g, extent = 1, same4all = TRUE)
kde.i <- kernel.area(kern.i, percent = c(50, 90, 95, 99), unin = "m", unout = "km2", standardize = FALSE)
data.frame(kde.i, row.names = c("50", "90", "95", "99"))
kud.allporcs <- kern.i
kde.allporcs <- kde.i

# make contours
cont95.all <- getverticeshr.estUDm(kern.i, percent = 99, unin = "m", unout = "km2", standardize = FALSE)
plot(cont95.all) ## looks too restrictive. try 99%, increasing bandwidth, or using MCP.

# try MCP
mcp.all <- mcp(sp.all.clip, percent = 100, unin="m", unout="km2")
plot(mcp.all, main="95% MCP, All")
mcp.all95.df <- data.frame(mcp.all95$id, rep('95', nrow(mcp.all95)), mcp.all95$area,
                           rep('all', nrow(mcp.all95)), rep('mcp', nrow(mcp.all95)))

### Try by max distance moved (other than outliers)
ids <- unique(porc.locs$id)
max.dist <- NULL

for (i in ids){
  max <- 0
  cur.dist <- NULL
  locs.i <- porc.locs[porc.locs$id == i,]
  for(j in 1:nrow(locs.i)){
    cur.pt <- locs.i[j,]
    for(k in 1:nrow(locs.i)){
      compare.pt <- locs.i[k,]
      cur.dist <- sqrt((cur.pt$utm_e - compare.pt$utm_e)^2+(cur.pt$utm_n - compare.pt$utm_n)^2)
      if(cur.dist > max){
        max <- cur.dist
      }
    }
  }
  max.dist <- rbind(max.dist, data.frame(max, i))
}

max.dist ## for each animal
max.dist <- max.dist[-c(4, 6),] ## remove 15.04 and 15.06
max <- max(max.dist$max) ## overall (in meters)

# Now buffer all the points
study.area <- gBuffer(sp.all.clip, width = max)
plot(study.area)

# ---------------------------------------------------------------
## 5. Manipulate list of overlap metrices & do comparisons
# ---------------------------------------------------------------

ids <- names(overlap)
overlap.all <- NULL
for(i in ids){
  overlap.i <- overlap[[i]] ## only keep those with both sum & win for comparisons
  if (nrow(overlap.i) > 2){
    all <- stack(overlap.i[,1])
    all$s2 <- 'all'
    sum <- stack(overlap.i[,2])
    sum$s2 <- 'sum'
    win <- stack(overlap.i[,3])
    win$s2 <- 'win'
    overlap.df <- rbind(all, sum, win)
    overlap.df$id <- i
    overlap.df$s1 <- substr(overlap.df$ind, 7, 9)
  }
  overlap.all <- rbind(overlap.all, overlap.df[,c(4, 5, 3, 1)])
}

all.all <- overlap.all[overlap.all$s1 == 'all' & overlap.all$s2 == 'all',]
sum.all <- overlap.all[overlap.all$s1 == 'sum' & overlap.all$s2 == 'all',]
win.all <- overlap.all[overlap.all$s1 == 'win' & overlap.all$s2 == 'all',]
sum.win <- overlap.all[overlap.all$s1 == 'sum' & overlap.all$s2 == 'win',]

overlap.means <- c(mean(sum.all$values), mean(win.all$values), mean(sum.win$values))
overlap.se <- c(sd(sum.all$values)/sqrt(nrow(sum.all)), sd(win.all$values)/sqrt(nrow(win.all)), sd(sum.win$values)/sqrt(nrow(sum.win)))
overlap.n <- c(nrow(sum.all), nrow(win.all), nrow(sum.win))
comparisons <- c('sum_all', 'win_all', 'sum_win')
overlaps <- data.frame('comp' = comparisons, 'mean' = overlap.means, 'se' = overlap.se, 'n' = overlap.n)

t.test(sum.all$values, win.all$values, paired = TRUE, alternative = 'greater')
t.test(sum.win$values, sum.all$values, paired = TRUE, alternative = 'less')
t.test(sum.win$values, win.all$values, paired = TRUE, alternative = 'less')

## Calculate overlap between males and females during summer & winter
## need to re-run UDs on same grid in order to compare; 'kerneloverlap' includes this step, but input 'grid' parameter
sum.sp <- SpatialPointsDataFrame(data.frame(sum.locs$utm_e, sum.locs$utm_n),
                                 data = data.frame(sum.locs$id),
                                 proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
win.sp <- SpatialPointsDataFrame(data.frame(win.locs$utm_e, win.locs$utm_n),
                                 data = data.frame(win.locs$id),
                                 proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

overlap.sum <- kerneloverlap(sum.sp, method = 'UDOI', percent = 95, conditional = TRUE, 
                             h = 60, grid = 535.7288, extent = 1) ## 535.7288 is 'g' from above code using 'sum.locs'
overlap.win <- kerneloverlap(win.sp, method = 'UDOI', percent = 95, conditional = TRUE, 
                             h = 60, grid = 176.9492, extent = 1) ## 535.7288 is 'g' from above code using 'win.locs'

mean(overlap.sum)
sd(overlap.sum) / sqrt(nrow(overlap.sum))
nrow(overlap.sum)
mean(overlap.win)
sd(overlap.win) / sqrt(nrow(overlap.win))
nrow(overlap.win)

var.test(overlap.sum, overlap.win) # p<0.05; cannot assume equal variance
t.test(overlap.sum, overlap.win, var.equal = FALSE, paired = FALSE)

## can't do means on these because there are lots of redundancies (will inflate n)
sum.melt <- melt(overlap.sum, varnames = c('row', 'col'), na.rm = TRUE) #need library(reshape)
sum.unique <- unique(t(apply(sum.melt, 1, sort)))
overlap.s <- sum.unique[(sum.unique[,2] != sum.unique[,3]),] ## good; 136 pair combinations (17 animals)
overlap.s.nonzero <- data.frame(overlap.s[(overlap.s[,1] != 0),])

win.melt <- melt(overlap.win, varnames = c('row', 'col'), na.rm = TRUE)
win.unique <- unique(t(apply(win.melt, 1, sort)))
overlap.w <- win.unique[(win.unique[,2] != win.unique[,3]),] ## good; 55 pair combinations (11 animals)
overlap.w.nonzero <-   data.frame(overlap.w[(overlap.w[,1] != 0),])

overlap.means <- c(mean(overlap.s.nonzero[,1]), mean(overlap.w.nonzero[,1]))
overlap.se <- c(sd(overlap.s.nonzero[,1])/sqrt(20), sd(overlap.w.nonzero[,1])/sqrt(7)) ## or should sample sizes be 58 (# of pairs) and 7?
n <- c(17, 7)
comparisons <- c('summer', 'winter')
overlaps <- data.frame('comp' = comparisons, 'mean' = overlap.means, 'se' = overlap.se, 'n' = n)

## add m/f
sex <- c('f', 'f', 'm', 'm', 'f', 'm', 'f', 'f', 'f', 'f', 'm', 'f', 'f', 'm', 'm', 'f', 'm', 'm', 'm')
sex_key <- data.frame('id' = as.numeric(rownames(overlap.sum)), 'sex' = sex) ## as.numeric to match with below... otherwise get '15.1'... figure this out

overlap.s.nonzero$sex1 <- sex_key[match(overlap.s.nonzero$X2, sex_key$id), 'sex'] 
overlap.s.nonzero$sex2 <- sex_key[match(overlap.s.nonzero$X3, sex_key$id), 'sex'] 
overlap.w.nonzero$sex1 <- sex_key[match(overlap.w.nonzero$X2, sex_key$id), 'sex']
overlap.w.nonzero$sex2 <- sex_key[match(overlap.w.nonzero$X3, sex_key$id), 'sex']

ff <- overlap.s.nonzero[overlap.s.nonzero$sex1 == 'f' & overlap.s.nonzero$sex2 == 'f',]
mean.s.ff <- mean(ff$X1)
se.s.ff <- sd(ff$X1) / sqrt(nrow(ff))
mm <- overlap.s.nonzero[overlap.s.nonzero$sex1 == 'm' & overlap.s.nonzero$sex2 == 'm',]
mean.s.mm <- mean(mm$X1)
se.s.mm <- sd(mm$X1) / sqrt(nrow(mm))
fm <- overlap.s.nonzero[overlap.s.nonzero$sex1 != overlap.s.nonzero$sex2,]
mean.s.fm <- mean(fm$X1)
se.s.fm <- sd(fm$X1) / sqrt(nrow(fm))
overlap.means1 <- c(mean.s.ff, mean.s.mm, mean.s.fm)
overlap.se1 <- c(se.s.ff, se.s.mm, se.s.fm)
comparisons <- c('ff', 'mm', 'fm')
overlaps1 <- data.frame('comp' = comparisons, 'mean' = overlap.means1, 'se' = overlap.se1)

## Now look at overlap between core areas only (50% KDEs)
overlap50.sum <- kerneloverlap(sum.sp, method = 'UDOI', percent = 50, conditional = TRUE, 
                               h = 60, grid = 535.7288, extent = 1) ## 535.7288 is 'g' from above code using 'sum.locs'
overlap50.win <- kerneloverlap(win.sp, method = 'UDOI', percent = 50, conditional = TRUE, 
                               h = 60, grid = 176.9492, extent = 1) ## 535.7288 is 'g' from above code using 'win.locs'

mean(overlap50.sum)
sd(overlap50.sum) / sqrt(nrow(overlap50.sum))
nrow(overlap50.sum)
mean(overlap50.win)
sd(overlap50.win)/ sqrt(nrow(overlap50.win))
nrow(overlap50.win)

var.test(overlap50.sum, overlap50.win) # p<0.05; cannot assume equal variance
t.test(overlap50.sum, overlap50.win, var.equal = FALSE, paired = FALSE)


## can't do means on these because there are lots of redundancies (will inflate n)
sum.melt50 <- melt(overlap50.sum, varnames = c('row', 'col'), na.rm = TRUE)
sum.unique50 <- unique(t(apply(sum.melt50, 1, sort)))
overlap.s50 <- sum.unique50[(sum.unique50[,2] != sum.unique50[,3]),] ## good; 136 pair combinations (17 animals)
overlap.s.nonzero50 <- data.frame(overlap.s50[(overlap.s50[,1] != 0),])

win.melt50 <- melt(overlap50.win, varnames = c('row', 'col'), na.rm = TRUE)
win.unique50 <- unique(t(apply(win.melt50, 1, sort)))
overlap.w50 <- win.unique50[(win.unique50[,2] != win.unique50[,3]),] ## good; 55 pair combinations (11 animals)
overlap.w.nonzero50 <-   data.frame(overlap.w50[(overlap.w50[,1] != 0),])

mean(overlap.s.nonzero50$X1)
sd(overlap.s.nonzero50$X1) / sqrt(19)
mean(overlap.w.nonzero50$X1)
sd(overlap.w.nonzero50$X1)/ sqrt(11)

var.test(overlap.s.nonzero50$X1, overlap.w.nonzero50$X1) # p>0.05; assume equal variance
t.test(overlap.s.nonzero50$X1, overlap.w.nonzero50$X1, var.equal = TRUE, paired = FALSE)


overlap.means50 <- c(mean(overlap.s.nonzero50[,1]), mean(overlap.w.nonzero50[,1]))
overlap.se50 <- c(sd(overlap.s.nonzero50[,1])/sqrt(nrow(overlap.s.nonzero50)), sd(overlap.w.nonzero50[,1])/sqrt(nrow(overlap.w.nonzero50))) ## or should sample sizes be 17 and 7?
comparisons <- c('summer', 'winter')
overlaps50 <- data.frame('comp' = comparisons, 'mean' = overlap.means50, 'se' = overlap.se50)

## add m/f
sex <- c('f', 'f', 'm', 'm', 'f', 'm', 'f', 'f', 'f', 'f', 'm', 'f', 'f', 'm', 'm', 'f', 'm', 'm', 'm')
sex_key <- data.frame('id' = as.numeric(rownames(overlap.sum)), 'sex' = sex) ## as.numeric to match with below... otherwise get '15.1'... figure this out

overlap.s.nonzero50$sex1 <- sex_key[match(overlap.s.nonzero50$X2, sex_key$id), 'sex'] 
overlap.s.nonzero50$sex2 <- sex_key[match(overlap.s.nonzero50$X3, sex_key$id), 'sex'] 
overlap.w.nonzero50$sex1 <- sex_key[match(overlap.w.nonzero50$X2, sex_key$id), 'sex']
overlap.w.nonzero50$sex2 <- sex_key[match(overlap.w.nonzero50$X3, sex_key$id), 'sex']

mean.s.ff50 <- mean(overlap.s.nonzero50$X1[overlap.s.nonzero50$sex1 == 'f' & overlap.s.nonzero50$sex2 == 'f'])
mean.s.mm50 <- mean(overlap.s.nonzero50$X1[overlap.s.nonzero50$sex1 == 'm' & overlap.s.nonzero50$sex2 == 'm'])
mean.s.fm50 <- mean(overlap.s.nonzero50$X1[overlap.s.nonzero50$sex1 != overlap.s.nonzero50$sex2])

### ok, 'overlaps' tells us means. so more overlapped during winter than summer...?
## how many overlapped at all? what was the highest/lowest overlap between a pair?
## next do males / females
