## Try weighted compositional analysis
## ala Millspaugh et al. 2006

#assumes you've already run "porc_wca_summer.R" and have locations / veg loaded

######################
## 1. First, load porcupine location data & veg data
######################

head(porc.locs) # this should have all VHF and GPS points (1 per day)

## Separate out winter locations
win.locs <- porc.locs[(porc.locs$date > "2015-10-31") & (porc.locs$date < "2016-03-01"), ]

## Keep only animals with >= 5 locations
n <- table(win.locs$id)
win.locs <- subset(win.locs, id %in% names(n[n >= 5]), drop=TRUE)
win.locs <- droplevels(win.locs)

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

ids <- unique(win.locs$id)

ud.list <- list()
ud.winter.list <- list()
ud.clipped.list <- list()
contour.list <- list()
contour.winter.list <- list()
kde.areas <- list()

for (i in ids){
      locs.i <- porc.locs[porc.locs$id == i,]
      locs.i$id_season <- rep(paste(i, "_all", sep = ""), nrow(locs.i))
      locs.win.i <- win.locs[win.locs$id == i,]
      locs.win.i$id_season <- rep(paste(i, "_win", sep = ""), nrow(locs.win.i))
      locs.all.i <- rbind(locs.i, locs.win.i)
      sp.i <- SpatialPointsDataFrame(data.frame(locs.all.i$utm_e, locs.all.i$utm_n),
                                     data=data.frame(locs.all.i$id_season),
                                     proj4string=CRS(proj4string(veg)))
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
      # calculate UD on both IDs ("all" and "winter") with same4all = TRUE
      kern.i <- kernelUD(xy = sp.i, h = 60, grid = g, extent = 1, same4all = TRUE)
      kde.i <- kernel.area(kern.i, percent = c(50, 90, 95, 99), unin = "m", unout = "km2", standardize = FALSE)
      data.frame(kde.i, row.names = c("50", "90", "95", "99"))
      kde.areas[[i]] <- kde.i
  
      # make 99% contours (full and winter)
      cont99.all.i <- getverticeshr.estUD(kern.i[[1]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)
      cont99.win.i <- getverticeshr.estUD(kern.i[[2]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)
  
      # clip winter UD to 99% contour from ALL points (not just winter), and veg extent
      win.ud.i <- (kern.i[[2]])[cont99.all.i,]
      win.ud.i <- win.ud.i[veg,]
  
      # save full UD, winter UD, and clipped UD:
      ud.list[[i]] <- kern.i[[1]]
      ud.winter.list[[i]] <- kern.i[[2]]
      ud.clipped.list[[i]] <- win.ud.i ##it's now a "SpatialPixelsDataFrame"
  
      # and save the contours:
      contour.list[[i]] <- cont99.all.i
      contour.winter.list[[i]] <- cont99.win.i # don't actually use this for anything later
}

## it's cool to look at a few here:
image(ud.clipped.list$`16.15`)
plot(veg, add=TRUE)
plot(contour.list$`16.15`, add=TRUE, border="blue", lwd=2)
plot(contour.winter.list$`16.15`, add=TRUE, border="green", lwd=2)

## output KDE areas
#write.csv(kde.areas, "csvs/kde_areas_winter_050916.csv")

######################
## 3. Then, create a list of tables with id, coord, and UD height for each porc
##    a. For UD height at each pixel
######################

ids <- unique(win.locs$id)
height.list.winter <- list()

for(i in ids){
      ud.i <- ud.clipped.list[[i]]
      ud.height.i <- ud.i$ud
      coords.i <- ud.i@coords
      ht.coords.i <- data.frame((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.coords.i) <- c("id", "height", "x", "y")
      height.list.winter[[i]] <- data.frame(ht.coords.i) 
}

## wireframe plots! better function to get lat/lon or put it on a map?
wireframe(height ~ x * y, data=height.list.winter$`16.15`, drape=TRUE, main="16.15 winter UD height")

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
######################

ids <- unique(win.locs$id)
final.list.winter <- list()
for (i in ids){
        ht.i <- height.list.winter[[i]]
        spdf.i <- SpatialPointsDataFrame(data.frame(ht.i$x, ht.i$y),
                                   data=data.frame(ht.i$id, ht.i$height),
                                   proj4string = CRS(proj4string(veg)))
        spdf.i@data$veg <- over(spdf.i, veg)$Class_2
        df.i <- data.frame(i, spdf.i@data$ht.i.height, spdf.i@coords, spdf.i@data$veg)
        colnames(df.i) <- c("id", "ud", "x", "y", "veg")
        df.i <- df.i[!is.na(df.i$veg),]
        final.list.winter[[i]] <- df.i
}

######################
## 5. For weighted compositional analysis: 
## sum raw UD values by veg type and divide the summed UD values by the 
## total UD value of all patches to obtain a UD-weighted estimate of use 
## for each habitat type for each individual animal 
## - (Millspaugh et al. 2004 p. 391)
######################

## need to reclassify "0" use values because log(0) = -Inf, which means that veg category
## will be excluded. Eads et al. 2012 use 0.30, "the minimum value that reduced type I error
## rates in simulation studies (see Bingham et al. 2007)." I'll use 1e-10...
## ** do it for the sum or as the raw use? check out that paper and do error sensitivity **
## example: for 15.02, if I change raw UD=0 to 1e-10, the log_ud_weights are meadow = -12, 
## pasture = -12, shrub swale = -13, and wooded swale = -13.
## but when I change ud sum (after "aggregate") from 0 to 1e-10, the weights are all -18.4
## I'll do it the second way for now (change summed UD for each category to 1e-10 if it's 0)

ids <- unique(win.locs$id)
tables.winter <- list()
for (i in ids){
        ud.i <- final.list.winter[[i]]
        table.i <- aggregate(ud ~ veg, data=ud.i, FUN = sum)
        table.i$ud[table.i$ud == 0] <- 1e-10 #do this here or before "aggregate"?
        table.i$ud_weight <- table.i$ud / sum(table.i$ud)
        table.i$log_ud_weight <- log(table.i$ud_weight)
        tables.winter[[i]] <- table.i
}

## now, need to calcluate "log-transformed availability data"

### CLIP POLYGONS Stolen from: https://philmikejones.wordpress.com/2015/09/01/clipping-polygons-in-r
## because "intersect" excludes some edge polygons and "gIntersection" doesn't retain polygon ids

ids <- unique(win.locs$id)
veg.99kdes <- list()
veg.areas <- list()
for (i in ids){
        cont99.i <- contour.list[[i]]
        clip.i <- gIntersection(cont99.i, veg, byid = T) #this is just a SpatialPolygons (no data)
        row.names(clip.i) <- gsub("homerange ", "", row.names(clip.i))    
        keep <- row.names(clip.i)
        clip.i <- spChFIDs(clip.i, keep) #changes feature IDs in the SP
        clip.data <- as.data.frame(veg@data[keep,]) #this is what we'll add as @data to the SPDF
        clip.spdf <- SpatialPolygonsDataFrame(clip.i, clip.data)  #this is fixed!
        clip.spdf <- clip.spdf[!is.na(clip.spdf@data$Class_2),] #get rid of NAs
        area.all <- gArea(clip.spdf, byid = TRUE) #units should be m^2
        veg.df.i <- data.frame(clip.spdf$Class_2, area.all)
        colnames(veg.df.i) <- c("veg", "area")
        veg.areas.i <- aggregate(area ~ veg, data=veg.df.i, FUN = sum)
        veg.areas.i$prop_area <- veg.areas.i$area / sum(veg.areas.i$area)
        veg.areas.i$log_avail <- log(veg.areas.i$prop_area)
        veg.99kdes[[i]] <- clip.spdf
        veg.areas[[i]] <- veg.areas.i
}

## good, no self-intersection errors!
## any missing polygons at edges?
plot(veg.99kdes$`16.17`) #all look great!

######################
## 6. Subtract differences in log-transformed availability data from the
##    log-transformed use data for each animal and then test for overall
##    selection using Wilks' lambda
######################

## combine used and avail in the same table
## we have "tables" (a list) and "veg.areas" (a list)

ids <- unique(win.locs$id)
full.table.winter <- NULL
final.table.winter <- NULL
for (i in ids){
      tables.i <- tables.winter[[i]]
      veg.areas.i <- veg.areas[[i]]       
      tables.i$log_avail <- veg.areas.i$log_avail
      tables.i$id <- rep(i, nrow(tables.i))
      tables.i$sel <- tables.i$log_ud_weight - tables.i$log_avail
      full.table.winter <- rbind(full.table.winter, tables.i)
      final.df <- data.frame(tables.i$id, tables.i$veg, tables.i$log_ud_weight,
                             tables.i$log_avail, tables.i$sel)
      colnames(final.df) <- c("id", "veg", "log_ud_wt", "log_avail", "sel")
      final.table.winter <- rbind(final.table.winter, final.df)
}

write.csv(final.table.winter, "csvs/wt_comp_analysis_winter_w-gps_050916.csv")

## box plot:
par(mar=c(5, 9, 3, 3), xpd=FALSE)
boxplot(sel ~ veg, data = final.table.winter, xaxt="n", las=2, ylim=c(-17, 5), horizontal=TRUE, outline=T)
axis(1)
title(xlab = "Differences in log ratio", line=3)
abline(v=0, lty = 1, col="red")

######################
## 7. Compute Wilk's lambda to test for overall selection
######################

## function "manova" or "Wilks.test"? the latter is in package "rrcov"
## "groups" are veg, and we are testing for differences between "log_ud_wt" and "log_avail"

groups <- as.factor(final.table.winter[,2])
x <- as.matrix(final.table.winter[,3:4])

## can do method "c" for mean and variance or "rank" for wilks lambda ranks
wilks.winter <- Wilks.test(x, grouping = groups, method="rank")
wilks.winter

######################
## 8. If use differes significantly from availability (p-value for Wilks lambda):
##    calculate the mean and st. dev. for the log-ratio differences,
##    and use these to rank each habitat type
## - Then, use t-test to assess difference between ranks and to determine where 
##   selection differed by habitat pairs (Millspaugh et al. 2006, p. 392)
######################

## calculate mean "sel" for each habitat type
## we have several categories with "-Inf" values for some animals
## (indicating that the UD height was 0, b/c log(0) is -Inf)... how to deal with these?
## Don't want to take them out because the use IS 0... change them all to 0?

veg_types <- unique(final.table.winter$veg)
means_table_win <- NULL
for (j in veg_types) {
        veg.j <- final.table.winter[final.table.winter$veg == j,]
        veg.j$log_ud_wt[veg.j$log_ud_wt == "-Inf"] <- 0 ## change "-Inf" values to 0s
        mean.sel.j <- mean(veg.j$sel)
        sd.sel.j <- sd(veg.j$sel)
        se.sel.j <- (sd.sel.j)/sqrt(nrow(veg.j))
        table.j <- data.frame(j, mean.sel.j, sd.sel.j, se.sel.j)
        colnames(table.j) <- c("veg", "mean", "sdev", "se")  
        means_table_win <- rbind(means_table_win, table.j)
}

## rank veg types and plot means:
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = means_table_win$mean + means_table_win$se,
              ymin = means_table_win$mean - means_table_win$se)

p <- ggplot(data = means_table_win, aes(x = reorder(veg, mean), y = mean))
p + geom_bar(stat = "identity", position = dodge, fill = factor(veg_names_col$veg_colors)) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) + ylab("Mean log-ratio differences (+/- SE)")

## where did the legend go?

## assign colors to each veg class for plotting
veg_names <- c("Beach", "Beachgrass dune", "Pasture", "Conifer forest", "Brackish marsh", 
               "Coastal scrub", "Meadow", "Freshwater marsh", "Wooded swale", "Shrub swale",
               "Fruit tree")
veg_colors <- c("khaki", "khaki3", "darkolivegreen3", "darkolivegreen4", "aquamarine", "khaki4",
                "yellow3", "cadetblue1", "aquamarine4", "darkseagreen3", "coral1")
veg_names_col <- data.frame(veg_names, veg_colors)

## t-tests to assess difference between ranks:



######################
######################
## cool figures but this is really messy (fix later):
sum.sp <- SpatialPointsDataFrame((data.frame(sum.locs$utm_e, sum.locs$utm_n)),
                                 data = data.frame(sum.locs$id),
                                 proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))


win.sp <- SpatialPointsDataFrame((data.frame(win.locs$utm_e, win.locs$utm_n)),
                                  data = data.frame(win.locs$id), 
                                  proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
ids <- unique(win.locs$id)
for (i in ids){
      veg.i <- veg.99kdes[[i]]
      mypath <- file.path("figures", "kdes_with_veg", "winter", "with_gps", paste(i, "_veg_99kde", ".png", sep = ""))
      png(file=mypath)
      mytitle = paste("99% KDE ", i, sep = "")
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(veg.i, main = mytitle)
      leg.txt <- sort(unique(veg$Class_2))
      leg.col <- c("khaki1", "khaki3", "aquamarine", "khaki4", "darkolivegreen4", "cadetblue1",
                   "coral1", "yellow3", "darkolivegreen3", "darkseagreen3", "aquamarine4")
      legend("topright", inset=c(-0.3,0), legend = leg.txt, pch = 15, col = leg.col, cex=0.9)
      plot(veg.i[veg.i$Class_2 == "Beach",], add=TRUE, col="khaki1")
      plot(veg.i[veg.i$Class_2 == "Beachgrass dune",], add=TRUE, col="khaki3")
      plot(veg.i[veg.i$Class_2 == "Brackish marsh",], add=TRUE, col="aquamarine")
      plot(veg.i[veg.i$Class_2 == "Coastal scrub",], add=TRUE, col="khaki4")
      plot(veg.i[veg.i$Class_2 == "Conifer forest",], add=TRUE, col="darkolivegreen4")
      plot(veg.i[veg.i$Class_2 == "Freshwater marsh",], add=TRUE, col="cadetblue1")
      plot(veg.i[veg.i$Class_2 == "Fruit tree",], add=TRUE, col="coral1")
      plot(veg.i[veg.i$Class_2 == "Meadow",], add=TRUE, col="yellow3")
      plot(veg.i[veg.i$Class_2 == "Pasture",], add=TRUE, col="darkolivegreen3")
      plot(veg.i[veg.i$Class_2 == "Shrub swale",], add=TRUE, col="darkseagreen3")
      plot(veg.i[veg.i$Class_2 == "Wooded swale",], add=TRUE, col="aquamarine4")
      plot(sum.sp[sum.sp$sum.locs.id == i,], add=TRUE, pch=16, cex=1, col="black")
      plot(win.sp[win.sp$win.locs.id == i,], add=TRUE, pch=16, cex=1, col="red")
      legend("bottomright", inset=c(-0.3,0), legend = c("Summer points", "Winter points"), pch=16, col=c("black", "red"), cex=0.9)
      dev.off() 
}

## may need to run this again to be able to plot again:
#dev.off()