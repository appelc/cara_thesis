## Calculate areas of vegetation polygons
## to see if this can help inform a biologically meaningful kernel bandwidth
## to get at patch-level selection

library(rgdal)

## load veg polygons shapefile
veg <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg categories CA", verbose=TRUE)
veg@proj4string
#proj4string(veg) <- CRS("+proj=utm +zone=10 +datum=NAD83") ## close

## calculate means for each veg class
## can change between "Class", "Class_2", & "Class_3" columns in veg.class for less detail
## (also need to change in first line of loop: class.i)

veg.class <- unique(veg$Class)
means <- NULL

for(i in veg.class){
     class.i <- subset(veg, Class == i)
     class.i <- droplevels(class.i@data)
     mean.i <- mean(class.i$Area_1)
     x <- data.frame(i, mean.i)
     means <- rbind(means, x)
}

## with all veg types:
means
barplot(means$mean.i, names.arg=means$i, ylab="area (m^2)", xlab="veg class", las=2)

## plot a line corresponding to a "patch" of a given bandwidth
## (assume patch is a circle with radius equal to bandwidth, in meters)
radius <- 60          # set bandwidth here
area <- pi*(radius^2)
abline(h=area)

## fix plot margins
par(mai=c(1.7,1,0.3,0.2))

##############

## try without certain patches from the start

veg_2 <- subset(veg, Class %in% c("Freshwater marsh, non-wooded", "Wooded swale", "Coastal meadow",
                "Shrub swale", "Cultivated fruit", "Dune forest", "Exotics-dominated meadow",
                "Herbaceous swale", "Coastal scrub", "Dune woodland"))

## mean/median of all patches
mean_2 <- mean(veg_2@data$Area_1)
median(veg_2@data$Area_1)

## what radius (h) does this mean correspond to (assuming patch=circle)?
sqrt((mean_2)/pi)

## plot means for each veg class
veg.class <- unique(veg_2$Class)
means <- NULL

for(i in veg.class){
  class.i <- subset(veg_2, Class == i)
  class.i <- droplevels(class.i@data)
  mean.i <- mean(class.i$Area_1)
  x <- data.frame(i, mean.i)
  means <- rbind(means, x)
}

barplot(means$mean.i, names.arg=means$i, ylab="area (m^2)", las=2)

# plot lines for hypothetical patch (radius=60 m) AND mean patch size
abline(h=area)
abline(h=mean_2, lty=2)

