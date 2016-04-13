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

means
barplot(means$mean.i, names.arg=means$i, ylab="area (m^2)", xlab="veg class", las=2)

## try without pasture b/c it's huge
means2 <- means[-10,]
barplot(means2$mean.i, names.arg=means2$i, ylab="area (m^2)", xlab="veg class", las=2)

## plot a line corresponding to a "patch" of a given bandwidth
## (assume patch is a circle with radius equal to bandwidth, in meters)
radius <- 60          # set bandwidth here
area <- pi*(radius^2)
abline(h=area)

## what's the avg without pasture?
mean.all.2 <- mean(means2$mean.i)
mean.all.2

## and without pasture, beach strand, brackish marsh, and beachgrass dune
means3 <- means[-c(5, 9, 10, 14),]
means3
barplot(means3$mean.i, names.arg=means3$i, ylab="area (m^2)", xlab="veg class", las=2)
abline(h=area)

## what's avg size without these?
mean.all.3 <- mean(means3$mean.i)
mean.all.3

## what radius (h) does this correspond to?
sqrt((mean.all.3)/pi)

## surprisingly close to 60 meters!!!
