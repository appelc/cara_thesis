#######################################
## Trying resource utilization functions (RUFs)
## with porcupine data
#######################################

## 1. First, load data
## 2. Then, need to extract the UD from "adehabitatHR" package
##    (set bandwidth, grid, and extent)
## 3. Then, create a table with id, coord, and UD height
##    a. Get UD height at each pixel
##    b. Get UD height at occurrence points only
## 4. Assign values of covariates (veg class, canopy height) to cells
##    a. For UD height at each pixel
##    b. For UD height at occurrece points only
## 5. Run RUF using package "ruf"
##    a. For UD height at each pixel
##    b. For UD height at occurrence points only

### FROM TIM:
# 1. Convert heights to reasonable numbers
# 2. Re-run Henrietta on 95%
# 2a. Run Henrietta wth Veg2
# 3. See if you can run on other porcupines
# 4. Theta??



library(adehabitatHR)
library(googlesheets)
library(raster)
library(rgdal)
library(ruf)

######################
## 1. First, load data
######################
gs_ls()
locs <- gs_title("Porc relocation data")
porc.locs <- data.frame(gs_read(ss=locs, ws="Relocations", is.na(TRUE), range=cell_cols(1:16)))
colnames(porc.locs) <- c("date", "id", "sess", "type", "time", "az", "utm_e", "utm_n", 
                             "obs", "loc", "pos", "notes", "xvar", "yvar", "cov", "error")
porc.locs <- subset(porc.locs, type %in% c("V","V*","P","P*","L"))
porc.locs$date <- as.Date(porc.locs$date, "%m/%d/%Y")
porc.locs$utm_e <- as.numeric(porc.locs$utm_e)
porc.locs$utm_n <- as.numeric(porc.locs$utm_n)

## Keep only animals with >= 5 locations
n <- table(porc.locs$id)
porc.locs <- subset(porc.locs, id %in% names(n[n >= 5]), drop=TRUE)
hr.data <- droplevels(porc.locs)

## Turn this into a Spatial Points Data Frame
porc.sp <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                                  data=data.frame(porc.locs$id),
                                  proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

## Load veg data
veg <- readOGR(dsn=".", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(porc.sp)

min.veg <- c(399298.8, 4635961.2)
max.veg <- c(401920.1, 4644228.8)

veg.x <- c(399298.8, 401920.1)
veg.y <- c(4635961.2, 4644228.8)

px.veg <- SpatialPoints(data.frame(veg.x, veg.y))
px.veg <- as(px.veg, "SpatialPixels") #specify cellsize, cells.dim 

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

## To figure out grid: (code from Ian)

# Cellsize for KDE estimates in meters
c = 10 ## this is the cell size I want (sq. meters?)

fake.kern <- kernelUD(xy=porc.sp, extent = 1, same4all = TRUE)
spdf <- raster(as(fake.kern[[1]],"SpatialPixelsDataFrame"))
eas <- diff(range(spdf@extent[1:2])) ## pulls out x min & max
nor <- diff(range(spdf@extent[3:4])) ## pulls out y min & max
if(eas > nor){
  g <- (eas/c)
} else {
  g <- (nor/c)
}

## Calculate KUD for all animals
kern.all <- kernelUD(xy=porc.sp, h = 60, grid = g, extent = 1, same4all = TRUE)
image(kern.all)

## Clip the KUD to the extent of the veg layer
kern.all <- kern.all[veg,]

######################
## 3. Then, create a table with id, coord, and UD height for each porc
## (too big to do them all in one .csv)
## a. Get UD height at each pixel
######################

ids <- names(all.kud)

ruf.list <- NULL

for(i in ids){
      ud.height.i <- kern.all[[i]]@data$ud
      coords.i <- kern.all[[i]]@coords
      ht.i <- cbind((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.i) <- c("id", "height", "x", "y")
      my.ruf.data.list[[i]] <- data.frame(ht.i) 
      ruf.list[[i]] <- data.frame(ht.i) 
}

######################
## b. Get UD height at occurrence points only
######################

## first, convert estUD to raster
## this combines step 2 (kernelUD)

## I used "writeOGR" and exported shapefiles because I was having trouble with raster(estUDm2spixdf()),
## But also used raster for "extract." Try GeoTIFF from writeOGR? This is kind of a mess...
## Can I just create separate raster objects for each animal instead of exporting them w/OGR?

ids <- unique(porc.locs$id)
ud_heights <- NULL

for(i in ids){
    i.locs <- subset(porc.locs, id == i)
    i.locs <- droplevels(i.locs)
    i.sp <- SpatialPointsDataFrame(data.frame(i.locs$utm_e, i.locs$utm_n),
                               data=data.frame(i.locs$id),
                               proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
    i.kud <- kernelUD(i.sp, h=60, extent=1, grid=1000)
    i.raster <- raster(estUDm2spixdf(i.kud))
    i.sp$udheight <- extract(i.raster, i.sp)
    name = paste(i, "raster")
    filepath <- file.path("C:","Users","Cara","Documents","cara_thesis", "rasters",
                        paste(i, "raster", ".shp", sep = ""))
    writeOGR(i.sp, dsn=filepath, layer=name, driver="ESRI Shapefile") ## not a raster!
    export <- data.frame(i.sp@data$i.locs.id, i.sp@data$udheight, i.sp@coords)
    ud_heights <- rbind(ud_heights, export)
    }

head(ud_heights)
colnames(ud_heights) <- c("id", "ud_height", "x", "y")

plot(ud_heights$y ~ ud_heights$x)
plot(i.raster) ## can I get raster objects for all of them separately?

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
##    a. For UD height at each pixel
######################

veg <- readOGR(dsn=".", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(porc.sp)

## do spatial join using package "sp"
## multiply "height" by 100 before converting to spdf **is 100 enough?**

## load one UD as an example
hen.ruf <- ruf.list[[1]]
hen.ruf$height2 <- hen.ruf$height*10000000000
hen.sp <- SpatialPointsDataFrame(data.frame(hen.ruf$x, hen.ruf$y),
                                  data=data.frame(hen.ruf$id, hen.ruf$height2),
                                  proj4string=CRS(proj4string(veg)))

## assign veg class to each cell (row)
hen.sp@data$veg <- over(hen.sp, veg)$Class_2

######################
## b. For UD height at occurrence points only
######################

head(ud_heights)

## the ud_heights object has points for ALL animals, with an "id" column

## multiply "height" by 100 before converting to spdf **is 100 enough?**
ud_heights$height2 <- (ud_heights$ud_height)*100
ud.sp <- SpatialPointsDataFrame(data.frame(ud_heights$x, ud_heights$y),
                                    data=data.frame(ud_heights$id, ud_heights$height2),
                                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

## assign veg class to each cell (row)
ud.sp@data$veg <- over(ud.sp, veg)$Class
head(ud.sp@data)

######################
## 5. Run RUF using package "ruf"
##    a. For UD height at each pixel
######################

## continue with spdf created above from ruf.list[[]]
## "ruf.fit" doesn't actually need a spdf...
## make a data.frame with only ud height, covariates, x, y

hen.df <- data.frame(hen.sp@data$hen.ruf.height2, hen.sp@coords, hen.sp$veg)
colnames(hen.df) <- c("ud", "x", "y", "veg")

## will "!is.na" fix the "subscript out of bounds" problem? 
## I Will need to "clip" the extent at some point, because a bunch of points within the UD have
## no covariate values! (Like ones out in the ocean or outside the area that I digitized for the
## veg polygons.) For now, just remove "NA" values for veg.

hen.df <- henr.df[!is.na(henr.df$veg),]
#tst <- subset(hen.df, ud >= (quantile(hen.df$ud, 0.05))) # there are more than 5% zeroes!
hen.df <- subset(hen.df, ud>0)

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)

## Estimate (unstandardized) coefficients
hen.fit <- ruf.fit(ud ~ factor(veg),
                    space = ~ x + y,
                    data=hen.df, theta=hval,
                    name="15.01",
                    standardized=F)

summary(hen.fit)

names(hen.fit)

######################
## 5. Run RUF using package "ruf"
######################

all.df <- data.frame(ud.sp$ud_heights.height2, ud.sp@coords, ud.sp$veg)
colnames(all.df) <- c("ud", "x", "y", "veg")
head(all.df)

all.df <- all.df[!is.na(all.df$veg),]

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)

## Estimate (unstandardized) coefficients
all.fit <- ruf.fit(ud ~ factor(veg),
                     space = ~ x + y,
                     data=all.df, theta=hval,
                     name="all porcupines",
                     standardized=F)

summary(all.fit)


## Estimate (standardized) coefficients
all.fit.2 <- ruf.fit(ud ~ factor(veg),
                    space = ~ x + y,
                    data=all.df, theta=hval,
                    name="all porcupines, standardized",
                    standardized=T)

summary(all.fit.2)

names(all.fit.2)

