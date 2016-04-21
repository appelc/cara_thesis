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

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

## Calculate KUD for all animals

## href just to compare
#all.kud <- kernelUD(porc.sp, h="href")
#image(all.kud) #href doesn't look great

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

kern.all <- kernelUD(xy=porc.sp, h = 60, grid = g, extent = 1, same4all = TRUE)
image(kern.all)

ba.result <- kerneloverlaphr(kern.all, meth="BA", conditional=TRUE)

## revisit extent parameter (actually calculate for cell ~3-6 meters)
all.kud <- kernelUD(porc.sp, h=60, extent=1, grid=1000)
image(all.kud)

######################
## 3. Then, create a table with id, coord, and UD height for each porc
## (too big to do them all in one .csv)
## a. Get UD height at each pixel
######################

ids <- names(all.kud)
#porc_uds <- NULL

my.ruf.data.list <- NULL

for(i in ids){
      ud.height.i <- all.kud[[i]]@data$ud
      coords.i <- all.kud[[i]]@coords
      ht.i <- cbind((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.i) <- c("id", "height", "x", "y")
      my.ruf.data.list[[i]] <- data.frame(ht.i)      
      # porc_uds <- rbind(porc_uds, ht.i)
      mypath <- file.path("C:","Users","Cara","Documents", "cara_thesis", "RUFs",
                           paste(i, "_ud", ".csv", sep = ""))
      write.csv(ht.i, file=mypath)
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

## load one ud.csv as an example
henrietta <- my.ruf.data.list[[1]]
henrietta$height2 <- henrietta$height*100
henr.sp <- SpatialPointsDataFrame(data.frame(henrietta$x, henrietta$y),
                                  data=data.frame(henrietta$id, henrietta$height2),
                                  proj4string=CRS(proj4string(veg)))

stevie <- read.csv("RUFs/15.07_ud.csv")
head(stevie)

roze <- read.csv("RUFs/15.08_ud.csv")
head(roze)

bowie <- read.csv("RUFs/15.05_ud.csv")
head(bowie)

## multiply "height" by 100 before converting to spdf **is 100 enough?**
stevie$height2 <- (stevie$height)*100
stevie.sp <- SpatialPointsDataFrame(data.frame(stevie$x, stevie$y),
                                    data=data.frame(stevie$id, stevie$height2),
                                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

roze$height2 <- (roze$height)*100
roze.sp <- SpatialPointsDataFrame(data.frame(roze$x, roze$y),
                                  data=data.frame(roze$id, roze$height2),
                                  proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

bowie$height2 <- (bowie$height)*100
bowie.sp <- SpatialPointsDataFrame(data.frame(bowie$x, bowie$y),
                                   data=data.frame(bowie$id, bowie$height2),
                                   proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

## assign veg class to each cell (row)
henr.sp@data$veg <- over(henr.sp, veg)$Class_3

stevie.sp@data$veg <- over(stevie.sp, veg)$Class
head(stevie.sp@data)

roze.sp@data$veg <- over(roze.sp, veg)$Class
head(roze.sp@data)

bowie.sp@data$veg <- over(bowie.sp, veg)$Class
head(bowie.sp@data)

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

## continue with spdf created above from csv.ud
## "ruf.fit" doesn't actually need a spdf...
## make a data.frame with only ud height, covariates, x, y

henr.df <- data.frame(henr.sp@data$henrietta.height2, henr.sp@coords, henr.sp$veg)
colnames(henr.df) <- c("ud", "x", "y", "veg")

stevie.df <- data.frame(stevie.sp$stevie.height2, stevie.sp@coords, stevie.sp$veg)
colnames(stevie.df) <- c("ud", "x", "y", "veg")
head(stevie.df)

roze.df <- data.frame(roze.sp$roze.height2, roze.sp@coords, roze.sp$veg)
colnames(roze.df) <- c("ud", "x", "y", "veg")
head(roze.df)

bowie.df <- data.frame(bowie.sp$bowie.height2, bowie.sp@coords, bowie.sp$veg)
colnames(bowie.df) <- c("ud", "x", "y", "veg")
head(bowie.df)

## will "!is.na" fix the "subscript out of bounds" problem? 
## I Will need to "clip" the extent at some point, because a bunch of points within the UD have
## no covariate values! (Like ones out in the ocean or outside the area that I digitized for the
## veg polygons.) For now, just remove "NA" values for veg.

henr.df <- henr.df[!is.na(henr.df$veg),]
stevie.df <- stevie.df[!is.na(stevie.df$veg),]
roze.df <- roze.df[!is.na(roze.df$veg),]
bowie.df <- bowie.df[!is.na(bowie.df$veg),]

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)

## Estimate (unstandardized) coefficients
henr.fit <- ruf.fit(ud ~ factor(veg),
                    space = ~ x + y,
                    data=henr.df, theta=hval,
                    name="15.01",
                    standardized=F)

summary(roze.fit)

## error in var(betas) + asycovbeta/con$nresamples : non-coformable arrays

# Estimate (standardized) coefficients
roze.fit <- ruf.fit(ud ~ factor(veg),
                    space = ~ x + y,
                    data=roze.df, theta=hval,
                    name="15.07 standardized",
                    standardized=T)
summary(stevie.fit)

names(stevie.fit)

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

