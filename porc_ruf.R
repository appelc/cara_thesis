#######################################
## Trying resource utilization functions (RUFs)
## with porcupine data
#######################################

## 1. First, load data
## 2. Then, need to extract the UD from "adehabitatHR" package
## 3. Then, create a table with id, coord, and UD height
##    a. Get UD height at each pixel
##    b. Get UD height at occurrence points only
## 4. Assign values of covariates (veg class, canopy height) to cells
##    a. For UD height at each pixel
##    b. For UD height at occurrece points only
## 5. Run RUF using package "ruf"

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
all.kud <- kernelUD(porc.sp, h="href")
image(all.kud) #href doesn't look great

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

for(i in ids){
      ud.height.i <- all.kud[[i]]@data$ud
      coords.i <- all.kud[[i]]@coords
      ht.i <- cbind((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.i) <- c("id", "height", "x", "y")
     # porc_uds <- rbind(porc_uds, ht.i)
      mypath <- file.path("C:","Users","Cara","Documents", "cara_thesis", "RUFs",
                           paste(i, "_ud", ".csv", sep = ""))
      write.csv(ht.i, file=mypath)
}

######################
## b. Get UD height at occurrence points only
######################

## first, convert estUD to raster
## need to do separately for each animal (I get an error in "raster" function that it only works
## when all animals have the same grid... can add "same4all=TRUE", but is this what we want?)
## NOt as easy to just do "raster" with one animal, because using something like "kud.all[1]" makes
## it not an "estUDm" object anymore (estUD-multiple animals) but an "estUD" (single animal), and
## the "estUDm2spixdf" will only convert an estUDm object to a SPDF, not an estUD object... phew.

## this should be a loop, but I can't figure out how to save the rasters separately
## could do writeOGR instead of raster function?
## otherwise, have to do this each time?

#ids <- names(all.kud)
i <- ids[17] # just change number here each time
#for(i in ids){
i.locs <- subset(porc.locs, id == i)
i.locs <- droplevels(i.locs)
i.sp <- SpatialPointsDataFrame(data.frame(i.locs$utm_e, i.locs$utm_n),
                               data=data.frame(i.locs$id),
                               proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
i.kud <- kernelUD(i.sp, h=60, extent=1, grid=1000)
i.raster <- raster(estUDm2spixdf(i.kud))
#}

## and rename it here each time
raster18 <- i.raster
plot(raster18)

# now, you can extract just at the coordinates from the raster
porc.sp$udheight <- extract(all.kud.raster, porc.sp)
test.udheight <- extract(raster01, porc.sp)

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
##    a. For UD height at each pixel
######################

veg <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(porc.sp)

## do spatial join using package "sp"

## load one ud.csv as an example
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
stevie.sp@data$veg <- over(stevie.sp, veg)$Class
head(stevie.sp@data)

roze.sp@data$veg <- over(roze.sp, veg)$Class
head(roze.sp@data)

bowie.sp@data$veg <- over(bowie.sp, veg)$Class
head(bowie.sp@data)

######################
## 5. Run RUF using package "ruf"
######################

## continue with spdf created above from csv.ud
## "ruf.fit" doesn't actually need a spdf...
## make a data.frame with only ud height, covariates, x, y

stevie.df <- data.frame(stevie.sp$stevie.height2, stevie.sp@coords, stevie.sp$veg)
colnames(stevie.df) <- c("ud", "x", "y", "veg")
head(stevie.df)

roze.df <- data.frame(roze.sp$roze.height2, roze.sp@coords, roze.sp$veg)
colnames(roze.df) <- c("ud", "x", "y", "veg")
head(roze.df)

bowie.df <- data.frame(bowie.sp$bowie.height2, bowie.sp@coords, bowie.sp$veg)
colnames(bowie.df) <- c("ud", "x", "y", "veg")
head(bowie.df)

## will this fix the "subscript out of bounds" problem? no...
## I Will need to "clip" the extent at some point, because a bunch of points within the UD have
## no covariate values! (Like ones out in the ocean or outside the area that I digitized for the
## veg polygons.) For now, just remove "NA" values for veg.

stevie.df <- stevie.df[!is.na(stevie.df$veg),]
roze.df <- roze.df[!is.na(roze.df$veg),]
bowie.df <- bowie.df[!is.na(bowie.df$veg),]

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)

## Estimate (unstandardized) coefficients
bowie.fit <- ruf.fit(ud ~ factor(veg),
                    space = ~ x + y,
                    data=bowie.df, theta=hval,
                    name="15.05",
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

