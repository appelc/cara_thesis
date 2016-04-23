#######################################
## Trying resource utilization functions (RUFs)
## with porcupine data
#######################################

## 1. First, load data
## 2. Then, extract the UD from "adehabitatHR" package (set bandwidth, grid, and extent)
## 3. Then, create a table with id, coord, and UD height
## 4. Assign values of covariates (veg class, canopy height) to cells
## 5. Run RUF using package "ruf"
## 
## 3a. - 5a.: For UD height at each pixel
## 3b. - 5b.: For UD height at occurrence points only

## *** Still need to clip grids to extent of study area / veg layer at some point ***
## *** and find paper to read about Matern correlation parameters (theta)

### FROM TIM:
# 1. Convert heights to reasonable numbers
## center & scale (x-mean)/sd ... but returns negative numbers for height=0
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
## 1. First, load porcupine location data & veg data
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
veg <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(porc.sp)

## Load veg extent boundary
#veg.dissolved <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg extent", verbose=TRUE)
#proj4string(veg.dissolved) <- proj4string(porc.sp)
#veg.ext <- as(veg.dissolved, "SpatialLines")

## Try study area b/c angles in veg extent too small
#study.area <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Area_extent", verbose=TRUE)
#proj4string(study.area) <- proj4string(porc.sp)
#study.ext <- as(study.area, "SpatialLines")

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

## Calculate grid & extent based on desired cell size
## Do for each animal separately 
## I tried using boundary=study.ext in kernelUD, but angles aren't acceptable
## (will need to clip grids later, I guess)

ids <- unique(porc.locs$id)
ud.list <- NULL

for (i in ids){
  i.locs <- subset(porc.locs, id == i)
  i.sp <- SpatialPointsDataFrame(data.frame(i.locs$utm_e, i.locs$utm_n),
                                    data=data.frame(i.locs$id),
                                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
  c = 5   ## desired cell size (meters)
  fake.kern <- kernelUD(xy = i.sp, extent = 1)
  spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
  eas <- diff(range(spdf@extent[1:2]))
  nor <- diff(range(spdf@extent[3:4]))
    if(eas > nor){
      g <- (eas/c)
    } else {
      g <- (nor/c)
    }
  kern.i <- kernelUD(xy = i.sp, h = 60, grid = g, extent = 1)
  ud.list[[i]] <- kern.i
}

image(ud.list[[11]]) ## check a few... looks good!

######################
## 3. Then, create a list of tables with id, coord, and UD height for each porc,
## and include a column for scaling the UD height: (x-mean)/sd
##    a. For UD height at each pixel
######################

ids <- unique(porc.locs$id)
height.list <- NULL

for(i in ids){
      ud.i <- ud.list[[i]]
      ud.height.i <- ud.i[[1]]@data$ud
      coords.i <- ud.i[[i]]@coords
      ht.coords.i <- data.frame((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.coords.i) <- c("id", "height", "x", "y")
      ht.coords.i$ht.scaled <- ((ht.coords.i$height) - (mean(ht.coords.i$height)))/(sd(ht.coords.i$height))
      height.list[[i]] <- data.frame(ht.coords.i) 
}

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
##    a. For UD height at each pixel
######################

## the "overlay" function takes a while, so first crop each grid
## subset based on top 95% (or 99%?) of UD values

## load one UD as an example
hen.ht <- height.list[[1]]
x <- quantile(hen.ht$ht.scaled, 0.05) ## or change to 0.01
y <- round(x[[1]], 14) ## next line only works if I round the quantile value
hen.ht95 <- hen.ht[hen.ht$ht.scaled >= y,]
hist(hen.ht95$ht.scaled) ## see the distribution

## now do them all!
ids <- unique(porc.locs$id)
ht95.list <- NULL

for (i in ids){
        ht.i <- height.list[[i]]
        x <- quantile(ht.i$ht.scaled, 0.05)
        y <- round(x[[1]], 14) ## next line only works if I round the quantile value
        ht95.i <- ht.i[ht.i$ht.scaled >= y,]
        ht95.list[[i]] <- data.frame(ht95.i)
}

## now assign veg class to each cell (row)
## (and eventually canopy height)

## this loop turns each data frame into a SPDF, does 'overlay' with veg class, 
## gets rid of cells where veg=NA (there shouldn't be many but they may mess up ruf.fit)
## then turns it BACK into a data frame for calculating RUF in next step

ids <- unique(porc.locs$id)
final.list <- NULL

for (i in ids){
        ht95.i <- ht95.list[[i]]  
        i.sp <- SpatialPointsDataFrame(data.frame(ht95.i$x, ht95.i$y),
                                       data=data.frame(ht95.i$id, ht95.i$ht.scaled),
                                       proj4string = CRS(proj4string(veg)))
        i.sp@data$veg <- over(i.sp, veg)$Class_2
        i.df <- data.frame(i, i.sp@data$ht95.i.ht.scaled, i.sp@coords, i.sp@data$veg)
        colnames(i.df) <- c("id", "ud", "x", "y", "veg")
        i.df <- i.df[!is.na(i.df$veg),]
        final.list[[i]] <- i.df
}

######################
## 5. Run RUF using package "ruf"
##    a. For UD height at each pixel
######################

## Need to switch from 64-bit R to 32-bit R now for the "ruf" package... 
## *hopefully* everything will be stored in the environment after restarting...

library(ruf)

## Now, "final.list" contains the data frames necessary to run ruf.fit
## (id, standardized ud height for top 95%, x, y, veg class)

hval <- c(0.2, 1.5)

ids <- unique(porc.locs$id)
ruf.list <- NULL
betas.list <- NULL

for (i in ids){
        i.df <- final.list[[i]]
        i.ruf <- ruf.fit(ud ~ factor(veg),
                         space = ~ x + y,
                         data = test.df, theta = hval,
                         name = "i",
                         standardized = F)
        ruf.list[[i]] <- i.ruf
        betas.list[[i]]  <- i.ruf$beta
}

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)


############################################################################3
## STEPS 3-5
## b. for height at occurrencep points only
############################################################################3

######################
## 3. Then, create a table with id, coord, and UD height for each porc
##    b. For UD height at occurrence points only
######################

## first, convert estUD to raster
## this combines step 2 (kernelUD)

## I used "writeOGR" and exported shapefiles because I was having trouble with raster(estUDm2spixdf()),
## But also used raster for "extract." Try GeoTIFF from writeOGR? This is kind of a mess...
## Can I just create separate raster objects for each animal instead of exporting them w/OGR?

##### TRY raster(as(...))

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
##    b. For UD height at occurrence points only
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
##    b. For UD height at occurrence points only
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
