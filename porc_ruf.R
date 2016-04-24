#######################################
## Trying resource utilization functions (RUFs)
## with porcupine data
#######################################

## 1. First, load data
## 2. Then, extract the UD from "adehabitatHR" package (set bandwidth, grid, and extent)
## 3. Then, create a table with id, coord, and UD height
## 4. Assign values of covariates (veg class, canopy height) to cells
## 5. Run RUF using package "ruf"

## 3a. - 5a.: For UD height at each pixel
## 3b. - 5b.: For UD height at occurrence points only

### FROM CARA (4/22):
## i.   is there any way to make boundary work in kernelUD? (not the same as clipping post-hoc?)
##      -if not, also clip to study area / veg layer extent when I clip to 95% contour!
## ii.  find paper about Matern correlation parameters (theta)
## iii. ways to standardize/normalize/scale the UD height (response variable)?
## iv.  any way to double-check that my cell sizes are all the same for each animal?
## v.   clean up / consolidate steps... e.g., don't calcluate "norm" and "log" until end of step 4

### FROM TIM (4/19):
## 1. Convert heights to reasonable numbers
##  -center & scale: (x-mean)/sd ... but returns negative numbers for height=0
##  -log: but returns negative numbers because all heights are between 0-1
##  -normalize: (x - min) / (max - min)
# 2. Re-run Henrietta on 95%
# 2a. Run Henrietta wth Veg2
# 3. See if you can run on other porcupines
# 4. Theta??

library(adehabitatHR)
library(googlesheets)
library(raster)
library(rgdal)
#install.packages("ruf",repos="http://www.stat.ucla.edu/~handcock")
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
veg <- readOGR(dsn="F:/Shapefiles/Veg map", layer="Veg categories CA", verbose=TRUE)
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

## Calculate grid & extent based on desired cell size (# meters on each side)
## Do this for each animal separately 
## I tried using "boundary=study.ext" in kernelUD, but angles aren't acceptable
## (will need to clip grids later; but does it matter for calculating kernelUD?)

ids <- unique(porc.locs$id)
ud.list <- NULL
contour.list <- list()

for (i in ids){
  locs.i <- subset(porc.locs, id == i)
  sp.i <- SpatialPointsDataFrame(data.frame(locs.i$utm_e, locs.i$utm_n),
                                    data=data.frame(locs.i$id),
                                    proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
  c = 5   ## desired cell size (meters)
  fake.kern <- kernelUD(xy = sp.i, extent = 1)
  spdf <- raster(as(fake.kern[[1]], "SpatialPixelsDataFrame"))
  eas <- diff(range(spdf@extent[1:2]))
  nor <- diff(range(spdf@extent[3:4]))
    if(eas > nor){
      g <- (eas/c)
    } else {
      g <- (nor/c)
    }
  kern.i <- kernelUD(xy = sp.i, h = 60, grid = g, extent = 1)
  hr95.i <- getverticeshr.estUDm(kern.i, percent = 95, unin = "m", unout = "km2", standardize = FALSE)
  ud.list[[i]] <- kern.i
  contour.list[[i]] <- hr95.i
}

image(ud.list[[17]]) ## check a few... looks good!
plot(contour.list[[17]], add = TRUE)

######################
## 3. Then, create a list of tables with id, coord, and UD height for each porc
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
      height.list[[i]] <- data.frame(ht.coords.i) 
}

## Turn each into a SPDF and clip to individual 95% contours (made above) 
## **consider doing 99% instead (Long et al.)**
## outputs as a list of SPDFs (still containing height, normalized height, etc.)
## now the range of values should be much better

spdf95.list <- list()

for(i in ids){
      ht.i <- height.list[[i]]
      sp.i <- SpatialPointsDataFrame(data.frame(ht.i$x, ht.i$y),
                                     data=data.frame(ht.i$id, ht.i$height),
                                     proj4string = CRS(proj4string(veg)))
      contour.i <- contour.list[[i]]
      spdf.i <- sp.i[contour.i,] # clip to 95% contour
      spdf.i <- spdf.i[veg,] # also clip to extent of veg layer
      spdf95.list[[i]] <- spdf.i
}

plot(spdf95.list[[17]]) ## look at a couple

## Alternative way to subset 95% based on 5% quantile values:
## doesn't really work because of so many zeroes
#ht95.list <- NULL
#for (i in ids){
#  ht.i <- height.list[[i]]
#  x <- quantile(ht.i$ht.log, 0.05)
#  y <- round(x[[1]], 14) ## next line only works if I round the quantile value
#  ht95.i <- ht.i[ht.i$ht.norm >= y,]
#  ht95.list[[i]] <- data.frame(ht95.i)
#}

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
## and include a column for normalizing the UD height: (x-min)/(max-min)
## (and/or log of height)
##    a. For UD height at each pixel
######################

## now assign veg class to each cell (row)
## (and eventually canopy height)

## this loop uses the SPDF created in the previous step, does 'overlay' with veg class, 
## gets rid of cells where veg=NA (there shouldn't be many but they may mess up ruf.fit),
## then turns it back into a data frame for calculating RUF in next step

ids <- unique(porc.locs$id)
final.list <- NULL

for (i in ids){
        spdf.i <- spdf95.list[[i]]
        spdf.i@data$veg <- over(spdf.i, veg)$Class_2
        df.i <- data.frame(i, spdf.i@data$ht.i.height, spdf.i@coords, spdf.i@data$veg)
        colnames(df.i) <- c("id", "ud", "x", "y", "veg")
        df.i <- df.i[!is.na(df.i$veg),]
        min <- min(df.i$ud)
        max <- max(df.i$ud)
        df.i$ud_norm <- ((df.i$ud) - min) / (max - min)
        df.i$ud_log <- log(df.i$ud)      
        final.list[[i]] <- df.i
}

######################
## 5. Run RUF using package "ruf"
##    a. For UD height at each pixel
######################

## Need to switch from 64-bit R to 32-bit R now for the "ruf" package... 
## *hopefully* everything will be stored in the environment after restarting...
## (actually, everything so far seems to run OK on 32-bit version)

library(ruf)

## Now, "final.list" contains the data frames necessary to run ruf.fit
## (id, normalized/log ud height for top 95%, x, y, veg class)

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)

ids <- unique(porc.locs$id)
ruf.list <- list()
betas.list <- list()

for (i in ids){
        df.i <- final.list[[i]]
        ruf.i <- ruf.fit(ud_norm ~ factor(veg),
                         space = ~ x + y,
                         data = df.i, theta = hval,
                         name = i,
                         standardized = F)
        ruf.list[[i]] <- ruf.i
        betas.list[[i]]  <- ruf.i$beta
        path <- file.path("F:", "RUF", paste(i, "_betas", ".csv", sep = ""))
        write.csv(betas.list[[i]], file=path)
}

write.csv(betas.list[[1]], file="15.01_ruf_betas_042316")

## For Henrietta, all the betas are positive using normalized UD heights. Does this make sense?
plot(betas.list[[1]])
## Look at distribution of normalized UD heights:
hist(final.list[[1]]$ud_norm)

############################################################################3
## STEPS 3-5
## b. for height at occurrence points only
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
  locs.i <- subset(porc.locs, id == i)
  locs.i <- droplevels(locs.i)
  sp.i <- SpatialPointsDataFrame(data.frame(locs.i$utm_e, locs.i$utm_n),
                                 data=data.frame(locs.i$id),
                                 proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))
  i.kud <- kernelUD(sp.i, h=60, extent=1, grid=1000)
  i.raster <- raster(estUDm2spixdf(i.kud))
  sp.i$udheight <- extract(i.raster, sp.i)
  name = paste(i, "raster")
  filepath <- file.path("C:","Users","Cara","Documents","cara_thesis", "rasters",
                        paste(i, "raster", ".shp", sep = ""))
  writeOGR(sp.i, dsn=filepath, layer=name, driver="ESRI Shapefile") ## not a raster!
  export <- data.frame(sp.i@data$locs.i.id, sp.i@data$udheight, sp.i@coords)
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
