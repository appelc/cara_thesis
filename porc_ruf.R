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
porc.locs$utm_e <- as.numeric(porc.locs$utm_e)
porc.locs$utm_n <- as.numeric(porc.locs$utm_n)
#porc.locs$date <- as.Date(porc.locs$date, "%m/%d/%Y") ## check date format before running line 51 or 52
#porc.locs$date <- as.Date(porc.locs$date, origin = as.Date("1899-12-30"))

## OPTIONAL: only keep summer locations (before Nov 1)
## could even incorporate this into the for-loop below
sum.locs <- subset(porc.locs, date < "2015-11-01")

## Keep only animals with >= 5 locations
n <- table(porc.locs$id)
porc.locs <- subset(porc.locs, id %in% names(n[n >= 5]), drop=TRUE)
porc.locs <- droplevels(porc.locs)

n <- table(sum.locs$id)
sum.locs <- subset(sum.locs, id %in% names(n[n >= 5]), drop=TRUE)
sum.locs <- droplevels(sum.locs)

## Turn these into a Spatial Points Data Frame
## I never actually use these!! Other than to set projection for "veg"...
porc.sp <- SpatialPointsDataFrame(data.frame(porc.locs$utm_e, porc.locs$utm_n),
                                  data=data.frame(porc.locs$id),
                                  proj4string=CRS("+proj=utm +zone=10 +datum=NAD83"))

sum.sp <- SpatialPointsDataFrame(data.frame(sum.locs$utm_e, sum.locs$utm_n),
                                 data=data.frame(sum.locs$id),
                                 proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

## Load veg data
veg <- readOGR(dsn="D:/GIS DATA/Veg map", layer="Veg categories CA", verbose=TRUE)
proj4string(veg) <- proj4string(porc.sp)

######################
## 2. Then, extract the UD from "adehabitatHR" package
###################### 

## Calculate grid & extent based on desired cell size (# meters on each side)
## For for each animal separately 

## Also calculate KUD based on summer points ONLY, but within grid of the extent
## for all of the points. Then clip to the 99% contour for all the points, as well
## as the veg layer extent.

ids <- unique(porc.locs$id)
ud.list <- list()
ud.summer.list <- list()
ud.clipped.list <- list()
contour.list <- list()
contour.summer.list <- list()

for (i in ids){
        locs.i <- subset(porc.locs, id == i)
        locs.i$id_season <- rep(paste(i, "_all", sep = ""), nrow(locs.i))
        locs.sum.i <- subset(sum.locs, id == i)
        locs.sum.i$id_season <- rep(paste(i, "_sum", sep = ""), nrow(locs.sum.i))
        locs.all.i <- rbind(locs.i, locs.sum.i)
        sp.i <- SpatialPointsDataFrame(data.frame(locs.all.i$utm_e, locs.all.i$utm_n),
                                    data=data.frame(locs.all.i$id_season),
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
  
        # calculate UD on both IDs ("all" and "summer") with same4all = TRUE
        kern.i <- kernelUD(xy = sp.i, h = 60, grid = g, extent = 1, same4all = TRUE)
        cont99.i <- getverticeshr.estUD(kern.i[[1]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)
        cont99.sum.i <- getverticeshr.estUD(kern.i[[2]], percent = 99, unin = "m", unout = "km2", standardize = FALSE)
        
        # clip summer UD to 99% contour from ALL points (not just summer), and veg extent
        sum.ud.i <- (kern.i[[2]])[cont99.i,]
        sum.ud.i <- sum.ud.i[veg,]
      
        # save full UD, summer UD, and clipped UD:
        ud.list[[i]] <- kern.i[[1]]
        ud.summer.list[[i]] <- kern.i[[2]]
        ud.clipped.list[[i]] <- sum.ud.i ##it's now a "SpatialPixelsDataFrame"
        
        # and save the contours:
        contour.list[[i]] <- cont99.i
        contour.summer.list[[i]] <- cont99.i
}

# Why is grid size not exactly 5? Try calculating "g" based on SPDF, not raster
# Extent of spatialpointsdataframe:  
#  eas <- diff(range(extent(sp.i)[1:2]))
#  nor <- diff(range(extent(sp.i)[3:4]))

######################
## 3. Then, create a list of tables with id, coord, and UD height for each porc
##    a. For UD height at each pixel
######################

ids <- unique(sum.locs$id)
height.list <- list()

for(i in ids){
      ud.i <- ud.clipped.list[[i]]
      ud.height.i <- ud.i$ud
      coords.i <- ud.i@coords
      ht.coords.i <- data.frame((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.coords.i) <- c("id", "height", "x", "y")
      height.list[[i]] <- data.frame(ht.coords.i) 
}

######################
## 4. Assign values of covariates (veg class, canopy height) to cells
## and include a column for normalizing the UD height: (x - min) / (max - min)
## (and/or log of UD height)
##    a. For UD height at each pixel
######################

## this loop creates a SPDF for each animal, does 'overlay' with veg class, 
## gets rid of cells where veg=NA (there shouldn't be many but they may mess up ruf.fit),
## then turns it back into a data frame for calculating RUF in next step

ids <- unique(sum.locs$id)
final.list <- list()

for (i in ids){
        ht.i <- height.list[[i]]
        sp.i <- SpatialPointsDataFrame(data.frame(ht.i$x, ht.i$y),
                                       data=data.frame(ht.i$id, ht.i$height),
                                       proj4string = CRS(proj4string(veg)))
        spdf.i@data$veg <- over(spdf.i, veg)$Class_2
        df.i <- data.frame(i, spdf.i@data$ht.i.height, spdf.i@coords, spdf.i@data$veg)
        colnames(df.i) <- c("id", "height", "x", "y", "veg")
        df.i <- df.i[!is.na(df.i$veg),]
        min <- min(df.i$height)
        max <- max(df.i$height)
        df.i$height_norm <- ((df.i$height) - min) / (max - min)
        df.i$height_log <- log(df.i$height)      
        final.list[[i]] <- df.i
}

######################
## 5. Run RUF using package "ruf"
##    a. For UD height at each pixel
######################

## Now, "final.list" contains the data frames necessary to run ruf.fit
## (id, normalized/log UD height for top 95%, x, y, veg class)

## Set initial estimates for range/smoothness
hval <- c(0.2, 1.5)

ids <- unique(sum.locs$id)
ruf.list <- list()
thetas.list <- list()
fit.list <- list()
betas.table <- NULL

for (i in ids){
        df.i <- final.list[[i]]
        ruf.i <- ruf.fit(height_norm ~ factor(veg),
                         space = ~ x + y,
                         data = df.i, name = i, standardized = F, theta = hval,
                         fixrange = FALSE, fixsmoothness = FALSE)
        ruf.list[[i]] <- ruf.i
        thetas.list[[i]] <- ruf.i$theta
        fit.list[[i]] <- ruf.i$fit
        betas.table <- bind_rows(betas.table, ruf.i$beta)
        #path <- file.path("F:", "RUF", paste(i, "_betas", ".csv", sep = ""))
        #write.csv(betas.list[[i]], file=path)
}

## Look at beta values for each animal (then get averages from betas.table)
plot(ruf.i$beta)
barplot(ruf.i$beta, names.arg=ruf.i$veg, ylab="beta coefficient", xlab="veg class", las=2)

## Look at distribution of normalized UD heights:
hist(final.list[[3]]$height_norm)


############################################################################
## STEPS 3-5
## b. for height at occurrence points only
############################################################################

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
