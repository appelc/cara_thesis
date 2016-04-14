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
all.kud <- kernelUD(porc.sp, h="href")
image(all.kud) #href doesn't look great

## revisit extent parameter (actually calculate for cell ~3-6 meters)
all.kud <- kernelUD(porc.sp, h=60, extent=1, grid=1000)
image(all.kud)

######################
## 3. Then, create a table with id, coord, and UD height
## (or do I want a different table for each porc?)
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

write.csv(porc_uds, "porc_uds.csv")

## too big... do them separately!

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

