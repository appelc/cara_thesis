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

library(googlesheets)
library(adehabitatHR)
library(ruf)
library(raster)

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

## Calculate KDE using kernelUD for all animals
all.kud <- kernelUD(porc.sp, h="href")
image(all.kud) #href doesn't look great

## revisit extent parameter (actually calculate for cell ~3-6 meters)
all.kud <- kernelUD(porc.sp, h=60, extent=1, grid=1000)
image(all.kud)
points(porc.sp)

# convert estUD to raster
all.kud.raster <- raster(estUDm2spixdf(all.kud))
###### ERROR: this function can only be used when the same grid was used for all animals
###### (can add "same4all=TRUE" to kernelUD code, but I don't think that's what we want)

######################
## 3. Then, create a table with id, coord, and UD height
## (or do I want a different table for each porc?)
## a. Get UD height at each pixel
######################

ids <- names(all.kud)
porc_uds <- NULL

for(i in ids){
      ud.height.i <- all.kud[[i]]@data$ud
      coords.i <- all.kud[[i]]@coords
      ht.i <- cbind((rep(i, length(ud.height.i))), ud.height.i, coords.i)
      colnames(ht.i) <- c("id", "height", "x", "y")
      porc_uds <- rbind(porc_uds, ht.i)
}
write.csv(porc_uds, "porc_uds.csv")
## too big... do them separately!

######################
## b. Get UD height at occurrence points only
######################

# or you can extract just at the coordinates from the raster
porc.sp$udheight <- extract(all.kud.raster, porc.sp)

