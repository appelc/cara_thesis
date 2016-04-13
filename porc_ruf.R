#######################################
## Trying resource utilization functions (RUFs)
## with porcupine data
##
## based on RUF lab 8 from Tim's class
#######################################

library(adehabitatHR)
library(ruf)
library(googlesheets)

## 1. First, load data
## 2. Then, need to extract the UD from "adehabitatHR" package
## 3. Then, create a table with id, coord, and ud height
## 4. Run RUF using package "ruf"

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

## try first with just a couple animals
hen.locs <- subset(porc.locs, id %in% c("15.01"))
ha.sp <- SpatialPointsDataFrame(data.frame(hen.locs$utm_e, hen.locs$utm_n),
                                 data=data.frame(hen.locs$id),
                                 proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))

## Calculate KDE using kernelUD first and then kernel.area
ha.kud <- kernelUD(ha.sp, h="LSCV", hlim=c(0.01,2)) #did not converge
plotLSCV(ha.kud)
ha.kud <- kernelUD(ha.sp, h="href") #no error

## Now try extent/grid combos
porc1e <- kernelUD(ha.sp, extent=0.2, grid=1000, h=60) 
image(porc1e)
points(ha.sp)

## Now do with all
all.kud <- kernelUD(porc.sp, h="href")
image(all.kud)

all.kud <- kernelUD(porc.sp, h=60, extent=0.2, grid=100)
image(all.kud)
points(porc.sp)



dups <- duplicated(porc.sp@coords) # find duplicates in coordinates
head(dups) # returns list of TRUE/FALSE
dup.loc <- porc.sp[dups,] # Create table with duplicates only
dup.loc

non.dups <- porc.sp[!dups,] ## Now remove the duplicates; "non.dups" is our new data name






# Set up graphical window to see multiple graphs
par(mfrow=c(2,2)) 
image(porc1a)
title(main="grid=20, extent=0.2")
image(porc1b)
title(main="grid=100, extent=0.2") ## I think I like this the best
image(porc1c)
title(main="grid=20, extent=2")
image(porc1d)
title(main="grid=100, extent=2")
image(porc1e)
title(main="grid=100, extent=1") ## compromise b/c error with 1b (extent too small)

# convert estUD to raster
porc1e.raster <- raster(estUDm2spixdf(porc1e))

# one way to extract all the values from estUD                        
ud.height1 <- ha.kud[[1]]@data$ud

# or you can extract just at the coordinates from the raster
ha.sp$udheight <- extract(porc1e.raster, ha.sp)



# The actual data.frame that you'd want to create should have
# the same columns (i.e. the estimates of the "height" of the utilization
# distribution; the x and y coordinates of each obesrvation; and the
# values of the independent variables at each coordinate.). The difference
# is that rather than including just 100 points, you'd include every single
# cell within the 95% (or 100%) KDE.

# Hint: You can extract the "height" of the UD, once you've created one
# in adehabitat, with kud@data$ud (where 'kud' is the name of the kernel
# density estimate you created). You can extract the coordinates with
# kud@coords. 