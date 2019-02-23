porc.new <- subset(porc.vhf, type %in% c('N'))
porc.new$utm_e <- as.numeric(porc.new$utm_e)
porc.new$utm_n <- as.numeric(porc.new$utm_n)
porc.new$date <- as.Date(porc.new$date, "%m/%d/%Y") 


## subset summer locations (before Nov 1 or after March 1) and winter (between Nov 1 and March 1)
sum.cutoff <- '2015-11-01' 
win.cutoff <- '2016-03-01'
sum.locs.new <- porc.new[(porc.new$date < sum.cutoff) | (porc.new$date >= win.cutoff), ]
win.locs.new <- porc.new[(porc.new$date >= sum.cutoff) & (porc.new$date < win.cutoff), ]

## export points as shapefile for making figure in ArcMap
porc.new.sp <- SpatialPointsDataFrame(data.frame(porc.new$utm_e, porc.new$utm_n),
                                       data = data.frame(porc.new),
                                       proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
writeOGR(porc.new.sp, dsn = 'Shapefiles/porc_new_110916', layer = 'porc_locs_110916', driver = 'ESRI Shapefile')

## summer points
porc.new.sp.sum <- SpatialPointsDataFrame(data.frame(sum.locs.new$utm_e, sum.locs.new$utm_n),
                                      data = data.frame(sum.locs.new),
                                      proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
writeOGR(porc.new.sp.sum, dsn = 'Shapefiles/porc_new_110916_sum', layer = 'porc_locs_110916_sum', driver = 'ESRI Shapefile')

## winter points
porc.new.sp.win <- SpatialPointsDataFrame(data.frame(win.locs.new$utm_e, win.locs.new$utm_n),
                                          data = data.frame(win.locs.new),
                                          proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
writeOGR(porc.new.sp.win, dsn = 'Shapefiles/porc_new_110916_win', layer = 'porc_locs_110916_win', driver = 'ESRI Shapefile')
