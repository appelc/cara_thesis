## get climate data across the porcupine's range
## precip and temp by month from Worldclim

## SHOULD DIGITIZE A DIFFERENT RANGE MAP... THIS ONE DOESN'T INCLUDE CALIFORNIA (OR TEXAS)!
## UPDATED HALL 1981? OR ROZE AND ILSE 2003?

library(dismo)
library(XML)
library(sp)
library(rgdal)
library(raster)
library(reshape2)
library(ggplot2)
library(directlabels)
library(ggrepel)
library(gridExtra)

# Get temp and precip data from World Clim (use 2.5 min so don't have to download tiles)
precip <- getData('worldclim', var = 'prec', res = 2.5)
  plot(precip) ## avg monthly precip (mm) for each month
temp <- getData('worldclim', var = 'tmean', res = 2.5)
  plot(temp) ## avg monthly temp (degrees C x 10) for each month

# Import the porcupine range map shapefile as a Spatial Points Data Frame
## (CAN SKIP THIS STEP)
  porc_range <- readOGR(dsn = 'D:/GIS DATA/Range maps', layer = 'porcupine_range_dissolved')
class(porc_range)

# Now crop climate data to the extent of the porcupine range map
porc_range <- spTransform(porc_range, precip@crs) ## make sure they're in the same coordinate reference system
  precip_cropped <- crop(precip, porc_range)
  precip_porc <- mask(precip_cropped, porc_range) #takes a while
  plot(precip_porc)

porc_range <- spTransform(porc_range, temp@crs) ## ensure same CRS
  temp_cropped <- crop(temp, porc_range)
  temp_porc <- mask(temp_cropped, porc_range) #takes a while
  plot(temp_porc)

# Extract values
precip_porc
## why only min and max?

##########
## PLOT: don't need porc range map for this (unless we want to avg across)

# SITE SOURCES: TDSP, Coltrane and Sinnott, Roze, Sweitzer and Berger, Morin et al., Ilse and Hellgren, 
# Pokallus and Pauli, Craig and Keller, Griesemer et al., Mally, Tenneson and Oring, Dodge and Barnes)
lats <- c(41.883, 61.72, 42.25, 40.733, 48.35, 29.618, 44.308, 42.75, 42.406, 46.563, 54.583, 47.197, 46.26)
lons <- c(-124.20, -150.02, -74.417, -119.333, -68.767, -100.451, -90.131, -112.333, -72.341, -114.083, -128.7, -95.164, -122.32)
sites <- c('CA', 'AK', 'NY', 'NV', 'QC', 'TX', 'WI', 'ID', 'MA', 'MT', 'BC', 'MN', 'WA') #abbreviated names of locations

pts <- SpatialPoints(coords = data.frame(x = lons, y = lats), proj4string = temp@crs)

precip_values <- extract(precip, pts)
  rownames(precip_values) <- sites
  precip_values <- (precip_values)/10 ## convert mm to cm
matplot(t(precip_values), type = 'l', xlab = 'Month', ylab = 'Mean precip (mm)', lwd = 2)
  legend('top', rownames(precip_values), col = seq_len(ncol(precip_values)), cex = 0.8, fill = seq_len(ncol(precip_values)))

temp_values <- extract(temp, pts)  
  rownames(temp_values) <- rownames(precip_values)
  temp_values <- (temp_values)/10 ## World Clim uses a scaling factor of 10
matplot(t(temp_values), type = 'l', xlab = 'Month', ylab = 'Mean temp (C)', lwd = 2)
  legend('bottom', rownames(temp_values), col = seq_len(ncol(temp_values)), cex = 0.8, fill = seq_len(ncol(temp_values)))
  
  
## in ggplot2:

colnames(precip_values) <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
precip_df <- data.frame(melt(precip_values))

## try with ID, MA, and MN removed (bc very similar to MT, NY, and WI, respectively) & MT for clarity
precip_df <- precip_df[precip_df$Var1 != 'ID' & precip_df$Var1 != 'MA' & precip_df$Var1 != 'MN' & precip_df$Var1 != 'MT',]

precip_plot <- qplot(x = Var2, y = value, color = Var1, data = precip_df, geom = 'line', group = Var1) + 
            geom_line(size = 1.3, colour = 'gray65') + #scale_colour_grey(start = 0.3, end = .7, guide = 'none') +  #enable scale_colour_gray and remove colour = 'gray60' to have diff shades of grayscale
            geom_line(data = precip_df[precip_df$Var1 == 'CA',], linetype = 'solid', colour = 'white', size = 1.3) +
            geom_line(data = precip_df[precip_df$Var1 == 'CA',], linetype = 'dashed', colour = 'black', size = 1.6) +
            ylab('Precipitation (cm)') +
            theme(axis.text = element_text(size=20, colour = 'black'),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 22, color = 'black', margin = margin(0,15,0,0)),
                    axis.line.x = element_blank(),
                    axis.line.y = element_line(size = 1, colour = 'black'),
                    axis.ticks.y = element_line(size = 1, colour = 'black'),
                    axis.ticks.x = element_blank(),
                    axis.ticks.length = unit(.3, 'cm'),
                    panel.background = element_rect(fill = 'white'),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
                    legend.position = 'none') +
            geom_text(label = 'A', aes(x = 1.3, y = 44), size = 12, colour = 'black') +
          #  geom_text_repel(data = precip_df[precip_df$Var2 == 'Dec',], aes(label = Var1), 
           #                 nudge_x = 1.3, size = 6, color = 'black', box.padding = unit(0.7, 'lines'),
            #                segment.alpha = 0.6, min.segment.length = unit(0, 'lines'), force = 1.5, 
             #               arrow = arrow(ends = 'last', angle = 20, length = unit(0.03, 'npc'), type = 'open')) +
            #geom_text_repel(data = precip_df[precip_df$Var2 == 'Jan' & precip_df$Var1 == 'MA' | 
           #                                   precip_df$Var2 == 'Jan' & precip_df$Var1 == 'ID' |
            #                                  precip_df$Var2 == 'Jan' & precip_df$Var1 == 'MT' |
             #                                 precip_df$Var2 == 'Jan' & precip_df$Var1 == 'NV' |
              #                                precip_df$Var2 == 'Jan' & precip_df$Var1 == 'MN',], 
               #             aes(label = Var1), nudge_x = -0.7, size = 6, color = 'black') +
            #geom_text_repel(data = precip_df[precip_df$Var2 == 'Dec' & precip_df$Var1 != 'MA' &
                #                               precip_df$Var2 == 'Dec' & precip_df$Var1 != 'ID' &
                 #                              precip_df$Var2 == 'Dec' & precip_df$Var1 != 'MT' &
                  #                             precip_df$Var2 == 'Dec' & precip_df$Var1 != 'NV' &
                   #                            precip_df$Var2 == 'Dec' & precip_df$Var1 != 'MN',], 
                    #        aes(label = Var1), nudge_x = 0.7, size = 6, color = 'black') +
            scale_x_discrete(expand = c(0.02, 0.3)) + 
            expand_limits(x = 13) + ## need both?
            expand_limits(y = -2)

precip_plot
## ID and MT very similar; NY and MA very similar; WI and MN very similar

#direct.label(precip_plot, 'angled.boxes') ##just testing

# use EITHER the geom_text_repel code above OR the following for labels:
## use package 'directlabel' to add labels (instead of legend, since this will be black & white)
qp.break <- qp.labels('y', 'bottom', 'top', make.tiebreaker('x', 'y')) #to ensure no overlap of labels
tiebreak <- list('last.points', 'calc.boxes', qp.break) 
direct.label(precip_plot, tiebreak)

## Enlarge the text size and spacing:
precip_dim <- list('last.points', cex = 1.5, colour = 'black', 'calc.boxes', dl.trans(h = 1.25*h), 'calc.borders', 'qp.break') ## I know colour='black' defeats the purpose but it's for pub
precip_plot_label <- direct.label(precip_plot, precip_dim)
precip_plot_label

## export 850 x 550 so that labels aren't hidden behind x-axis

## same for temp:
colnames(temp_values) <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
temp_df <- data.frame(melt(temp_values))

## try with ID, MA, and MN removed (bc very similar to MT, NY, and WI, respectively) & MT for clarity
temp_df <- temp_df[temp_df$Var1 != 'ID' & temp_df$Var1 != 'MA' & temp_df$Var1 != 'MN' & temp_df$Var1 != 'MT',]

temp_plot <- qplot(x = Var2, y = value, color = Var1, data = temp_df, geom = 'line', group = Var1) + 
          geom_line(size = 1.3, colour = 'gray65') + #scale_colour_grey(start = 0.3, end = .7, guide = 'none') +
          geom_line(data = temp_df[temp_df$Var1 == 'CA',], linetype = 'solid', colour = 'white', size = 1.3) +
          geom_line(data = temp_df[temp_df$Var1 == 'CA',], linetype = 'dashed', colour = 'black', size = 1.6) +
          ylab('Temperature (°C)') + 
          theme(axis.text = element_text(size = 20, colour = 'black'),
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 22, color = 'black', margin = margin(0,15,0,0)),
                axis.line.x = element_line(size = 1, colour = 'black'),
                axis.line.y = element_line(size = 1, colour = 'black'),
                axis.ticks = element_line(size = 1, colour = 'black'),
                axis.ticks.length = unit(.3, 'cm'),
                panel.background = element_rect(fill = 'white'),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
                legend.position = 'none') +
          geom_text(label = 'B', aes(x = 1.3, y = 28), size = 12, colour = 'black') +
      #    geom_text_repel(data = temp_df[temp_df$Var2 == 'Dec',], aes(label = Var1),
       #                   nudge_x = 1.3, size = 6, color = 'black', box.padding = unit(0.7, 'lines'),
        #                  arrow = arrow(ends = 'last', angle = 20, length = unit(0.03, 'npc'), type = 'open'),
         #                 segment.alpha = 0.5, min.segment.length = unit(0, 'lines'), force = 1) +
          #geom_text_repel(data = temp_df[temp_df$Var2 == 'Jan' & temp_df$Var1 == 'MA' | 
           #                                  temp_df$Var2 == 'Jan' & temp_df$Var1 == 'ID' |
            #                                 temp_df$Var2 == 'Jan' & temp_df$Var1 == 'MT' |
             #                                temp_df$Var2 == 'Jan' & temp_df$Var1 == 'NV' |
              #                               temp_df$Var2 == 'Jan' & temp_df$Var1 == 'MN',], 
               #           aes(label = Var1), nudge_x = -0.7, size = 6, color = 'black') +
          #geom_text_repel(data = temp_df[temp_df$Var2 == 'Dec' & temp_df$Var1 != 'MA' &
           #                                  temp_df$Var2 == 'Dec' & temp_df$Var1 != 'ID' &
            #                                 temp_df$Var2 == 'Dec' & temp_df$Var1 != 'MT' &
             #                                temp_df$Var2 == 'Dec' & temp_df$Var1 != 'NV' &
              #                               temp_df$Var2 == 'Dec' & temp_df$Var1 != 'MN',], 
               #           aes(label = Var1), nudge_x = 0.7, size = 6, color = 'black') +
          scale_x_discrete(expand=c(0.02, 0.3)) +  
          expand_limits(x = 13) + ## need both?
          expand_limits(y = -14)
temp_plot

direct.label(temp_plot, tiebreak)
temp_dim <- list('last.points', cex = 1.5, colour = 'black', 'calc.boxes', dl.trans(h = 1.25*h), 'calc.borders', 'qp.break') ## I know 'colour='black' defeats the purpuse but it's for pub
temp_plot_label <- direct.label(temp_plot, temp_dim)
temp_plot_label

## combine figures (export at 700 wide x 960 high) 

grid.arrange(precip_plot_label, temp_plot_label, ncol = 1) 

grid.arrange(precip_plot, temp_plot, ncol = 1) ## no labels

#####################
## temp at our study area?

temp_values
  # warmest month (Aug) = 15.5
  # coldest month (Dec) = 8.7

sum(precip_values[1,])
  # 191.5


# https://stackoverflow.com/questions/27796583/how-to-add-colour-matched-legend-to-a-r-matplot
# also helpful: https://gis.stackexchange.com/questions/227585/how-to-use-r-to-extract-data-from-worldclim
  
  
  
  