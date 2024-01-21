rm(list=ls()) # Cleaning the shit off
gc()

library(DSSAT)
library(lubridate)
library(reshape)
library(ggplot2)
library(rasterVis)
library(RColorBrewer)
library(dplyr)
library(extrafont) 
loadfonts(device = "win")

library(reshape2)
library(grid)
library(rgdal)
library(raster)
library(sp)
library(maptools)

library(ggpolypath)
library(xlsx)

source(sprintf("%sHelper_pdates.R", "~/DATABANK/Floating/Dates/")) # my own library that i created 


setwd("~/DATABANK/Floating/Dates/")

admin <- readOGR("SHP/NGA_adm0.shp")    #Download Shp files https://www.diva-gis.org/datadown
region <- readOGR("SHP/NGA_adm1.shp")

loc.data <- read.csv("Ng_unique_names.csv")
loc.data.cut <- loc.data[c("Long","Lat","Location")]
colnames(loc.data.cut)[3] <- "site"
head(loc.data.cut)

M060.out <- read_output("Outputs/M060.OUT")
M060.yield <- dssat2yld(M060.out) # I created this function to shorten the code. mits the same as the old one but it repeated many times so this is better
colnames(M060.yield) <- c("site","year","M060")

M067.out <- read_output("Outputs/M067.OUT")
M067.yield <- dssat2yld(M067.out) # 
colnames(M067.yield) <- c("site","year","M067")

M074.out <- read_output("Outputs/M074.OUT")
M074.yield <- dssat2yld(M074.out) # 
colnames(M074.yield) <- c("site","year","M074")

M081.out <- read_output("Outputs/M081.OUT")
M081.yield <- dssat2yld(M081.out) # 
colnames(M081.yield) <- c("site","year","M081")

M088.out <- read_output("Outputs/M088.OUT")
M088.yield <- dssat2yld(M088.out) # 
colnames(M088.yield) <- c("site","year","M088")

M095.out <- read_output("Outputs/M095.OUT")
M095.yield <- dssat2yld(M095.out) # 
colnames(M095.yield) <- c("site","year","M095")

M102.out <- read_output("Outputs/M102.OUT")
M102.yield <- dssat2yld(M102.out) # 
colnames(M102.yield) <- c("site","year","M102")

M109.out <- read_output("Outputs/M109.OUT")
M109.yield <- dssat2yld(M109.out) # 
colnames(M109.yield) <- c("site","year","M109")

M116.out <- read_output("Outputs/M116.OUT")
M116.yield <- dssat2yld(M116.out) # 
colnames(M116.yield) <- c("site","year","M116")

M123.out <- read_output("Outputs/M123.OUT")
M123.yield <- dssat2yld(M123.out) # 
colnames(M123.yield) <- c("site","year","M123")

M130.out <- read_output("Outputs/M130.OUT")
M130.yield <- dssat2yld(M130.out) # 
colnames(M130.yield) <- c("site","year","M130")

M137.out <- read_output("Outputs/M137.OUT")
M137.yield <- dssat2yld(M137.out) # 
colnames(M137.yield) <- c("site","year","M137")

M144.out <- read_output("Outputs/M144.OUT")
M144.yield <- dssat2yld(M144.out) # 
colnames(M144.yield) <- c("site","year","M144")

M151.out <- read_output("Outputs/M151.OUT")
M151.yield <- dssat2yld(M151.out) # 
colnames(M151.yield) <- c("site","year","M151")

M158.out <- read_output("Outputs/M158.OUT")
M158.yield <- dssat2yld(M158.out) # 
colnames(M158.yield) <- c("site","year","M158")

M165.out <- read_output("Outputs/M165.OUT")
M165.yield <- dssat2yld(M165.out) # 
colnames(M165.yield) <- c("site","year","M165")

M172.out <- read_output("Outputs/M172.OUT")
M172.yield <- dssat2yld(M172.out) # 
colnames(M172.yield) <- c("site","year","M172")

M179.out <- read_output("Outputs/M179.OUT")
M179.yield <- dssat2yld(M179.out) # 
colnames(M179.yield) <- c("site","year","M179")

M186.out <- read_output("Outputs/M186.OUT")
M186.yield <- dssat2yld(M186.out) # 
colnames(M186.yield) <- c("site","year","M186")

M193.out <- read_output("Outputs/M193.OUT")
M193.yield <- dssat2yld(M193.out) # 
colnames(M193.yield) <- c("site","year","M193")

M200.out <- read_output("Outputs/M200.OUT")
M200.yield <- dssat2yld(M200.out) # 
colnames(M200.yield) <- c("site","year","M200")

M207.out <- read_output("Outputs/M207.OUT")
M207.yield <- dssat2yld(M207.out) # 
colnames(M207.yield) <- c("site","year","M207")

M214.out <- read_output("Outputs/M214.OUT")
M214.yield <- dssat2yld(M214.out) # 
colnames(M214.yield) <- c("site","year","M214")

M221.out <- read_output("Outputs/M221.OUT")
M221.yield <- dssat2yld(M221.out) # 
colnames(M221.yield) <- c("site","year","M221")

M228.out <- read_output("Outputs/M228.OUT")
M228.yield <- dssat2yld(M228.out) # 
colnames(M228.yield) <- c("site","year","M228")

M235.out <- read_output("Outputs/M235.OUT")
M235.yield <- dssat2yld(M235.out) # 
colnames(M235.yield) <- c("site","year","M235")

M242.out <- read_output("Outputs/M242.OUT")
M242.yield <- dssat2yld(M242.out) # I created this function to shorten the code. mits the same as the old one but it repeated many times so this is better
colnames(M242.yield) <- c("site","year","M242")

M249.out <- read_output("Outputs/M249.OUT")
M249.yield <- dssat2yld(M249.out) # 
colnames(M249.yield) <- c("site","year","M249")

M256.out <- read_output("Outputs/M256.OUT")
M256.yield <- dssat2yld(M256.out) # 
colnames(M256.yield) <- c("site","year","M256")

M263.out <- read_output("Outputs/M263.OUT")
M263.yield <- dssat2yld(M263.out) # 
colnames(M263.yield) <- c("site","year","M263")

M270.out <- read_output("Outputs/M270.OUT")
M270.yield <- dssat2yld(M270.out) # 
colnames(M270.yield) <- c("site","year","M270")

M277.out <- read_output("Outputs/M277.OUT")
M277.yield <- dssat2yld(M277.out) # 
colnames(M277.yield) <- c("site","year","M277")

M284.out <- read_output("Outputs/M284.OUT")
M284.yield <- dssat2yld(M284.out) # 
colnames(M284.yield) <- c("site","year","M284")

M291.out <- read_output("Outputs/M291.OUT")
M291.yield <- dssat2yld(M291.out) # 
colnames(M291.yield) <- c("site","year","M291")

M298.out <- read_output("Outputs/M298.OUT")
M298.yield <- dssat2yld(M298.out) # 
colnames(M298.yield) <- c("site","year","M298")

M305.out <- read_output("Outputs/M305.OUT")
M305.yield <- dssat2yld(M305.out) # 
colnames(M305.yield) <- c("site","year","M305")

M312.out <- read_output("Outputs/M312.OUT")
M312.yield <- dssat2yld(M312.out) # 
colnames(M312.yield) <- c("site","year","M312")


all.dates.combs <- cbind(M060.yield$M060,M067.yield$M067,M074.yield$M074,M081.yield$M081,M088.yield$M088,M095.yield$M095,M102.yield$M102,M109.yield$M109,
                         M116.yield$M116,M123.yield$M123,M130.yield$M130,M137.yield$M137,M144.yield$M144,M151.yield$M151,M158.yield$M158,M165.yield$M165,
                         M172.yield$M172,M179.yield$M179,M186.yield$M186,M193.yield$M193,M200.yield$M200,M207.yield$M207,M214.yield$M214,M221.yield$M221,
                         M228.yield$M228,M235.yield$M235,M242.yield$M242,M249.yield$M249,M256.yield$M256,M263.yield$M263,M270.yield$M270,M277.yield$M277,
                         M284.yield$M284,M291.yield$M291,M298.yield$M298,M305.yield$M305,M312.yield$M312)

colnames(all.dates.combs) <- c("site","year","M60","M67","M74","M81","M88","M95","M102","M109","M116","M123","M130","M137","M144","M151","M158","M165",
                               "M172","M179","M186","M193","M200","M207","M214","M221","M228","M235","M242","M249","M256","M263","M270","M277","M284","M291","M298","M305","M312")

#ANALYSIS 
all.dates.combs$bpdat <- apply(all.dates.combs, 1, which.max) # identify the max yield value
all.dates.combs$bpdat <- colnames(all.dates.combs)[all.dates.combs$bpdat] # identify the pdat with max yield
all.dates.combs$bpdat <- as.numeric(sub("^M", "", all.dates.combs$bpdat)) # convert it into a numeric value
all.dates.bdat <- all.dates.combs[,c("site", "year", "bpdat")] # remove trash 
head(all.dates.bdat, 39)

site.dates.combs <- merge(loc.data.cut, all.dates.bdat, by= "site")
head(site.dates.combs, 12)

#site.dates.combs.81 <- subset(site.dates.combs, year >= 1995 & year <=2014)
site.dates.combs.81 <- subset(site.dates.combs, year == 1981)
head(site.dates.combs.81)

centroids.df <- as.data.frame(coordinates(region))
names(centroids.df) <- c("lon", "lat")
centroids.df$label <- region$NAME_1[match(rownames(centroids.df), region$NAME_1)]
centroids.df$label <- region$NAME_1


#plot with ggplot
ggplot() +
  geom_tile(data = site.dates.combs.81, aes(x = Long, y = Lat, fill = bpdat)) + 
  scale_fill_gradientn("Planting date\n[DOY]", #values = c(0, 0.3,0.6, 1),  
                       colours=c(brewer.pal(9,'YlGn'))) + 
  geom_polygon(data=admin, aes(x = long, y=lat, group=group), fill = NA, col="black", size = 0.8) + 
  geom_polygon(data=region, aes(x = long, y=lat, group=group), fill = NA, col="grey42", size = 0.2) +
  #geom_text(data = centroids.df, aes(lon, lat, label = label), col="darkred", size = 3) +
  geom_hline(yintercept=seq(4,14,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  geom_vline(xintercept=seq(1,15,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  theme(axis.text=element_text(family="Calibri", size = rel(1.2)), 
        axis.title.x=element_text(margin=margin(t=0.4, unit="cm"), size=13,family="Calibri"),
        axis.title.y=element_text(margin=margin(r=0.4, unit="cm"), size=13, family="Calibri"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title =element_text(size = 13, face = 'bold', family="Calibri", margin=margin(b=0.4, unit="cm")),
        legend.key.width = unit(0.35, "cm"), 
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(family="Calibri", size = rel(1.2)),
        legend.position = "right")  +
  scale_x_continuous(breaks= seq(1,15,by=2), name= "Longitude (?E)") + 
  scale_y_continuous(breaks= seq(4,14,by=2), name="Latitude (?N)") + 
  coord_map(xlim = c(2,15), ylim = c(4,14))

#ggsave("./OUT/Current.jpg",width=4, height=3.5,units="in") 

head(site.dates.combs, 12)

site.81.91.01.11 <- subset(site.dates.combs, year == 1981|year == 1991|year == 2001|year == 2011)
unique(site.81.91.01.11$year)
head(site.81.91.01.11, 12)

#Plot
ggplot() +
  geom_tile(data = site.81.91.01.11, aes(x = Long, y = Lat, fill = bpdat)) + # the year is to allow for multiple years to be plotted in a panel, but for now its just one panel
  scale_fill_gradientn("Planting date\n[DOY]", values = c(0, 0.3,0.6, 1),  # values can twist the legend
                       colours=c(brewer.pal(9,'YlGn')), #colorramp for the colouring.The number is for levels. Can be changed to ""Greens", "PuBu", etc
                       na.value = 'grey') + #For na.a
  geom_polygon(data=admin, aes(x = long, y=lat, group=group), fill = NA, col="black", size = 0.8) + # for the country map
  geom_polygon(data=region, aes(x = long, y=lat, group=group), fill = NA, col="grey42", size = 0.2) + # for the country map
  #geom_text(data = centroids.df, aes(lon, lat, label = label), col="darkred", size = 3) +
  geom_hline(yintercept=seq(4,14,by=2), linetype="dashed", lwd=0.3, color = "grey62") + # for y-axis
  geom_vline(xintercept=seq(1,15,by=2), linetype="dashed", lwd=0.3, color = "grey62") + # for x-axis
  theme(axis.text=element_text(family="Calibri", size = rel(1.2)), #Text fonts 
        axis.title.x=element_text(margin=margin(t=0.4, unit="cm"), size=7,family="Calibri"),
        axis.title.y=element_text(margin=margin(r=0.4, unit="cm"), size=7, family="Calibri"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title =element_text(size = 13, face = 'bold', family="Calibri", margin=margin(b=0.4, unit="cm")),
        legend.key.width = unit(0.35, "cm"), 
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(family="Calibri", size = rel(1.2)),
        legend.position = "right")  +
  # geom_polygon(data=Ng.grid, aes(x = long, y=lat, group=group), fill = NA, col="grey80", size = 0.05) + # for the country map
  facet_grid(~year) +  # ndo pane mari yese apa
  scale_x_continuous(breaks= seq(1,15,by=2), name= "Longitude (?E)") + #for the map to fit, to remove background greys on x axis
  scale_y_continuous(breaks= seq(4,14,by=2), name="Latitude (?N)") + #for the map to fit, to remove background greys on x axis
  coord_map(xlim = c(2,15), ylim = c(4,14)) 





#####_____________________________________________________________________________________________________
print("Inini Ndinonzi Inini(_)")











