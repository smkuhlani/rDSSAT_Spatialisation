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
library(rcolors)
library(ggpolypath)

options(java.parameters = "-Xmx8000m")
library(xlsx)


#IMPORTANT!!!!!!!
# source(sprintf("%sHelper_pdates.R", "C:/Users/Admin/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/Instructions to Run n Batch/Processing/Dates/")) # my own library that i created 
# source(sprintf("%sHelper_pdates.R", "C:/Users/Admin/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/Instructions to Run n Batch/Processing/Dates/"))
# 
# setwd("C:/Users/Admin/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/Instructions to Run n Batch/Processing/Dates/")

# source(sprintf("%sHelper_pdates.R", "~/DATABANK/Nigeria/")) # my own library that i created 
source(sprintf("%sHelper_pdates.R", "D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/Processing/Dates/"))

# setwd("~/DATABANK/Nigeria/")
setwd("D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/Processing/Dates/")

admin <- readOGR("SHP/NGA_adm0.shp")    #Download Shp files https://www.diva-gis.org/datadown
region <- readOGR("SHP/NGA_adm1.shp")
agrozones <- readOGR("SHP/Nigeria_AEZ.shp")
ng.hollow <- readOGR("SHP/NGA_hollow.shp")
ng.grid <- readOGR("SHP/NGA_grid.shp")
ng.box<- readOGR("SHP/Box_ng.shp")

loc.data <- read.csv("Ng_unique.csv")
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
M242.yield <- dssat2yld(M242.out) 
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

all.dates.combs.m <- cbind(M060.yield$site,M060.yield$year,M060.yield$M060,M067.yield$M067,M074.yield$M074,M081.yield$M081,M088.yield$M088,M095.yield$M095,M102.yield$M102,M109.yield$M109,
                           M116.yield$M116,M123.yield$M123,M130.yield$M130,M137.yield$M137,M144.yield$M144,M151.yield$M151,M158.yield$M158,M165.yield$M165,
                           M172.yield$M172,M179.yield$M179,M186.yield$M186,M193.yield$M193,M200.yield$M200,M207.yield$M207,M214.yield$M214,M221.yield$M221,
                           M228.yield$M228,M235.yield$M235,M242.yield$M242,M249.yield$M249,M256.yield$M256,M263.yield$M263,M270.yield$M270,M277.yield$M277,
                           M284.yield$M284,M291.yield$M291,M298.yield$M298,M305.yield$M305,M312.yield$M312)

all.dates.combs.m<-as.data.frame(all.dates.combs.m)
colnames(all.dates.combs.m) <- c("site","year","M60","M67","M74","M81","M88","M95","M102","M109","M116","M123","M130","M137","M144","M151","M158","M165",
                                 "M172","M179","M186","M193","M200","M207","M214","M221","M228","M235","M242","M249","M256","M263","M270","M277","M284","M291","M298","M305","M312")

head(all.dates.combs.m)
#all.dates.combs.m$max.yield <- max(all.dates.combs.m[,3:39])
#write.csv(all.dates.combs.m, 'Posterior/NG_sowing_dates.csv')

#ANALYSIS 
all.dates.combs.m$bpdat <- apply(all.dates.combs.m, 1, which.max) # identify the max yield value
all.dates.combs.m$bpdat <- colnames(all.dates.combs.m)[all.dates.combs.m$bpdat] # identify the pdat with max yield
all.dates.combs.m$bpdat <- as.numeric(sub("^M", "", all.dates.combs.m$bpdat)) # convert it into a numeric value
all.dates.bdat <- all.dates.combs.m[,c("site", "year", "bpdat")] # remove trash 
head(all.dates.bdat, 39)


#write.csv(site.dates.combs, 'Posterior/Series_best_dates.csv')

# xSelect one year for visualisation 
head(all.dates.combs.m)

all.dates.combs.m$max_yield <- all.dates.combs.m[3:39][cbind(seq_len(nrow(all.dates.combs.m)), max.col(all.dates.combs.m[3:39]))]
all.dates.max.y <- all.dates.combs.m[,c("site", "year", "max_yield")] # remove trash 
head(all.dates.max.y, 39)

#combine with the lonlat data
site.dates.combs <- merge(loc.data.cut, all.dates.max.y, by= "site")
head(site.dates.combs, 12)
tail(site.dates.combs, 12)


site.dates.combs.00 <- subset(site.dates.combs, year == 2017)
head(site.dates.combs.00)
site.dates.combs.00$max_yield <- as.numeric(site.dates.combs.00$max_yield)

#required for putting labels on the map
centroids.df <- as.data.frame(coordinates(agrozones))
names(centroids.df) <- c("lon", "lat")
centroids.df$label <- agrozones$AEZ2[match(rownames(centroids.df), agrozones$AEZ2)]
centroids.df$label <- agrozones$AEZ2

#plot with ggplot
ggplot() +
  geom_tile(data = site.dates.combs.00, aes(x = Long, y = Lat, fill = max_yield)) + 
  scale_fill_gradientn("Maize yield\n[Kg/ha]", #values = c(0, 0.3,0.6, 1),  
                       colours=rcolors$GreenYellow,
                       na.value = 'grey') +
  geom_polygon(data=admin, aes(x = long, y=lat, group=group), fill = NA, col="black", size = 1.5) + 
  geom_polygon(data=ng.grid, aes(x = long, y=lat, group=group), fill = NA, col="grey82", size = 0.2) +
  geom_polygon(data=agrozones, aes(x = long, y=lat, group=group), fill = NA, linetype="dashed", col="grey42", size = 0.2) +
  geom_polypath(data=ng.hollow, aes(x = long, y=lat, group=group), fill = "white", col="grey62", size = 0.2) + # for the country map
  geom_text(data = centroids.df, aes(lon, lat, label = label), col="darkred", size = 3) +
  geom_hline(yintercept=seq(4,14,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  geom_vline(xintercept=seq(1,15,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  theme(axis.text=element_text(family="DejaVu Sans", size = rel(0.8)),
        axis.title.x=element_text(margin=margin(t=0.4, unit="cm"), size=11,family="DejaVu Sans"),
        axis.title.y=element_text(margin=margin(r=0.4, unit="cm"), size=11, family="DejaVu Sans"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title =element_text(size = 12, face = 'plain', family="DejaVu Sans", margin=margin(b=0.4, unit="cm")),
        legend.key.width = unit(0.35, "cm"), 
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(family="DejaVu Sans", size = rel(1.2)),
        panel.background = element_rect(fill = "white"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 10),
        legend.position = "right")  +
  scale_x_continuous(breaks= seq(1,15,by=2), name= "Longitude") + 
  scale_y_continuous(breaks= seq(4,14,by=2), name="Latitude") + 
  coord_map(xlim = c(1,15), ylim = c(4,14))

ggsave("Posterior//2017_MaximumMaizeYield.jpg",width=8, height=4,units="in")

#Create rasters from the data 

head(site.dates.combs.00)
maize.mean.d <- site.dates.combs.00[ , c("Long", "Lat", "max_yield")]
head(maize.mean.d)

xy <- maize.mean.d[,1:2]
ng.mz.dssat.pnts = SpatialPointsDataFrame(maize.mean.d, xy)
ng.mz.dssat.pnts = SpatialPoints(ng.mz.dssat.pnts[,c("Lat","Long")],proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
ng.mz.dssat.pnts.ras <- rasterFromXYZ(as.data.frame(ng.mz.dssat.pnts)[, c("Long", "Lat", "max_yield")])
projection(ng.mz.dssat.pnts.ras) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#ng.mz.dssat.pnts.ras.res <- resample(ng.mz.dssat.pnts.ras, ng.mz.apsim.pnts.ras)
#ng.mz.dssat.pnts.ras.x <- mask(ng.mz.dssat.pnts.ras.res, ng.bound)
2777

ng.mz.dssat.pnts.ras.std <- (ng.mz.dssat.pnts.ras/ng.mz.dssat.pnts.ras@data@max) * 2777 #standardise the yield according to observed yields

plot(ng.mz.dssat.pnts.ras.std, main = "Maize yield")

plot(admin, add = T)
plot(ng.grid, add=T)

#writeRaster(ng.mz.dssat.pnts.ras, "OUT/Ras/curr_yield_w5e5.tif", format="GTiff", overwrite=TRUE)

#Read spam data
spam17.maize.all <- raster("~/DATABANK/SPAM/spam2017v2r1_ssa_yield.geotiff/spam2017V2r1_SSA_Y_MAIZ_R.tif")
ng.spam17.maize <- crop(spam17.maize.all,ng.mz.dssat.pnts.ras)
plot(ng.spam17.maize)
plot(admin, add = T)
plot(ng.grid, add=T)

#resample to 25
res.factor = res(ng.mz.dssat.pnts.ras) / res(ng.spam17.maize)
obs.maize.spam17.res.mean <- aggregate(ng.spam17.maize, fact=res.factor[1], expand=FALSE,fun=mean,na.rm=T)
#obs.maize.spam17.res.mode <- aggregate(ng.spam1717.maize, fact=res.factor[1], expand=FALSE,fun=modal, na.rm=T)
#obs.maize.spam17.res.median <- aggregate(ng.spam1717.maize, fact=res.factor[1], expand=FALSE,fun=median,na.rm=T)

obs.maize.spam17.res <- raster::resample(obs.maize.spam17.res.mean, ng.mz.dssat.pnts.ras)
obs.maize.spam17.res.x <- mask(obs.maize.spam17.res, admin)
plot(obs.maize.spam17.res.x)
plot(admin, add=T)
plot(ng.grid, add=T)

dim(ng.mz.dssat.pnts.ras)
dim(obs.maize.spam17.res.x)

sim17.obs.all <- stack (ng.mz.dssat.pnts.ras.std, obs.maize.spam17.res.x)

coul <- colorRampPalette(brewer.pal(12, "RdYlGn"))
cpal <- c('chocolate4','tomato', 'tan2','gold2','palegreen','darkseagreen', 'forestgreen')

cols <- colorRampPalette(brewer.pal(9,"YlGn"))
#jpeg("~/AGRICA/Uganda/Output/dssat_mapspam17.jpg",width = 6000, height = 4000, res=600)
p <- levelplot(sim17.obs.all, col.regions=rcolors$MPL_RdBu, att='class', ylab=list(label = "Latitude", vjust = -0.3),
               xlab=list(label = "Longitude", vjust = -0.5), layout = c(2, 1),
               colorkey=list(labels=list(cex=0.9), space="right", width= 0.8,height = 0.5), 
               main=list('Maize yield in Nigeria',side=1,line=0.5), axes = FALSE, margin = F,
               names.attr=c("DSSAT Yield", "SPAM Yield")) +
  latticeExtra::layer(sp.polygons(ng.grid, col="white", alpha=0.6, lwd= 2.5)) + 
  latticeExtra::layer(sp.polygons(ng.hollow, fill='white')) +
  latticeExtra::layer(sp.polygons(admin, lwd= 0.6, col="darkorange4", fill='transparent')) + 
  latticeExtra::layer(sp.polygons(agrozones, fill='transparent', col="grey20", alpha=0.6, lwd= 0.4)) +
  #latticeExtra::layer(sp.polygons(water, fill='cyan', col="cyan", alpha=0.6, lwd= 0.4)) + 
  latticeExtra::layer(sp.text(coordinates(agrozones), txt = agrozones$AEZ2, cex= 0.5, pos = 1)) +
  latticeExtra::layer(grid.text('N', x=0.1, y=0.9, rot=0, gp = gpar(fontfamily="Sans", fontface="italic", cex=0.8)))
p + latticeExtra::layer(grid.text(expression(Delta), x=0.1, y=0.85, rot=0, gp = gpar(fontfamily="Times New Roman", cex=1.2)))
#dev.off()

ng.mz.dssat.spam17 <- as.data.frame(raster::extract(obs.maize.spam17.res.x, ng.mz.dssat.pnts))
colnames(ng.mz.dssat.spam17) <- "spam17"
head(ng.mz.dssat.spam17)

ng.mz.dssat.ext <- as.data.frame(raster::extract(ng.mz.dssat.pnts.ras.std, ng.mz.dssat.pnts))
colnames(ng.mz.dssat.ext) <- "dssat17"
head(ng.mz.dssat.ext)


spam17.dssat <- cbind(ng.mz.dssat.ext, ng.mz.dssat.spam17)
spam17.dssat$spam17 <- round(spam17.dssat$spam17,0)
#spam17.dssat.cln <- spam17.dssat[spam17.dssat$spam17 > 10, ]
spam17.dssat.cln <- na.omit(spam17.dssat)
head(spam17.dssat.cln)

#spam17.dssat.cln$max_tfrm <- ((spam17.dssat.cln$max_yield/max(spam17.dssat.cln$max_yield)) * max(spam17.dssat.cln$spam17))
#spam17.dssat.cln <- spam17.dssat.cln[spam17.dssat.cln$max_tfrm > 10, ]
head(spam17.dssat.cln)
dim(spam17.dssat.cln)

library(ggpmisc)
#jpeg("~/DATABANK/PDSSAT/Uganda/OUT/dssat_mapspam17_cor.jpg",width = 6000, height = 4000, res=600)
ggplot(spam17.dssat.cln, aes(x = dssat17, y = spam17)) +
  #stat_poly_line() + geom_text(label=spam17.dssat.cln$site, size=2, color="darkred") +
  #stat_poly_eq() +
  geom_point() + 
  xlab("Simulated yield (Kg/ha)") + ylab("spam17 yield (kg/ha)")+
  #ylim(c(0, 7050)) + xlim(c(0, 7050)) +  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, 
                                    linetype="solid"),
        panel.grid.minor = element_line(color = "steelblue4", size = 0.5,linetype = "dashed"))
#dev.off()
#write.csv(spam.dssat, "OUT/spam_dssat_ug_l5.csv")

library(hydroGOF)
gof17.dssat <- gof(sim=spam17.dssat.cln$dssat17, obs=spam17.dssat.cln$spam17) 
gof17.dssat

gof17.dssat.name <- cbind(row.names(gof17.dssat), gof17.dssat)
gof17.dssat.name                        

#Evalauate per AEZ
all.coords.ext <- as.data.frame(ng.mz.dssat.pnts@coords)
all.dssat.spam.site.pt <- cbind(all.coords.ext[,1:2], spam17.dssat)
head(all.dssat.spam.site.pt)

coordinates(all.dssat.spam.site.pt) = ~Long + Lat 
projection(all.dssat.spam.site.pt) <- projection(agrozones)

all.dssat.spam.aez <- sp::over(all.dssat.spam.site.pt, agrozones)
head(all.dssat.spam.aez)
all.dssat.spam.aez.match <- cbind(all.dssat.spam.aez[,c("AEZ")], spam17.dssat)
colnames(all.dssat.spam.aez.match)[1] <- "AEZ"
head(all.dssat.spam.aez.match)


all.dssat.spam.site.all <- all.dssat.spam.site.pt
colnames(all.dssat.spam.site.all) <- c("Long","Lat","Simulated","Reference")
all.dssat.spam.site.long <- reshape2::melt(all.dssat.spam.site.all, id.vars=c("Long", "Lat"))

head(all.dssat.spam.site.long)

#plot with ggplot
ggplot() +
  geom_tile(data = all.dssat.spam.site.long, aes(x = Long, y = Lat, fill = value)) + 
  scale_fill_gradientn("Yield\n[kg/ha]", limits=c(0,4100), values = c(0, 0.3,0.6, 1),  
                       colours=rev(rcolors$BlueYellowRed),
                       na.value = 'grey') +
  geom_polygon(data=admin, aes(x = long, y=lat, group=group), fill = NA, col="black", size = 1.5) + 
  geom_polygon(data=ng.grid, aes(x = long, y=lat, group=group), fill = NA, col="grey82", size = 0.2) +
  geom_polygon(data=agrozones, aes(x = long, y=lat, group=group), fill = NA, linetype="dashed", col="grey42", size = 0.2) +
  geom_polypath(data=ng.hollow, aes(x = long, y=lat, group=group), fill = "white", col="grey62", size = 0.2) + # for the country map
  geom_text(data = centroids.df, aes(lon, lat, label = label), col="salmon3", size = 2.5) +
  geom_hline(yintercept=seq(4,14,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  geom_vline(xintercept=seq(1,15,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  theme(axis.text=element_text(family="DejaVu Sans", size = rel(0.8)),
        axis.title.x=element_text(margin=margin(t=0.4, unit="cm"), size=11,family="DejaVu Sans"),
        axis.title.y=element_text(margin=margin(r=0.4, unit="cm"), size=11, family="DejaVu Sans"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title =element_text(size = 12, face = 'plain', family="DejaVu Sans", margin=margin(b=0.4, unit="cm")),
        legend.key.width = unit(0.35, "cm"), 
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(family="DejaVu Sans", size = rel(1.2)),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 10),
        legend.position = "right")  +
  facet_grid(~variable) + 
    scale_x_continuous(breaks= seq(1,15,by=2), name= "Longitude") + 
  scale_y_continuous(breaks= seq(4,14,by=2), name="Latitude") + 
  coord_map(xlim = c(1,15), ylim = c(4,14))

ggsave("Figures/simulated_vs_spam.jpg",width=10, height=4,units="in")

#"Derived Savanna", "Humid Forest", "Southern Guinea Savanna", "Northern Guinea Savanna" 
#"Semi-arid/Sudan Savanna", Arid/Sahel", "Mid Altitude", "Water bodies"           

#For each zone
arid.sahel.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Arid/Sahel")
der.savana.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Derived Savanna")
humid.for.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Humid Forest")
mid.altit.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Mid Altitude")
north.gun.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Northern Guinea Savanna")
semi.arid.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Semi-arid/Sudan Savanna")
south.gun.acc <- all.dssat.spam.aez.match %>% filter(AEZ == "Southern Guinea Savanna")
#GOF
arid.sahel.gof <- gof(sim=arid.sahel.acc$dssat17, obs=arid.sahel.acc$spam17) 
arid.sahel.gof
der.savana.gof <- gof(sim=der.savana.acc$dssat17, obs=der.savana.acc$spam17) 
der.savana.gof
humid.for.gof <- gof(sim=humid.for.acc$dssat17, obs=humid.for.acc$spam17)
humid.for.gof
mid.altit.gof <- gof(sim=mid.altit.acc$dssat17, obs=mid.altit.acc$spam17) 
mid.altit.gof
north.gun.gof <- gof(sim=north.gun.acc$dssat17, obs=north.gun.acc$spam17)
north.gun.gof
semi.arid.gof <- gof(sim=semi.arid.acc$dssat17, obs=semi.arid.acc$spam17)
semi.arid.gof
south.gun.gof <- gof(sim=south.gun.acc$dssat17, obs=south.gun.acc$spam17) 
south.gun.gof

write.csv(spam17.dssat.cln, "Posterior/AEZ/GOF/spam17_dssat_cln.csv")
write.csv(arid.sahel.acc, "Posterior/AEZ/GOF/arid_sahel_acc.csv")
write.csv(der.savana.acc, "Posterior/AEZ/GOF/der_savana_acc.csv")
write.csv(humid.for.acc, "Posterior/AEZ/GOF/humid_for_acc.csv")
write.csv(mid.altit.acc, "Posterior/AEZ/GOF/mid_altit_acc.csv")
write.csv(north.gun.acc, "Posterior/AEZ/GOF/north_gun_acc.csv")
write.csv(semi.arid.acc, "Posterior/AEZ/GOF/semi_arid_acc.csv")
write.csv(south.gun.acc, "Posterior/AEZ/GOF/south_gun_acc.csv")

###-----------------------------------------------------------------------------------------------------
print("Inini Ndinonzi Inini")


