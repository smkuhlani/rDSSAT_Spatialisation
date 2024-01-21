head(site.dates.combs)

a1.dec <- site.dates.combs[site.dates.combs$year <= 1990,]
a2.dec <- site.dates.combs[site.dates.combs$year > 1990 & site.dates.combs$year <= 2000, ]
a3.dec <- site.dates.combs[site.dates.combs$year >  2000 & site.dates.combs$year >= 2010, ]
a4.dec <- site.dates.combs[site.dates.combs$year >  2010 & site.dates.combs$year >= 2015, ]

head(a1.dec)
#Mode function
# Mode1 <- function(x) {
#   ux <- unique(x)
#   if(!anyDuplicated(x)){
#     NA_character_ } else { 
#       tbl <-   tabulate(match(x, ux))
#       toString(ux[tbl==max(tbl)])
#     }
# }
# 
# a1.dec.mode <- a1.dec %>%
#   group_by(site) %>%
#   summarise(mode_pdat = Mode1(bpdat))
# head(a1.dec.mode)

Mode2 <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

a1.dec.mode <- as.data.frame(a1.dec %>%
                               group_by(site) %>%
                               summarise(a1_mode = Mode2(bpdat)))
head(a1.dec.mode)

a2.dec.mode <- as.data.frame(a2.dec %>%
                               group_by(site) %>%
                               summarise(a2_mode = Mode2(bpdat)))
head(a2.dec.mode)

a3.dec.mode <- as.data.frame(a3.dec %>%
                               group_by(site) %>%
                               summarise(a3_mode = Mode2(bpdat)))
head(a3.dec.mode)

a4.dec.mode <- as.data.frame(a4.dec %>%
                               group_by(site) %>%
                               summarise(a4_mode = Mode2(bpdat)))
head(a4.dec.mode)

modal.all2 <- merge(a1.dec.mode, a2.dec.mode, by= "site")
modal.all3 <- merge(modal.all2, a3.dec.mode, by= "site")
modal.all4 <- merge(modal.all3, a4.dec.mode, by= "site")
head(modal.all4)

modal.all4$a1a4 <- modal.all4$a4_mode - modal.all4$a1_mode
modal.all4$a1a3 <- modal.all4$a3_mode - modal.all4$a1_mode
modal.all4$a1a2 <- modal.all4$a2_mode - modal.all4$a1_mode

modal.a1a4 <- modal.all4[,c("site","a1a4")]
modal.a1a4$dec <- rep("a1a4", nrow(modal.a1a4))
colnames(modal.a1a4) <- c("site","delt" ,"dec")
modal.a1a4.site <- merge(modal.a1a4, loc.data.cut, by= "site")

modal.a1a3 <- modal.all4[,c("site","a1a3")]
modal.a1a3$dec <- rep("a1a3", nrow(modal.a1a3))
colnames(modal.a1a3) <- c("site","delt" ,"dec")
modal.a1a3.site <- merge(modal.a1a3, loc.data.cut, by= "site")

modal.a1a2 <- modal.all4[,c("site","a1a2")]
modal.a1a2$dec <- rep("a1a2", nrow(modal.a1a2))
colnames(modal.a1a2) <- c("site","delt" ,"dec")
modal.a1a2.site <- merge(modal.a1a2, loc.data.cut, by= "site")

modal.site <- rbind(modal.a1a4.site, modal.a1a3.site, modal.a1a2.site)
head(modal.site)

#Plot
ggplot() +
  geom_tile(data = modal.site, aes(x = Long, y = Lat, fill = delt)) + 
  scale_fill_gradientn("Change in date\n [Days]",  limits=c(-150,150), colours=c(brewer.pal(9,"YlGn")),
                       na.value = 'grey')+
  geom_polygon(data=admin, aes(x = long, y=lat, group=group), fill = NA, col="black", size = 0.8) +
  geom_polygon(data=region, aes(x = long, y=lat, group=group), fill = NA, col="grey42", size = 0.2) + 
  #geom_text(data = centroids.df, aes(lon, lat, label = label), col="darkred", size = 3) +
  geom_hline(yintercept=seq(4,14,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  geom_vline(xintercept=seq(1,15,by=2), linetype="dashed", lwd=0.3, color = "grey62") + 
  theme(axis.text=element_text(family="Calibri", size = rel(1.2)), 
        axis.title.x=element_text(margin=margin(t=0.4, unit="cm"), size=7,family="Calibri"),
        axis.title.y=element_text(margin=margin(r=0.4, unit="cm"), size=7, family="Calibri"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title =element_text(size = 13, face = 'bold', family="Calibri", margin=margin(b=0.4, unit="cm")),
        legend.key.width = unit(0.35, "cm"), 
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(family="Calibri", size = rel(1.2)),
        legend.position = "right")  +
  # geom_polygon(data=Ng.grid, aes(x = long, y=lat, group=group), fill = NA, col="grey80", size = 0.05) + 
  facet_grid(~dec) +
  scale_x_continuous(breaks= seq(1,15,by=2), name= "Longitude (?E)") + 
  scale_y_continuous(breaks= seq(4,14,by=2), name="Latitude (?N)") + 
  coord_map(xlim = c(2,15), ylim = c(4,14)) 
