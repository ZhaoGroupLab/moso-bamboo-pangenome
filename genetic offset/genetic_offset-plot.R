library(tidyr)
library(ggplot2)
library(raster)
library("sf")
library(rgdal)
library(maptools)
library(gridExtra)
library(scales)
library(gstat)
library(sp)

offest_color= c("#2892C7","#57A0BA","#78ADAC","#97BD9E","#B5CF8F",
                "#CFDB8A","#E8DE82","#f7cb79","#F7B76D","#f59e5f",
                "#F58653","#F7754D","#FF3333")
names(offest_color)=c("0-0.025","0.025-0.0275","0.0275-0.03","0.03-0.0325",
                      "0.0325-0.035","0.035-0.040","0.040-0.045","0.045-0.050",
                      "0.050-0.055","0.055-0.060","0.060-0.065","0.065-0.070","0.070+")
country_shp=sf::st_read("D:/offset_serious/7-Genetic_offset/gis/行政边界/china/1.shp")
raster_ACCESS126=stack("future_tif/126/wc2.1_2.5m_bioc_ACCESS-CM2_ssp126_2061-2080.tif")
raster_ACCESS585=stack("future_tif/585/wc2.1_2.5m_bioc_ACCESS-CM2_ssp585_2061-2080.tif")
creategroup <- function(gf){
  colnames(gf)=c("x","y","offset")
  gf$level=ifelse(gf$offset<=0.025,"0-0.025",
                  ifelse(gf$offset<=0.0275,"0.025-0.0275",
                         ifelse(gf$offset<=0.03,"0.0275-0.03",
                                ifelse(gf$offset<=0.0325,"0.03-0.0325",
                                       ifelse(gf$offset<=0.035,"0.0325-0.035",
                                              ifelse(gf$offset<=0.04,"0.035-0.040",
                                                     ifelse(gf$offset<=0.045,"0.040-0.045",
                                                            ifelse(gf$offset<=0.05,"0.045-0.050",
                                                                   ifelse(gf$offset<=0.055,"0.050-0.055",
                                                                          ifelse(gf$offset<=0.06,"0.055-0.060",
                                                                                 ifelse(gf$offset<=0.065,"0.060-0.065",
                                                                                        ifelse(gf$offset<=0.07,"0.065-0.070","0.070+"
                                                                                        ))))))))))))
  return(gf)
}
#### ACCESS-CM2_ssp126_2061-2080
resdf<-as.data.frame(raster_ACCESS126,xy=TRUE)%>%drop_na()
resdf=creategroup(resdf)
p_8026=ggplot()+
  #geom_sf(fill="#DADADA",data=country_shp)+
  geom_raster(aes(x=x,y=y,fill=factor(level)),data=resdf,stat="identity")+
  geom_sf(fill="transparent",data=country_shp,size=0.5)+
  coord_sf(xlim=c(115, 135),ylim=c(39.4, 54.7))+
  scale_fill_manual(values=offest_color[unique(resdf$level)])+
  scale_x_continuous(limits = c(115, 135))+
  scale_y_continuous(limits = c(39.4, 54.7))+
  labs(x="Longitude",y="Latitude") +
  theme_bw()+
  theme(text=element_text(family="serif"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key.size=unit(0.1,'inch'),
        legend.title = element_text(size=10.5, color="black",vjust=0.5, hjust=0.5),
        legend.position = "none",
        legend.background=element_rect(colour= "grey" ,fill= "white" ,size=0.6),
        legend.text= element_text(size=7.3, color="black",vjust=0.5, hjust=0.5),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(colour = "white",size=0.1,linetype = 4),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))+
  guides(fill=F)
ggsave(p_8026,file="2-ACCESS126.pdf",height=8,width=8)

#### ACCESS-CM2_ssp585_2061-2080
resdf<-as.data.frame(raster_ACCESS585,xy=TRUE)%>%drop_na()
resdf=creategroup(resdf)
p_8070=ggplot()+
  #geom_sf(fill="#DADADA",data=country_shp)+
  geom_raster(aes(x=x,y=y,fill=factor(level)),data=resdf,stat="identity")+
  geom_sf(fill="transparent",data=country_shp,size=0.5)+
  coord_sf(xlim=c(115, 135),ylim=c(39.4, 54.7))+
  scale_fill_manual(values=offest_color[unique(resdf$level)])+
  scale_x_continuous(limits = c(115, 135))+
  scale_y_continuous(limits = c(39.4, 54.7))+
  labs(x="Longitude",y="Latitude") +
  theme_bw()+
  theme(text=element_text(family="serif"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key.size=unit(0.1,'inch'),
        legend.title = element_text(size=10.5, color="black",vjust=0.5, hjust=0.5),
        legend.position = "none",
        legend.background=element_rect(colour= "grey" ,fill= "white" ,size=0.6),
        legend.text= element_text(size=7.3, color="black",vjust=0.5, hjust=0.5),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(colour = "white",size=0.1,linetype = 4),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))+
  guides(fill=F)
ggsave(p_8070,file="2-ACCESS585.pdf",height=8,width=8)