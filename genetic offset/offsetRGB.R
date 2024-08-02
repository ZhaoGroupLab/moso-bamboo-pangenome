######offset
library(sf)
library(ggplot2)

pops <- read.csv("popInfo.csv")

clipped_shapefile <- st_read("D:/offset_serious/7-Genetic_offset/gis/行政边界/研究区域_3.shp")


#forward_offset=read.csv('D:/offset_serious/7-Genetic_offset/gis/0926/126/ACCESS-CM2_ssp126_2061-2080_forward.csv')
#reverse_offset=read.csv('D:/offset_serious/7-Genetic_offset/gis/0926/126/ACCESS-CM2_ssp126_2061-2080_reverse.csv')

creategroup <- function(tiff){
  colnames(tiff)=c("x","y","bio_value")
  var=max(tiff$bio_value)-min(tiff$bio_value)
  min=min(tiff$bio_value)
  tiff$level=
    ifelse(tiff$bio_value>=min+var/20*19,'255',
           ifelse(tiff$bio_value>=min+var/20*18,'253',
                  ifelse(tiff$bio_value>=min+var/20*17,'250',
                         ifelse(tiff$bio_value>=min+var/20*16,'247',
                                ifelse(tiff$bio_value>=min+var/20*15,'244',
                                       ifelse(tiff$bio_value>=min+var/20*14,'241',
                                              ifelse(tiff$bio_value>=min+var/20*13,'237',
                                                     ifelse(tiff$bio_value>=min+var/20*12,'235',
                                                            ifelse(tiff$bio_value>=min+var/20*11,'230',
                                                                   ifelse(tiff$bio_value>=min+var/20*10,'225',
                                                                          ifelse(tiff$bio_value>=min+var/20*9,'220',
                                                                                 ifelse(tiff$bio_value>=min+var/20*8,'215',
                                                                                        ifelse(tiff$bio_value>=min+var/20*7,'210',
                                                                                               ifelse(tiff$bio_value>=min+var/20*6,'205',
                                                                                                      ifelse(tiff$bio_value>=min+var/20*5,'200',
                                                                                                             ifelse(tiff$bio_value>=min+var/20*4,'180',
                                                                                                                    ifelse(tiff$bio_value>=min+var/20*3,'120',
                                                                                                                           ifelse(tiff$bio_value>=min+var/20*2,'80','0'
                                                                                                                           ))))))))))))))))))
  return(tiff)
}

#data=merge(forward_offset[c("x1",'y1','local','forwardOffset')],reverse_offset[c("x1",'y1','reverseOffset')],by=c('x1','y1'))


data = data126
local=creategroup(data[c("x1",'y1','local')])
colnames(local)=c('x','y','local','local_red')
forward=creategroup(data[c("x1",'y1','forwardOffset')])
colnames(forward)=c('x','y','forward','forward_green')
reverse=creategroup(data[c("x1",'y1','reverseOffset')])
colnames(reverse)=c('x','y','reverse','reverse_blue')
temp=merge(local,forward,by=c('x','y'))
data=merge(temp,reverse,by=c('x','y'))
a=c()
for (i in 1:dim(data)[1]){
  color=rgb(as.numeric(data$local_red[i]),as.numeric(data$forward_green[i]),as.numeric(data$reverse_blue[i]),maxColorValue = 255 ) 
  a<-append(a,color)
}
data$color=a
color=levels(factor(data$color))
write.csv(data,"126RGB.csv")

pdf("126offset.pdf",width=6,height=6)
ggplot()+
  geom_tile(aes(x = x, y = y, fill = factor(color)), data = data) +
  geom_point(aes(x = long, y = lat), pops, fill = 'white', shape = 21, size = 3) +
  geom_sf(data = clipped_shapefile_1, color = "black", fill = "transparent", size = 1.5) +
  scale_fill_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()
############local
pdf('126local_forward_scatter.pdf',width=5,height=5)
ggplot(data, aes(x=local, y=forward,color=factor(color)))  + geom_point(size=1.5) +
  xlim(0,max(data$local))+ylim(0,max(data$forward))+
  scale_color_manual(values=color)+geom_abline(mapping=NULL,45)+
  theme_bw()+
  theme(legend.position = "none")
dev.off()
pdf('126local_reverse_scatter.pdf',width=5,height=5)
ggplot(data, aes(x=local, y=reverse,color=factor(color)))  + 
  xlim(0,max(data$local))+ylim(0,max(data$reverse))+
  geom_abline(mapping=NULL,45)+
  geom_point(size=1.5) +
  scale_color_manual(values=color)+
  theme_bw()+
  theme(legend.position = "none")
dev.off()
pdf('126forward_reverse_scatter.pdf',width=5,height=5)
ggplot(data, aes(x=forward, y=reverse,color=factor(color)))  + geom_point(size=1.5) +
  xlim(0,max(data$forward))+ylim(0,max(data$reverse))+
  scale_color_manual(values=color)+geom_abline(mapping=NULL,45)+
  theme_bw()+
  theme(legend.position = "none")
dev.off()