#                                 WORKSHOP - BIOINFIRMATICS
#                           Silvia Fajardo and Andressa Dias Leao


# Step 1
# We installed the packages before starting 
install.packages("ggplot2")
install.packages("gridExtra")
#setwd is in the results file
library("ggplot2")
library("gridExtra")
# set the directory, it is saved on the data file
setwd("~/Desktop/WS-repeat/Data") 

# Step 2
#upload satelite data
chr="Chr2"
satelite<-read.table("~/Desktop/WS-repeat/Data/X.tropicalis_Chr_2.fa.2.5.7.80.10.50.200.dat.nr.tab.map",header=FALSE)
satelite
map1<-setNames(satelite,c("chr","PosX","PosY","sat_span","sat_number","sat_type"))
# select the data corresponding to the chromosome 2
map2<-map1[map1$chr==chr,]
#to check if everythink is ok
head(map1)
head(map2)


#Step 3 
#Set up the plot. label

xlab <-"Xt Chr. 2 coordinate in bp"
ylab <-"Nb of satellites per 150 kbp"

#Step 3.1 to build the lines, we used the Line plot function

ggplot(map1)+geom_line(aes(y=sat_number,x=PosX,color=sat_type),size=1)+ylab("Nb of satellites per 150 kbp")+xlab("Xt Chr. 2 coordinate in bp")+ggtitle("satellite density on X.tropicalis chromosome 2")
#ggsave(paste0("Figures/XT.",chr,".density.pdf"),width=16,height=8)

#Step 3.2 to analyze separately the Nb of satellites we used the Facet function

plot_trf_dat<-ggplot(map1)+geom_line(aes(y= sat_number,x= PosX,color=sat_type),size=1)+ylab("Nb of satellites per 150 kbp")+xlab("Xt Chr. 2 coordinate in bp")+ggtitle("satellite density on X.tropicalis chromosome 2")
plot_trf_dat+facet_wrap(~sat_type)

# Step 3.3 Best pick
#a- The base plot is made to see the area under each line
ggplot(map1,aes(y= sat_number,x= PosX,fill=sat_type))+geom_area()

#to change the name, we define this plot as being plot0
plot0<-ggplot(map1,aes(y= sat_number,x= PosX,fill=sat_type))+geom_area()

#b- We change the appearance of the drawing to increase readability. 
#This makes the plot look more defined
ggplot(map1,aes(y= sat_number,x= PosX,fill=sat_type))+geom_area(colour="black",size=0.2)
#c- to differentiate from plot0, we define this plot as being plot1
plot1<-ggplot(map1,aes(y= sat_number,x= PosX,fill=sat_type))+geom_area(colour="black",size=0.2)
#and then we have to check plot1
plot1

#d- We add names for the axis x and y and changed the title
plot1+ylab("Nb of satellites per 150 kbp")+xlab("Chr. 2 XENOPUS coordinate in bp")+ggtitle("SATELLITE DENSITY ON X. tropicalis CHROMOSOME 2")
plot2<-plot1+ylab("Nb of satellites per 150 kbp")+xlab("Chr. 2 XENOPUS coordinate in bp")+ggtitle("SATELLITE DENSITY ON X. tropicalis CHROMOSOME 2")

#e- color change! We add a personalized grey scheme picked from color brewer and remove the text of the fill legend
plot2+scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0"))
plot3<-plot2+scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0"),guide_legend(title=""))

#shunt the previous step
plot3<-plot2

#f- to make it nicer to read, we change the appearance of text labels: bold fonts and center the title
plot3+theme(axis.text=element_text(face="bold",size=12),axis.title=element_text(face="bold",size=14),title=element_text(face="bold",size=16))
plot4<-plot3+theme(axis.text=element_text(face="bold",size=12),axis.title=element_text(face="bold",size=14),title=element_text(face="bold",size=12),plot.title = element_text(hjust = 0.5))
#g- visual is important so we also decided to modify the location of the legend to increase the space allocated to the plot, and modify the background!!!
plot4+theme(legend.position=c(0.9,0.9),panel.background=element_rect(fill="white",colour="black"),panel.grid = element_line(colour = "grey80",size=0.2))
plot5<-plot4+theme(legend.position=c(0.9,0.9),legend.text=element_text(face="bold"),panel.background=element_rect(fill="white",colour="black"),panel.grid = element_line(colour = "grey80",size=0.2))
#g- We modify the appearance of the x scale to make it more readable
plot5+scale_x_continuous(labels=scales::comma)
plot6<-plot5+scale_x_continuous(labels=scales::comma)

#h- We decided to save the plot as a pdf file
ggsave("Figure_C9_10S_satellite_distribution.pdf",width=16,height=8)

# To see the differences between the base plot and the final plot we used grid arrange. 
grid.arrange(plot0, plot6, ncol=2)


#Step 3.4 also we made a Boxplot to see the density between the satellites types.
head(map1)
ggplot(map1,aes(x=sat_type,y=sat_number,fill=sat_type))+
  geom_boxplot()+scale_y_log10()+xlab("satellite type")+
  ylab("satellite density")+scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0"))

#Step 3.4 Finally we made an histogram to plot the satellite density od each satellite type
ggplot(map1)+geom_histogram(mapping=aes(sat_number,fill=sat_type),color="black")+xlab("satellite type")+xlab("satellite density")+ylab("Count")+scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0"))+theme(axis.text=element_text(face = "bold"))
