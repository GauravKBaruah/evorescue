
rm(list=ls())
# R script for figure 4 parameter space
#libraries needed for this R script
source('01_functions.R')
require(deSolve)
require(cowplot) 
library(dplyr)
library(readr)
library(beepr)
library(ggplot2)
library(viridis)
library(igraph)
library(bipartite)
library(akima)
library(sna)
library(network)
library(GGally)


#data for figure 3
load("evo_rescue_parameter_space.RData")

#for the first network
web1<- fact1 %>% filter(webname == "datasets_1/M_PL_052.csv", h2=="evo_trait_and_variance")


w1.hvar<-web1 

#offsetting temperature to temperature shift by difference with original 22 degree C
w1.hvar$temperature_shift<-w1.hvar$temperature_shift -22
temp<-with(w1.hvar,interp(w1.hvar$temperature_shift,w1.hvar$w_c ,
                          w1.hvar$richness,xo=seq(min(w1.hvar$temperature_shift),max(w1.hvar$temperature_shift),length=9), 
                          yo=seq(min(w1.hvar$w_c),max(w1.hvar$w_c), length=10), duplicate ='mean'))


temp<-interp2xyz(temp, data.frame=TRUE)
colnames(temp)<-c("x","y","Proportion_survived")
w1<-ggplot(temp , aes(x=x,y=y,z=round(Proportion_survived,3)))+
  geom_raster(aes(fill=Proportion_survived),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0.1,0.92), high = 'yellow', low = 'red')+
  theme_cowplot()+ ylim(c(0.05,0.51))+ylab(expression(omega[c]))+labs(x = "Temperature shift",
                                           color = "Proportion survived") + ggtitle("Evolution of trait and variance")
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
w1



# Same netwokr but now only evolution of mean trait
web2<- fact1 %>% filter(webname == "datasets_1/M_PL_052.csv", h2=="evo_trait")
w2<-web2 

w2$temperature_shift<-w2$temperature_shift -22

temp2<-with(w2,interp(w2$temperature_shift,w2$w_c ,
                      w2$richness,xo=seq(min(w2$temperature_shift),max(w2$temperature_shift),length=9), 
                          yo=seq(min(w2$w_c),max(w2$w_c), length=10), duplicate ='mean'))


temp2<-interp2xyz(temp2, data.frame=TRUE)
colnames(temp2)<-c("x","y","Proportion_survived")
w2<-ggplot(temp2 , aes(x=x,y=y,z=round(Proportion_survived,3)))+
  geom_raster(aes(fill=Proportion_survived),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0.1,0.92), high = 'yellow', low = 'red')+
  
  xlab("Temperature shift")+ylab(expression(omega[c])) + ggtitle("Evolution of mean trait")+
  theme_cowplot()+ ylim(c(0.05,0.51))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
w2



#same network but without any kind of evolution
web3<- fact1 %>% filter(webname == "datasets_1/M_PL_052.csv", h2=="no_evo")
w3<-web3 

w3$temperature_shift<-w3$temperature_shift -22

temp3<-with(w3,interp(w3$temperature_shift,w3$w_c ,
                      w3$richness,xo=seq(min(w3$temperature_shift),max(w3$temperature_shift),length=9), 
                      yo=seq(min(w3$w_c),max(w3$w_c), length=10), duplicate ='mean'))


temp3<-interp2xyz(temp3, data.frame=TRUE)
colnames(temp3)<-c("x","y","Proportion_survived")
w3<-ggplot(temp3 , aes(x=x,y=y,z=round(Proportion_survived,3)))+
  geom_raster(aes(fill=Proportion_survived),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0.1,0.92), high = 'yellow', low = 'red')+
  
  xlab("Temperature shift")+ylab(expression(omega[c])) + ggtitle("No evolutionary dynamics")+
  theme_cowplot()+ ylim(c(0.05,0.51))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
w3




########################################## 2nd network #################################################


#with evolution of trait variance and mean trait
web4<- fact1 %>% filter(webname == "datasets_1/M_PL_003.csv", h2=="evo_trait_and_variance")


w4.hvar<-web4 

w4.hvar$temperature_shift<-w4.hvar$temperature_shift -22
temp4<-with(w4.hvar,interp(w4.hvar$temperature_shift,w4.hvar$w_c ,
                          w4.hvar$richness,xo=seq(min(w4.hvar$temperature_shift),max(w4.hvar$temperature_shift),length=9), 
                          yo=seq(min(w4.hvar$w_c),max(w4.hvar$w_c), length=10), duplicate ='mean'))


temp4<-interp2xyz(temp4, data.frame=TRUE)
colnames(temp4)<-c("x","y","Proportion_survived")
w4<-ggplot(temp4 , aes(x=x,y=y,z=round(Proportion_survived,3)))+
  geom_raster(aes(fill=Proportion_survived),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0.1,0.92), high = 'yellow', low = 'red')+
  theme_cowplot()+ ylim(c(0.05,0.51))+ylab(expression(omega[c]))+labs(x = "Temperature shift",
                                           color = "Proportion survived") + 
  ggtitle("Evolution of trait and variance")
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
w4


#only evolution of mean trait
web2<- fact1 %>% filter(webname == "datasets_1/M_PL_003.csv", h2=="evo_trait")
w5<-web2 

w5$temperature_shift<-w5$temperature_shift -22

temp5<-with(w5,interp(w5$temperature_shift,w5$w_c ,
                      w5$richness,xo=seq(min(w5$temperature_shift),max(w5$temperature_shift),length=9), 
                      yo=seq(min(w5$w_c),max(w5$w_c), length=10), duplicate ='mean'))


temp5<-interp2xyz(temp5, data.frame=TRUE)
colnames(temp5)<-c("x","y","Proportion_survived")
w5<-ggplot(temp5 , aes(x=x,y=y,z=round(Proportion_survived,3)))+
  geom_raster(aes(fill=Proportion_survived),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0.1,0.92), high = 'yellow', low = 'red')+
  
  xlab("Temperature shift")+ylab(expression(omega[c])) + ggtitle("Evolution of mean trait")+
  theme_cowplot()+ ylim(c(0.05,0.51))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
w5



#no evolution at all
web3<- fact1 %>% filter(webname == "datasets_1/M_PL_052.csv", h2=="no_evo")
w6<-web3 

w6$temperature_shift<-w6$temperature_shift -22

temp6<-with(w6,interp(w6$temperature_shift,w6$w_c ,
                      w6$richness,xo=seq(min(w6$temperature_shift),max(w6$temperature_shift),length=9), 
                      yo=seq(min(w6$w_c),max(w6$w_c), length=10), duplicate ='mean'))


temp6<-interp2xyz(temp6, data.frame=TRUE)
colnames(temp6)<-c("x","y","Proportion_survived")
w6<-ggplot(temp6 , aes(x=x,y=y,z=round(Proportion_survived,3)))+
  geom_raster(aes(fill=Proportion_survived),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0.1,0.92), high = 'yellow', low = 'red')+
  
  xlab("Temperature shift")+ylab(expression(omega[c])) + ggtitle("No evolutionary dynamics")+
  theme_cowplot()+ ylim(c(0.05,0.51))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
w6









################################# plotting networks #######################################


mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]




g<-adj.mat("datasets_1/M_PL_003.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}
net = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
ggnet2(net,  mode="circle",  color ="color", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)


deg<-c(degree.plants,degree.animals)
net_2<-ggnet2(net, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")







###################################### network 2 #############################3




g<-adj.mat("datasets_1/M_PL_052.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}
net = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
ggnet2(net,  mode="circle",  color ="color", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)


deg<-c(degree.plants,degree.animals)
net_1<-ggnet2(net, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")


net_1

ggpubr::ggarrange(net_1,w1,w2,w3,
                  net_2,w4,w5,w6,ncol=4,nrow=2)









