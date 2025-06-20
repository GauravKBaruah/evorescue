rm(list=ls())

#R script for figure 2
#### loading the functions script where all the functions are coded.
source('01_functions.R')

#all the important libraries maybe there are others
#which could be loaded in one command. I am lazy to find out.
require(deSolve) 
require(cowplot) 
library(dplyr)
library(readr)
library(beepr)
library(ggplot2)
library(viridis)
library(igraph)
library(bipartite)
library(GGally)
library(network)
library(sna)
library(ggplot2)

set.seed(32)
theme_set(theme_bw()) 

mydir = 'datasets_1'

#-------------------------------------------------------------------------------
webfiles = list.files(path = mydir, pattern = "*.csv", full.names = TRUE)

  g<-adj.mat(webfiles[which(webfiles =="datasets_1/M_PL_061_31.csv")]) #ork web names
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.animals<-degree.plants<-numeric()
  
  #degree of plants and anichmals
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  
  degree<-c(degree.animals, degree.plants)

  #initial state values
  N <- runif( (Aspecies+Plantspecies) , 1,1)  ## initial species densities
  na <-  runif( (Aspecies) , 1,1) #initial densities
  np<- runif( (Plantspecies) , 1,1)
  muA<- runif(Aspecies, 19, 24) #initial mean phenotypic optimum trait values
  muP<- runif(Plantspecies, 19, 24)
  
  #max simulation time
  tmax<-2000
  temp<-22  #local environmental temperature
  nestedness<-nestedness_NODF(g) #nestedness NODF of the network
  C<-Connectance(g) #connectance of the network
  web.name<-"datasets_1/M_PL_061_31.csv" #name of the network
  dganimals<-degree.animals #animal degree
  dgplants<-degree.plants #plant degree
  sig <-runif((Aspecies+Plantspecies),0.009,0.03) #vector of trait variances 
  envA<-runif(Aspecies,0.009, 0.03) #environmental variances -animals
  envP<- runif(Plantspecies, 0.009, 0.03)# environmental variances plants
  sa<- sig[1:Aspecies]
  sp<-sig[(Aspecies+1): (Aspecies+Plantspecies)]
  pa<-envA+sa #phenotypic variance of animals
  pp<-envP+sp #phenotypic variance of plants
  trait_evolve<-1 #evolution of mean trait
  var_evolve<-0  #no evolution of variance 
  w<-1 #width of mutualistic interaction
  w_c<-0.2 #width of species competition
  mut.strength<-1 #avg. mutualistic interaction, gamma_0
  params <- list(time=time,matrix=g,Temp=temp,tmax=tmax,
                 sig=sig,variation="high",
                 web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=var_evolve,
                 mut.strength=mut.strength,C=C,nestedness=nestedness,h2="evo_trait",
                 web.name=web.name,dganimals=dganimals,degree=degree,w_c=w_c,
                 dgplants=dgplants,envA=envA,envP=envP)
  
  t1<-plot_snapshot(Na=na, Np=np, m = c(muA,muP), 
                sigma = c(pa,pp), Temp=22, limits = c(19,25),res = 701)+
    annotate("text",x=22, y =19, label="time == 0",parse=T)
  
  t1

  
  #parameters
  bw  <-2  #width of the temperature tolerance curve
  aw  <- 0 #not-needed
  gi <- 1.5 #peak of the temperature tolerance curve
  ki <-0.4#mortality rate
  ic <-c(na, np, muA, muP, sa , sp )
  species_index <- which(degree==max(degree))
  #parameter list
  params <- list(time=time,matrix=g,Temp=temp,tmax=tmax,species_index=species_index,
                 sig=sig,bw=bw,gi=gi,ki=ki,variation="high",
                 web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=var_evolve,
                 mut.strength=mut.strength,C=C,nestedness=nestedness,h2="evo_trait",
                 web.name=web.name,dganimals=dganimals,degree=degree,w_c=0.2,
                 dgplants=dgplants,envA=envA,envP=envP)
  species_index <- which(degree==max(degree))
  #time for simulation
  tmax<-2000
  #solver
  OUT<-ode(func=eqs_type2_new, y=ic, parms=params,
           times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params) 
  
  OUT %>% plot_density()
  
  #subset at final timpoint
  sol_1<-  OUT %>% filter(time == tmax)
  Na1<- (sol_1 %>% filter(type=="N"))$v
  Np1<-(sol_1 %>% filter(type =="P"))$v
  ma1<-(sol_1 %>% filter(type == "muA"))$v
  mp1<-(sol_1 %>% filter(type == "muP"))$v
  sa1<-(sol_1 %>% filter(type == "sA"))$v
  sp1<-(sol_1 %>% filter(type == "sP"))$v
  
  #ensuring if density < 1e-5 it ensures it goes to 0.05
  Na1[which(Na1<1e-5)]<-0.05
  Np1[which(Np1<1e-5)]<-0.05
  
  #phenotypic variance = genetic variance+environmental variance
  Pa1<-sa1+envA
  Pp1<-sp1+envP
  #snapshot at max time point
  t2<- plot_snapshot(Na=Na1, Np=Np1, m = c(ma1,mp1), 
                     sigma = c(Pa1,Pp1), Temp=22, limits = c(19,26),res = 1001)+
    annotate("text",x=22, y =50, label="time == 2000",parse=T)  
  
t2


########################### Temperature shift - evolution trait ######################


#parameters
bw  <-2  #width of the temperature tolerance curve
aw  <- 0 #not-needed
gi <- 1.5 #peak of the temperature tolerance curve
ki <-0.4#mortality rate
w<- 1
ic <-c(Na1, Np1, ma1, mp1, sa , sp )
mut.strength<-1
species_index <- which(degree==max(degree))

var_evolve<-0 #variance evolve flag
trait_evolve<-1 #trait mean evolve flag
#parameter list
params <- list(time=time,matrix=g,Temp=26,tmax=tmax,species_index=species_index,sig=sig,bw=bw,gi=gi,ki=ki,variation="high",
               web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=var_evolve,
               mut.strength=mut.strength,C=C,nestedness=nestedness,h2="evo_trait",
               web.name=web.name,dganimals=dganimals,degree=degree,w_c=0.2,
               dgplants=dgplants,envA=envA,envP=envP)

species_index <- which(degree==max(degree))

tmax<-2000
#solver dynamics
OUT_evotrait_net1<-ode(func=eqs_type2_new, y=ic, parms=params,
                       times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params) 


#plotting timeseries
plot_density(OUT_evotrait_net1)
OUT$time
OUT_evotrait_net1$time<-OUT_evotrait_net1$time+tmax

new_dat_evo_trait<-rbind(OUT,OUT_evotrait_net1)

r3<-plot_density_timeseries(new_dat_evo_trait)+
  geom_vline(xintercept = 2000, linetype="dashed",size=1, color="grey5")+
  scale_color_viridis_d(alpha = 1)

r3
gvar3<-plot_density_timeseries_genvar(new_dat_evo_trait)+
  geom_vline(xintercept = 2000, linetype="dashed",size=1, color="grey5")+
  scale_color_viridis_d(alpha = 1)


#subsetting simulated data at final timepoint
Na3<- (OUT_evotrait_net1 %>% filter(time == 4e3,type=="N"))$v
Np3<-(OUT_evotrait_net1 %>% filter(time == 4e3,type =="P"))$v
ma3<-(OUT_evotrait_net1 %>% filter(time == 4e3,type == "muA"))$v
mp3<-(OUT_evotrait_net1 %>% filter(time == 4e3,type == "muP"))$v
sa3<-(OUT_evotrait_net1 %>% filter(time == 4e3,type == "sA"))$v
sp3<-(OUT_evotrait_net1 %>% filter(time == 4e3,type == "sP"))$v

#phenotypic variance
Pa3<-sa3+envA 
Pp3<-sp3+envP

#snapshot at tmax
t_evo_trait<- plot_snapshot(Na=Na3, Np=Np3, m = c(ma3,mp3), 
                            sigma = c(Pa1,Pp1), Temp=26, limits = c(19,31),res = 1001)+
  annotate("text",x=26, y =12, label="time == 4000",parse=T)  
t_evo_trait 




#preparing data for plotting of the plant-pollinator network
g<-adj.mat(webfiles[which(webfiles == "datasets_1/M_PL_061_31.csv")])
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and animals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}


#vertex names 
net = network(g, bipartite = T, directed = FALSE)
deg <- c(degree.plants, degree.animals)

# groups and base color
net %v% "groups" <- ifelse(1:network.size(net) <= dim(g)[1], "plants", "animals")
net %v% "color" <- ifelse(net %v% "groups" == "plants", "#E69F00", "#E69F00")


# Choosing  nodes by index which has densities greater than 1e-4.
highlight_index <- c(which(Np3>1e-4), which(Na3>1e-4))


highlight_color <- net %v% "color"  #  with base group colors
highlight_color[highlight_index] <- "firebrick"  # override with firebrick for highlights

# A label vector
highlight_label <- rep("", network.size(net))
highlight_label[highlight_index] <- network.vertex.names(net)[highlight_index]

net %v% "highlight.color" <- highlight_color
net %v% "label" <- highlight_label

# Plotting web
w_evo<-ggnet2(net, 
                mode = "circle", 
                size = 5, 
                color = "highlight.color", 
                label.color = "black",
                label.size = 5,
                edge.size = 2,
                edge.color = "gray40",
                edge.alpha = 0.8)


########################### Temperature shift - evolution trait + evo variance ######################3


bw  <-2
aw  <- 0
gi <- 1.5
ki <-0.4#mortality rate
w<- 1
ic <-c(Na1, Np1, ma1, mp1, sa , sp )
mut.strength<-1
species_index <- which(degree==max(degree))

var_evolve<-1
trait_evolve<-1
params <- list(time=time,matrix=g,Temp=26,tmax=tmax,species_index=species_index,
              sig=sig,bw=bw,gi=gi,ki=ki,variation="high",
               web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=var_evolve,
               mut.strength=mut.strength,C=C,nestedness=nestedness,h2="evo_trait_and_variance",
               web.name=web.name,dganimals=dganimals,degree=degree,w_c=0.2,
               dgplants=dgplants,envA=envA,envP=envP)

species_index <- which(degree==max(degree))

tmax<-2000
OUT_evotrait_var_net1<-ode(func=eqs_type2_new, y=ic, parms=params,
                           times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params) 


plot_density(OUT_evotrait_var_net1)

OUT$time
OUT_evotrait_var_net1$time<-OUT_evotrait_var_net1$time+tmax

new_dat_evo_trait_var<-rbind(OUT,OUT_evotrait_var_net1)


#just plotting timeseries
r4<-plot_density_timeseries(new_dat_evo_trait_var)+
  geom_vline(xintercept = 2000, linetype="dashed",size=1, color="grey5")+
  scale_color_viridis_d(alpha = 1)

r4

#subsetting for final timepoint
Na4<- (OUT_evotrait_var_net1 %>% filter(time == 4e3,type=="N"))$v
Np4<-(OUT_evotrait_var_net1 %>% filter(time == 4e3,type =="P"))$v
ma4<-(OUT_evotrait_var_net1 %>% filter(time == 4e3,type == "muA"))$v
mp4<-(OUT_evotrait_var_net1 %>% filter(time == 4e3,type == "muP"))$v
sa4<-(OUT_evotrait_var_net1 %>% filter(time == 4e3,type == "sA"))$v
sp4<-(OUT_evotrait_var_net1 %>% filter(time == 4e3,type == "sP"))$v

#phenotypic variance
Pa4<-sa4+envA
Pp4<-sp4+envP

#snapshot
t_evo_trait_var<- plot_snapshot(Na=Na4, Np=Np4, m = c(ma4,mp4), 
                                sigma = c(Pa4,Pp4), Temp=26, limits = c(19,29),res = 1001)+
  annotate("text",x=26, y =80, label="time == 4000",parse=T)  

t_evo_trait_var 


#plotting web for this model of evolution of trait and genetic variance
g<-adj.mat(webfiles[which(webfiles == "datasets_1/M_PL_061_31.csv")]) #ork web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}


#vertex names 
net = network(g, bipartite = T, directed = FALSE)
deg <- c(degree.plants, degree.animals)

# Assign groups and base color
net %v% "groups" <- ifelse(1:network.size(net) <= dim(g)[1], "plants", "animals")
net %v% "color" <- ifelse(net %v% "groups" == "plants", "#E69F00", "#E69F00")

# Choosing  nodes by index which has densities greater than 1e-4.
highlight_index <- c(which(Np4>1e-4), which(Na4>1e-4))


highlight_color <- net %v% "color" 
highlight_color[highlight_index] <- "firebrick"  # override with red for highlights


highlight_label <- rep("", network.size(net))
highlight_label[highlight_index] <- network.vertex.names(net)[highlight_index]


net %v% "highlight.color" <- highlight_color
net %v% "label" <- highlight_label

# ---- PLOT ----
w_evo_var_trait<-ggnet2(net, 
                mode = "circle", 
                size = 5, 
                color = "highlight.color", 
                label.color = "black",
                label.size = 5,
                edge.size = 2,
                edge.color = "gray40",
                edge.alpha = 0.8)




############################### temperature shift-- no evo ########################

bw  <-2
aw  <- 0
gi <- 1.5
ki <-0.4#mortality rate
w<- 1
ic <-c(Na1, Np1, ma1, mp1, sa , sp )
mut.strength<-1
species_index <- which(degree==max(degree))

var_evolve<-0
trait_evolve<-0
params <- list(time=time,matrix=g,Temp=26,tmax=tmax,species_index=species_index,sig=sig,bw=bw,gi=gi,ki=ki,variation="high",
               web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=var_evolve,
               mut.strength=mut.strength,C=C,nestedness=nestedness,h2="no_evo",
               web.name=web.name,dganimals=dganimals,degree=degree,w_c=0.2,
               dgplants=dgplants,envA=envA,envP=envP)

species_index <- which(degree==max(degree))

tmax<-2000
OUT_noevo_net1<-ode(func=eqs_type2_new, y=ic, parms=params,
         times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params) 


OUT_noevo_net1 %>% plot_density()

OUT$time
OUT_noevo_net1$time<-OUT_noevo_net1$time+2000

new_dat_noevo_trait<-rbind(OUT,OUT_noevo_net1)

r2<-plot_density_timeseries(new_dat_noevo_trait)+
  geom_vline(xintercept = 2000, linetype="dashed",size=1, color="grey5")+
  scale_color_viridis_d(alpha = 1)
r2


Na2<- (OUT_noevo_net1 %>% filter(time == 4e3,type=="N"))$v
Np2<-(OUT_noevo_net1 %>% filter(time == 4e3,type =="P"))$v
ma2<-(OUT_noevo_net1 %>% filter(time == 4e3,type == "muA"))$v
mp2<-(OUT_noevo_net1 %>% filter(time == 4e3,type == "muP"))$v
sa2<-(OUT_noevo_net1 %>% filter(time == 4e3,type == "sA"))$v
sp2<-(OUT_noevo_net1 %>% filter(time == 4e3,type == "sP"))$v

Pa2<-sa2+envA
Pp2<-sp2+envP
t_noevo<- plot_snapshot(Na=Na2, Np=Np2, m = c(ma2,mp2), 
                   sigma = c(Pa2,Pp2), Temp=26, limits = c(19,29),res = 1001)+
  annotate("text",x=22, y =5, label="time == 4000",parse=T)  

t_noevo 

Na2[which(Na2<1e-4)]<-0
Np2[which(Np2<1e-4)]<-0


g<-adj.mat(webfiles[which(webfiles == "datasets_1/M_PL_061_31.csv")]) #ork web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}


#vertex names 
net = network(g, bipartite = T, directed = FALSE)
deg <- c(degree.plants, degree.animals)

# Assign groups and base color
net %v% "groups" <- ifelse(1:network.size(net) <= dim(g)[1], "plants", "animals")
net %v% "color" <- ifelse(net %v% "groups" == "plants", "#E69F00", "#E69F00")

# ---- HIGHLIGHT TWO NODES BY INDEX ----
# Choose two nodes by index (e.g., 5 and 10)
highlight_index <- c(which(Np2>0), which(Na2>0))

# Create a custom color vector
highlight_color <- net %v% "color"  # start with base group colors
highlight_color[highlight_index] <- "firebrick"  # override with red for highlights

# Create a label vector
highlight_label <- rep("", network.size(net))
highlight_label[highlight_index] <- network.vertex.names(net)[highlight_index]

# Assign highlight attributes
net %v% "highlight.color" <- highlight_color
net %v% "label" <- highlight_label

# ---- PLOT of the network highligted by nodes which are either dead or alive 
#based on the color either firebrick or orange
w_noevo<-ggnet2(net, 
       mode = "circle", 
       size = 5, 
       color = "highlight.color", 
       label.color = "black",
       label.size = 5,
       edge.size = 2,
       edge.color = "gray40",
       edge.alpha = 0.8)



####### plottting all models and networks in one #########
ggpubr::ggarrange(w_noevo, t2, t_noevo, r2,ncol=4,
                 w_evo, t2, t_evo_trait, r3,
                 w_evo_var_trait,t2, t_evo_trait_var,r4, nrow=3,labels=c("A","B","C",
                                                             "D","E","F",
                                                             "G","H","I","J","K","L"))
