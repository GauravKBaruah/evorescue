#R script that has the functions that are needed for simulating population and evolutionary dynamics
#of plant-pollinatore networks
#Also has functions that are required to organise data and plotting

require(statmod)
cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)

#function that allows for transferring a .csv file to a adjacency matrix
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}

# function that takes in an adjacency matrix and gives the minimum path it takes to 
#get to a hyper generalist species in a network
min_no_path_to_target_node<- function(graph_matrix, start_node, target_node ){
  
  
  g1 <- graph_from_biadjacency_matrix(graph_matrix, weighted = T)
  #all_paths <- all_simple_paths(g1, from = start_node, to = target_node, mode = "out")
  shortest_path<- length((shortest_paths(g1, from = start_node, to = target_node))$vpath[[1]])-1

  return(shortest_path)
}

# function for numerically solve the mutualistic interaction type-2 functional curves
#na, #np : vector of animal and plant densities
#a_index: animal species i
# w : gaussian width of interaction
#degree.animal: degree of the animal
# Pa : phenotypic variance
#mat : adjacency matrix
#mut_strength : 1, gamma_0 in the paper
#sigma : phenotypic variance llist
#h : handling time
#this function uses Gaussian quadrature to numerically solve the type 2 mutualistic functional curve
#see https://en.wikipedia.org/wiki/Gaussian_quadrature for this.
type_2_animals<-function(m,sigma,w,h,np,na,mut.strength,a_index,
                         points,mat,degree.animal,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma[a_index], sigma =sqrt(sigma$sa[a_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma[a_index],sigma =sqrt(sigma$sa[a_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(np),ncol=points)
  w2<-matrix(0,nrow=length(np),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
  N_strength<-m_strength<-g_strength<-numeric()

  for(j in 1:points){  
    for(k in 1:length(np)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
      
      numer_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength)*sum((z1[j]-m$ma[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- np[k]*mat[k,a_index]*1/Pa[a_index]^2*(mut.strength)*sum(((z1[j]-m$ma[a_index])^2-Pa[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- np[k]*mat[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }

  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}
# function for numerically solve the mutualistic interaction type-2 functional curves
#na, #np : vector of animal and plant densities
#a_index: animal species i
# w : gaussian width of interaction
#degree.animal: degree of the animal
# Pa : phenotypic variance
#mat : adjacency matrix
#mut_strength : 1, gamma_0 in the paper
#sigma : phenotypic variance llist
#h : handling time
#this function uses Gaussian quadrature to numerically solve the type 2 mutualistic functional curve
#see https://en.wikipedia.org/wiki/Gaussian_quadrature for this.
type_2_plants<-function(m,sigma,w,h,np,na,mut.strength,p_index,
                        points,mat,degree.plant,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp[p_index],sigma =sqrt(sigma$sp[p_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(na),ncol=points)
  w2<-matrix(0,nrow=length(na),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
  N_strength<-m_strength<-g_strength<-numeric()

  for(j in 1:points){  
    for(k in 1:length(na)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
      
      numer_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength)*sum((z1[j]-m$mp[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- na[k]*mat[p_index,k]*1/Pa[p_index]^2*(mut.strength)*sum(((z1[j]-m$mp[p_index])^2-Pa[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- na[k]*mat[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}


#function for the ODE solver: that simulates the eco-evo variance dynamics
#of plant-pollinator networks
#time: tmax, maximum timepoint for simulation
#state: initial state values i.e., Na, Np, ma,mp,va,vp, all state variables for densities, 
#mean trait values, variance of species.
eqs_type2_new<- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1] #no of plants
  Np<-state[(A+1):(A+P)] #initial state variable for plants
  Na<-state[1:A] #initial density for animals
  s <- pars$sig #phenotypic variances
  ma<-state[(A+P+1):(A+P+A)] #muA : mean trait values for animals
  mp<-state[(A+P+A+1):(A+P+A+P)] # mean trait values for plants
  gvarA<-state[(A+P+A+P+1):(A+P+A+P+A)] #genetic variance for animals
  gvarP<-state[(A+P+A+P+A+1):(A+P+A+P+A+P)]#genetic variance for plants
  
  aij<-bij<-vij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-vji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-vj<-vi<-numeric()
  #*cutoff(gvarA/(1e-7))
  varA <- ( (gvarA + pars$envA))
  varP <- ( (gvarP + pars$envP))
  
  Temp<- pars$Temp# local environmental temperature

  dmA <- outer(muA, muA, FUN="-") ## difference matrix of trait means
  dmP <- outer(muP, muP, FUN="-") ## difference matrix of trait means
  
  svA <- outer(varA, varA, FUN="+") ## sum matrix of trait variances
  svP <- outer(varP, varP, FUN="+") ## sum matrix of trait variances
  
  alphaA <- exp(-dmA^2/(2*svA+pars$w_c^2))*pars$w_c/sqrt(2*svA+pars$w_c^2) ## alpha matrix
  alphaP <-  exp(-dmP^2/(2*svP+pars$w_c^2))*pars$w_c/sqrt(2*svP+pars$w_c^2) ## alpha matrix

  diag(alphaA)<-1 #intra specific competition animals
  diag(alphaP)<-1 #intra specific competition plants
  
  betaA <- alphaA*2*varA*(-dmA)/(2*svA+pars$w_c^2)^1.5 ## beta matrix
  betaP <- alphaP*2*varP*(-dmP)/(2*svP+pars$w_c^2)^1.5 ## beta matrix
  diag(betaA)<-0
  diag(betaP)<-0

  xetaA  <- alphaA/(2*svA+pars$w_c^2)^2*(dmA^2-(2*svA+pars$w_c^2)) 
  xetaP <- alphaP/(2*svP+pars$w_c^2)^2*(dmP^2-(2*svP+pars$w_c^2))
  diag(xetaA)<-0 
  diag(xetaP)<-0
  

  # temperature tolerance function on growth rate for animals on density, and on mean trait, and variance
  ba<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+varA))*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA)) - pars$ki
  bar_ba<- pars$gi/(pars$bw)*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA))*(varA)*(pars$bw)*(Temp -muA)/((pars$bw)^2 + varA)^1.5
  xar_ba<- (2*pars$gi*(2*muA^2-4*Temp*muA+2*Temp^2-2*pars$bw^2-varA)*exp(-(Temp-muA)^2/(2*pars$bw^2+varA)))/(sqrt(pars$bw^2+varA)*(2*pars$bw^2+varA)^2)
  
  bp<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+ varP))*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP)) - pars$ki
  bar_bp<-pars$gi/(pars$bw)*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP))*(varP)*(pars$bw)*(Temp - muP)/((pars$bw)^2 + varP)^1.5 
  xar_bp<- (2*pars$gi*(2*muP^2-4*Temp*muP+2*Temp^2-2*pars$bw^2-varP)*exp(-(Temp-muP)^2/(2*pars$bw^2+varP)))/(sqrt(pars$bw^2+varP)*(2*pars$bw^2+varP)^2)
  
  #for loop for calculating numerically the double integral of type 2 functional curve of mutualistic interaction
  #we don not ahave an analytical result, maybe there is better way to do this, but i didnt delve deep to find out.
  for(r in 1:A){
  
      m.temp<-list(ma=muA,mp=muP)
      sigma1<-list(sa=varA,sp=varP)
      #type 2 function curve (see above the function description)
      temp1<-type_2_animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np,na=Na,mut.strength=pars$mut.strength,a_index=r,
                     points=6,mat=pars$matrix,degree.animal=pars$dganimals[r],Pa=varA)
    ai[r]<-temp1$G
    bi[r]<-temp1$B
    vi[r]<-temp1$V
  }
  #for loop for calculating numerically the double integral of type 2 functional curve of mutualistic interaction
  #we don not ahave an analytical result, maybe there is better way to do this, but i didnt delve deep to find out.
  for(k in 1:P){
   
      m2.temp<-list(ma=muA,mp=muP)
      sigma2<-list(sa=varA,sp=varP)
      #type 2 function curve (see above the function description)
      temp2<-type_2_plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np,na=Na,
                             mut.strength=pars$mut.strength,p_index=k,
                             points=6,mat=pars$matrix, 
                             degree.plant =pars$dgplants[k],
                             Pa=varP)

    aj[k]<-temp2$G
    bj[k]<-temp2$B
    vj[k]<-temp2$V
  }
# state dynamics : dn(A)/dt, dn(P)/dt, dm(A)/dt, dm(P)/dt, dv(P)/dt, dv(A)/dt
  dndt_a<- Na*(ba-alphaA%*%Na+ai)*cutoff(Na/(1e-6))  #population dynamics
  dndt_p<- Np*(bp-alphaP%*%Np+aj)*cutoff(Np/(1e-6))  #population dynamics
  dudt_A<- pars$trait_evolve*gvarA/varA*(bar_ba- betaA%*%Na+ bi) #mean trait dynamics
  dudt_P<- pars$trait_evolve*gvarP/varP*(bar_bp -betaP%*%Np+ bj) #mean trait dynamics
  dudt_vA<- pars$var_evolve*1/2*(gvarA/varA)^2*(xar_ba - xetaA%*%Na + vi)*cutoff(gvarA/(1e-7)) # genetic variance dynamics
  dudt_vP<- pars$var_evolve*1/2*(gvarP/varP)^2*(xar_bp - xetaP%*%Np + vj)*cutoff(gvarP/(1e-7)) # genetic variance dynamics
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P,dudt_vA,dudt_vP)))
}


create_matrix<-function(SA,SP){
  web <- matrix(rbinom(SA*SP, 1, prob=1),nrow=SA,ncol =SP)
  return(web)
  
}


#this is a function that was used for the Theobiota Cluster that takes the ODE solver function,
# and final state variables from the ODE solver functions.
#also gives the categories of Rescue, Decline, Decline-stasis, Stasis, Growth.
#this is a messy function, but one that only uses output data from the cluster analysis to caategories species in different
#categories, calculates nestedness, connectance,trait-lag, overall competition strenght etc.
cluster_run_func<-function(params, ic, tmax ){
  
  OUT<-ode(func=eqs_type2_new, y=ic, parms=params,
           times=seq(0, tmax, by=10)) %>% 
    organize_results(pars = params) 
  (Na_r<-(OUT  %>% filter(time == tmax,type=="N"))$v)
  (Np_r<-(OUT %>% filter(time == tmax,type=="P"))$v)
  (mua_r<-(OUT  %>% filter(time == tmax,type=="muA"))$v)
  (mup_r<-(OUT %>% filter(time == tmax,type=="muP"))$v)
  (va_r<-(OUT  %>% filter(time == tmax,type=="sA"))$v)
  (vp_r<-(OUT %>% filter(time == tmax,type=="sP"))$v)
  
   Na_r[(Na_r < 1e-7)] <- 0
   Np_r[(Np_r < 1e-7)] <- 0
   
   A<-length(Na_r)
   P<-length(Np_r)
   na<-ic[1:length(Na_r)]
   np<-ic[(length(Na_r)+1) : (length(Na_r)+length(Np_r) )]

   N_initial<-c(na,np)   
  N_dat_min<-(OUT  %>% filter(  time >= 10 & time <= tmax , type=="N") %>% select(species, v) %>%  
    group_by(species) %>% summarise(min_n = min(v)))$min_n
  P_dat_min<-(OUT  %>% filter(  time >= 10 & time <= tmax , type=="P") %>% select(species, v) %>%  
    group_by(species) %>% summarise(min_p = min(v)))$min_p
 
  max_Na_t <-(OUT  %>% filter( time == tmax, type=="N"))$v 
  max_Np_t <-(OUT  %>% filter( time == tmax, type=="P"))$v 
  
  N_min <-c(N_dat_min,P_dat_min)
  N_final <-c(max_Na_t,max_Np_t)
 
  N_dat_min[(N_dat_min< 1e-6)]<-0
  P_dat_min[(P_dat_min< 1e-6)]<-0
  
  max_Na_t[(max_Na_t< 1e-6)]<-0
  max_Np_t[(max_Np_t< 1e-6)]<-0
 
  #Rescue metric
  index_N<- which((N_min+0.1) < N_final & (N_min+0.1) < N_initial & N_initial>1e-4 & N_final > 1e-4)
  fraction_evo_rescue <- length(index_N)/(A+P)

  N_evorescue<-N_initial
  N_evorescue[index_N] <- 1
  N_evorescue[which(N_evorescue != 1)] <- 0
  
  #this following variables utilise the descriptive figure S11 to categories species into Decline, Stasis, Growth, Decline-stasis.
  index_N_declines<- which(N_initial>N_min & N_min >= N_final & N_final < 1e-4 & N_initial>1e-4)  #Decline metric
  index_stasis <- which(N_initial>1e-4 & N_initial < (N_min+0.1) & (N_min+0.1) > N_final & N_final > 1e-4 ) #stasis metric
  index_growth<-which(N_initial>1e-4 & N_initial < N_final & (N_min + 0.1) > N_initial & N_final>(N_min+0.1) & N_final > 1e-4) #growth metric 
  index_decline_stasis<- which( N_initial > 1e-4 &  N_initial>N_final & N_final < (N_min + 0.1) & N_initial> (N_min+0.1) & N_final > 1e-4 ) #decline-stasis metric
  
  
  ## fraction of all the species into the categories
  fraction_decline <- length(index_N_declines)/(A+P)
  fraction_stasis <- length(index_stasis)/(A+P)
  fraction_growth<- length(index_growth)/(A+P)
  fraction_decline_stasis<- length(index_decline_stasis)/(A+P)
  
  N_decline <- N_initial
  N_decline[index_N_declines]<-1
  N_decline[which(N_decline!=1)]<-0
 
  
  N_stasis<-N_initial
  N_stasis[index_stasis]<-1
  N_stasis[which(N_stasis!=1)]<-0


  N_growth<-N_initial
  N_growth[index_growth]<-1
  N_growth[which(N_growth!=1)]<-0

  N_decline_stasis<-N_initial
  N_decline_stasis[index_decline_stasis]<-1
  N_decline_stasis[which(N_decline_stasis!=1)]<-0

  final_density<-c(Na_r,Np_r)
  g<-params$matrix
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] 
  proportion_evo_rescue<-fraction_evo_rescue
  na<-ic[1:length(Na_r)]
  np<-ic[(length(Na_r)+1) : (length(Na_r)+length(Np_r) )]
  muA<-ic[(Aspecies+Plantspecies+1):(Aspecies+Plantspecies+Aspecies)]
  muP<-ic[(Aspecies+Plantspecies+Aspecies+1):(Aspecies+Plantspecies+Aspecies+Plantspecies)]
  Pa<-ic[(Aspecies+Plantspecies+Aspecies+Plantspecies+1):(Aspecies+Plantspecies+Aspecies+Plantspecies+Aspecies)]
  Pp<-ic[(Aspecies+Plantspecies+Aspecies+Plantspecies+Aspecies+1):(Aspecies+Plantspecies+Aspecies+Plantspecies+Aspecies+Plantspecies)]
  
  
  #quantifying competition strength:
  dmA <- outer(muA, muA, FUN="-") ## difference matrix of trait means
  dmP <- outer(muP, muP, FUN="-") ## difference matrix of trait means
  
  svA <- outer(Pa, Pa, FUN="+") ## sum matrix of trait variances
  svP <- outer(Pp, Pp, FUN="+") ## sum matrix of trait variances
  
  
  comp_str_A<- exp(-dmA^2/(2*svA+params$w_c^2))*params$w_c/sqrt(2*svA+params$w_c^2)
  comp_str_P<- exp(-dmP^2/(2*svP+params$w_c^2))*params$w_c/sqrt(2*svP+params$w_c^2)
  cumulative_comp_A<-rowMeans(comp_str_A)
  cumulative_comp_P<-rowMeans(comp_str_P)
  cumulative_comp<-c(cumulative_comp_A,cumulative_comp_P)
  
  #trait lag
  trait_lag_A <- params$Temp - muA  
  trait_lag_P <- params$Temp - muP
  trait_lag<-c(trait_lag_A,trait_lag_P)
  
  #Absolute trait value at equilibrium
  eq_muA<-muA
  eq_muP<-muP
  eq_mu<-c(muA,muP)

  network_metric <- networklevel(params$matrix)
  graph_mat<-graph_from_biadjacency_matrix(t(g),weighted = T)
  
  extinctions <- length(which( final_density < 1e-4))/(Aspecies+Plantspecies)
  change_in_density <- sum(final_density)-sum(c(na,np))
  richness <- length(which(final_density>1e-4) )/(Aspecies+Plantspecies)
  change_richness <- length(which(final_density>1e-4) ) -length(which(c(na,np)>1e-4) ) 
  community_biomass<- sum(final_density)
  change_community_biomass<- sum(final_density) - sum(c(na,np))
  NODF <- nestedness_NODF(g)
  connectance <- Connectance(g)
  network_size <- (Aspecies+Plantspecies)
  modularity<-network_metric[[7]]
  avg_betweeness <- mean(igraph::betweenness(graph_mat))
  sd_betweeness <- sd(igraph::betweenness(graph_mat))
  
  sa_sp<-OUT %>% filter(type %in%c("sA","sP"))
   dd<-sa_sp %>%
    group_by(type, species) %>%
    summarize(max_v = max(v, na.rm = TRUE))


  out  <-data.frame(webname=params$web.name,
                    h2=params$h2,
                    mean_trait_lag=mean(trait_lag),
                    fraction_decline=fraction_decline,
                    fraction_growth=fraction_growth,
                    fraction_evo_rescue=fraction_evo_rescue,
                    fraction_stasis=fraction_stasis,
                    fraction_decline_stasis=fraction_decline_stasis,
             change_in_density=change_in_density,
             proportion_evo_rescue=proportion_evo_rescue,
             temperature_shift= params$Temp,
             w_c=params$w_c,
             ki=params$ki,
             richness=richness,
             change_richness=change_richness,
             community_biomass=community_biomass,
             change_community_biomass=change_community_biomass,
             Nestedness=NODF,
             connectance=connectance,
             network_size=network_size,
             avg_betweeness=avg_betweeness,
             sd_betweeness=sd_betweeness,
             modularity=modularity)
  
  species_index<-params$species_index
  final_density_1<-final_density
  final_density_1[which(final_density_1 < 1e-4)]<-0
  final_density_1[which(final_density_1 > 1e-4)]<-1
  extinct_binomial<-final_density_1
  proportion_survived<-sum(final_density_1)/(Aspecies+Plantspecies)
  rownames(g) <- as.character(seq(1:nrow(g)))
  Ls<-c(letters,LETTERS, seq(500,800,1))
  colnames(g) <- Ls[c(1:ncol(g))]
  index_names_of_g<-c(Ls[c(1:ncol(g))],as.character(seq(1:nrow(g))))   
  target_generalist<-index_names_of_g[species_index]
  min_no_path_to_generalist<-numeric()
  for(i in 1:(Aspecies+Plantspecies)){
    
    min_no_path_to_generalist[i]<- min_no_path_to_target_node(graph_matrix = g,start_node = index_names_of_g[i],
                                                              target_node = target_generalist)
    
    
  }
  # making sure the index of the generalist path to itself is zero
  min_no_path_to_generalist[species_index] <- 0
  #min_no_path_to_generalist[is.na(min_no_path_to_generalist)]<- Aspecies+Plantspecies
  
  percentage_change_biomass<- ((c(Na_r,Np_r) - c(na,np))/c(na,np))*100
  change_biomass<-c(Na_r,Np_r)-c(na,np)
 
 species_level_data <-data.frame(
    rep(seq(1,(nrow(g)+ncol(g)),1)), #number of species
    c(muA,muP),
    (c(mua_r, mup_r) - c(muA,muP)),
    (c(va_r,vp_r) - c(Pa,Pp)),
    c(params$envA,params$envP),
    as.numeric(c(Na_r,Np_r )),
    as.numeric(proportion_survived),
    as.numeric(fraction_evo_rescue),
    as.numeric(fraction_decline),
    as.numeric(fraction_growth),
    as.numeric(fraction_stasis),
    as.numeric(fraction_decline_stasis),
    dd$max_v,
    cumulative_comp,
    trait_lag,
    N_decline,
    N_stasis,
    N_growth,
    N_evorescue,
    N_decline_stasis,
    params$mut.strength[1],
    params$ki,
    params$Temp,
    params$w,
    params$w_c,
    percentage_change_biomass,
    change_biomass,
    igraph::betweenness(graph_mat),
    as.numeric(min_no_path_to_generalist),
    rep(as.character(params$web.name),each=(Aspecies+Plantspecies)),
    params$degree,
    as.character(params$variation),
    as.numeric(rep(nestedness_NODF(params$matrix), each=((Aspecies+Plantspecies)) )),
    as.numeric(rep(Connectance(params$matrix), each=((Aspecies+Plantspecies)) )),
    as.numeric(rep( (Aspecies+Plantspecies),each=((Aspecies+Plantspecies)))),
    as.character(params$h2),
    rep(as.character(params$forcing),each=(Aspecies+Plantspecies)),
    as.numeric(modularity))
  
  
  
  colnames(species_level_data)<-c("Species",
				 "mean_trait",
                                  "change_in_mean_trait",
                                  "change_in_mean_gvariance",
                                  "environmental_variance",
                                  "density", 
                                  "proportion_persisted",
                                  "proportion_rescued",
                                  "proportion_declined",
                                  "proportion_growth",
                                  "proportion_stasis",
                                  "proportion_decline_stasis",
				 "max_genetic_variance",
                                  "competition_strength",
                                  "trait_lag",
                                  "decline",
                                  "stasis",
                                  "growth",
                                  "evorescue",
			          "decline_stasis",
                                  "mutualism_strength", 
                                  "mortality",
                                  "Temperature_shift",
                                  "width_mutualism",
                                  "width_competition",
                                  "percent_change_biomass",
				  "change_biomass",
                                  "betweeness_centrality",
                                  "minimum_path_generalist",
                                  "webname",
                                  "Degree",
                                  "variation",
                                  "Nested`ness", 
                                  "Connectance",
                                  "Network_size",
                                  "h2",
                                  "forcing",
                                  "modularity")

  final_output<- list(output=out, species_level_data=species_level_data)
  return(final_output)
  
}



## Organize simulation results into tidy table
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma
organize_results <- function(sol, pars) {
  S <- length(pars$sig) ## number of species
  A<-dim(pars$matrix)[2] # no. of animals
  P<-dim(pars$matrix)[1] # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  temp<- temp 
  names(temp)[1] <- "time"
  names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)

  names(temp)[(A+2):(A+1+P)] <- paste0("P_", 1:P) ## name trait mean columns
  names(temp)[(A+P+2):(A+P+A+1)] <-  paste0("muA_", 1:(A))
  names(temp)[(A+P+A+2):(A+P+A+P+1)] <- paste0("muP_", 1:(P))
  names(temp)[(A+P+A+P+2):(A+P+A+P+A+1)] <- paste0("sA_", 1:(A))
  names(temp)[(A+P+A+P+A+2):(A+P+A+P+A+P+1)]<-paste0("sP_", 1:(P))
   
  temp <- temp %>%
   tidyr::gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
    tidyr::separate(variable, c("type", "species"), sep="_") %>%
    #spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species,v) %>% ## rearrange columns
    mutate(species=as.integer(species), w=pars$gamma,
           Nestedness=pars$nestedness, Connectance=pars$C,
           theta=pars$theta,Web.name=pars$web.name) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}



## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
## used to produce figure 1.
plot_density<- function(dat) {
  dat %>%
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
    scale_y_continuous(name="population density") +
    theme(legend.position="none") + facet_wrap(.~type,scales = "free") %>%
    return
}

## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
plot_density_timeseries<- function(dat) {
  dat %>% filter(type == c("N","P")) %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
    scale_y_continuous(name="population density") +
    theme_cowplot()+
    theme(legend.position="none") + 
    facet_wrap(.~type,scales = "free") %>%
    return
}

## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot for genetic variable dynamics
plot_density_timeseries_genvar<- function(dat) {
  dat %>% filter(type == c("sA","sP")) %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
    scale_y_continuous(name="Genetic variance") +
    theme_cowplot()+
    theme(legend.position="none") + 
    facet_wrap(.~type,scales = "free") %>%
    return
}

# plots the density distribution  at a particular timepoint. This function was used to produce figure 1.
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# moment: mean
# limits: limits of the mean trait axis which in the study are -1,1
# code adapted from Barabas & D'Andrea 2016 Eco. Letts.
plot_snapshot <- function(Na, Np, m, sigma, Temp, limits=c(-1, 1), res=1001) {
  S_a <-length(Na) ## number of species
  S_p <- length(Np)
  ma<- m[1:(S_a)]
  mp<- m[(S_a+1):(S_a+S_p)]
  sigma_a <-sigma[1:(S_a)]
  sigma_p <- sigma[(S_a+1):(S_a+S_p)]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:S_a, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=(S_a+1):(S_a+S_p), trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:S_a) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (j in 1:S_p) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==((S_a)+j))] <- Np[j]*dnorm(traits_p$trait[(traits_p$species==((S_a)+j))], 
                                                                     mp[j], sigma_p[j]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA

  landscape <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= (1.5/2*exp(-(Temp-trait)^2/(2*2^2)))-0.2) %>%      
    mutate(r=ifelse(r<=0, NA, r)) %>%
    mutate(r=r*max(c(traits_a$density,traits_p$density),na.rm = T))
  
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d(alpha = 1)+
    facet_wrap(.~species_group, nrow = 2)+
    theme_classic()+
    theme(legend.title = element_text(size = 14), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14))+
    geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
              colour="black", alpha=1, na.rm=TRUE) +
    scale_x_continuous(name=expression("Trait value (" ~ degree*C ~ ")"), limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    scale_color_viridis_d(alpha=1)+
    theme(legend.position="none") %>%
    return 
}



#computes the raw NODF taken from Song et al 2017 J. Animal Ecology
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}



# measures connectance of a web network
Connectance<-function(web)
{
  return(sum(web)/(ncol(web)*nrow(web)))}



#function to estimate per capita growth of plants
#z: trait values
#Np, Na: densities
#pars: parameter values - list
#Temp: temperature
#muA,muP: mean trait values of animals and plants
#varA, varP: variance of animals and variance of plants 
#it returns: for z values, and per-capita growth rate i.e, r(z)^{P}

fitness_f_P<-function(z,Np,Na, pars,Temp, muA,muP, varA, varP){
  
  b <- pars$gi/pars$bw*exp(-(Temp-z)^2/(2*pars$bw)) 
  
  m2.temp<-list(ma=muA,mp=muP)
  sigma2<-list(sa=varA,sp=varP)
  mut_temp<-azz<-numeric()
  for(i in 1:length(z)){
    azz[i] <-   sum(pars$w_c/sqrt(pars$w_c^2+ 2*varP)*exp(-(z[i]-muP)^2/(pars$w_c^2+2*varP))*Np) 
  }
  mut_temp<- type_2_plants_s1fig(m=m2.temp,sigma=sigma2,w=pars$w,h=0.1,np=Np,na=Na,
                                 mut.strength=pars$mut.strength,
                                 points=length(z),mat=pars$matrix,degree.plant=mean(pars$dgplants),z = z)$G
  f_g <- b - mean(Np) - azz + mut_temp
  
  return(list(z=z,F_g=f_g))  
}


#function to estimate per capita growth of animals
#z: trait values
#Np, Na: densities
#pars: parameter values - list
#Temp: temperature
#muA,muP: mean trait values of animals and plants
#varA, varP: variance of animals and variance of plants 
#it returns: for z values, and per-capita growth rate i.e, r(z)
fitness_f_A<-function(z,Np,Na, pars,Temp,muA,muP, varA, varP){
  
  b <- pars$gi/pars$bw*exp(-(Temp-z)^2/(2*pars$bw)) 
  
  m2.temp<-list(ma=muA,mp=muP)
  sigma2<-list(sa=varA,sp=varP)
  mut_temp<-azz<-numeric()
  for(i in 1:length(z)){
    azz[i] <-   sum(pars$w_c/sqrt(pars$w_c^2+ 2*varA)*exp(-(z[i]-muA)^2/(pars$w_c^2+2*varA))*Na) 
  }
  mut_temp<- type_2_animals_s1fig(m=m2.temp,sigma=sigma2,w=pars$w,h=0.1,np=Np,na=Na,
                                  mut.strength=pars$mut.strength,
                                  points=length(z),mat=pars$matrix,degree.animal=mean(pars$dganimals),z = z)$G
  
  
  
  
  f_g <- b - mean(Na) - azz + mut_temp
  
  return(list(z=z,F_g=f_g))  
}



# plots the fitness landscape at a particular timepoint. 
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# Temp: local environmental temperature
# limits: limits of the mean trait axis which in the study
plot_snapshot_per_cap_growth <- function(Na, Np, m, sigma, Temp, pars,limits=c(-1, 1), res=301) {
  S_a <- length(Na) ## number of species
  S_p <- length(Np)
  ma<- m[1:(S_a)]
  mp<- m[(S_a+1):(S_a+S_p)]
  sigma_a <-sigma[1:(S_a)]
  sigma_p <- sigma[(S_a+1):(S_a+S_p)]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:S_a, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=(S_a+1):(S_a+S_p), trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
 Na<-Na/(sum(Na))
  Np<-Np/sum(Np)
  for (i in 1:S_a) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:S_p) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==((S_a+1)+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==((S_a+1)+i))], 
                                                                     mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  f_P<-fitness_f_P(z = traitaxis,Np = Np,Na = Na,pars = pars,Temp=Temp,muA = ma,muP = mp,
                   varA =sigma_a,varP=sigma_p)$F_g
  f_A<-fitness_f_A(z = traitaxis,Np = Np,Na = Na,pars = pars,Temp=Temp,muA = ma,muP = mp,
                   varA =sigma_a,varP = sigma_p)$F_g
  
  landscape_a <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= f_A)
  
  landscape_p <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= f_P)
  
  landscape<-data.frame(rbind(landscape_a,landscape_p), 
                        species_group=c(rep("Animals", nrow(landscape_a)),
                                        rep("Plants", nrow(landscape_p))))
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(landscape) + ## generate plot
    #geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    # geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
    #            alpha=0.15, colour=NA)
    geom_line(data=landscape, aes(x=trait, y=r),
              colour="firebrick", alpha=1, na.rm=TRUE) +
    scale_fill_viridis_d(alpha = 1)+
    facet_wrap(.~species_group, nrow = 2,scales = "free")+
    theme_classic()+
    theme(legend.title = element_text(size = 14), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14))+
    geom_hline(yintercept = 0,linetype="dashed", color="black")+
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous( name="Per-capita growth") +
    scale_color_viridis_d(alpha=1)+
    theme(legend.position="none") %>%
    return 
}



# type 2 numerical approximation function that uses Gaussian quadrature that is used to produce figure Figure 2
type_2_animals_s1fig<-function(m,sigma,w,h,np,na,mut.strength,
                               points,mat,degree.animal,z){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  z1<- z 
  z2<-matrix(0,nrow=length(np),ncol=points)
  w2<-matrix(0,nrow=length(np),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
  N_strength<-m_strength<-g_strength<-numeric()
  
  
  for(j in 1:points){  
    for(k in 1:length(np)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
      
      numer_a[j,k]<- np[k]*sum(mat[k,])*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- np[k]*sum(mat[k,])*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
       }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))
    }
  return(list(G= N_strength))
  
}


# type 2 numerical approximation function that uses Gaussian quadrature that is used to produce figure Figure 2
type_2_plants_s1fig<-function(m,sigma,w,h,np,na,mut.strength,
                              points,mat,degree.plant,z){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  
  z1<- z #gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
 
  z2<-matrix(0,nrow=length(na),ncol=points)
  w2<-matrix(0,nrow=length(na),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
  N_strength<-m_strength<-g_strength<-numeric()
  
  for(j in 1:points){  
    for(k in 1:length(na)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
      
      numer_a[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
  
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))
  }

  return(list(G=N_strength))
  
}



