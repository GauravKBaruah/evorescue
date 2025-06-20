
rm(list=ls())
# R script for Figure 6, and Figure S1, and Figure S5
#required functions and libraries
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

#load simulated data
load("final_dat_evorescue.RData")

#subsetting for competition width of 0.2
species_data1<-sp_data %>% filter(width_competition == 0.2)


#renaming some factors of models
species_data1$h2<-plyr::revalue(species_data1$h2, c("no_evo"= "No evolution",
                                                    "evo_trait" = "trait evolution",
                                                    "evo_trait_and_variance"="evolution trait and variance" ))

#structure of the data-stretched
str(species_data1)


#loading PCA libraries
library(corrplot)
library(factoextra)

#some variables in the data is not named properly such as `nestedness`
princpdata<-species_data1 %>% select(modularity,`Nested\`ness`,Connectance,Network_size)

#principal component analysis
res.pca <- prcomp(princpdata, scale = TRUE)
print(res.pca)
eig.val<-get_eigenvalue(res.pca)
eig.val
pc2<-fviz_eig(res.pca, col.var="blue")

# Color by cos2 values: quality on the factor map
PCA_PLOT<-fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#this is for figure S1 inn the appendix
ggpubr::ggarrange(PCA_PLOT,pc2)


# scores of the pca axis that are not uncorrelated now also known as loadings
PC1_scores<-res.pca$x[,1] 
PC2_scores<-res.pca$x[,2]

#checking the relationship between Network size and PCA score
plot(PC2_scores,species_data1$Network_size)


#creating two new columns with PC1 and PC2 axis scores - 
species_data1$PC1<-PC1_scores
species_data1$PC2<- PC2_scores



#now arranging the data for the final plot of the manuscript 
subset_data_decline <- subset(species_data1, select = c(decline, PC1, PC2, width_competition, 
                                                competition_strength, Connectance, h2, 
                                                Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean decline grouped by multiple variables
grouped_data_decline <- aggregate(decline ~ h2 + PC1 + PC2 + Network_size + 
                            width_competition + `Nested\`ness` + Connectance + 
                            Temperature_shift + webname, 
                          data = subset_data_decline, 
                          FUN = mean, na.rm = TRUE)


s1<-ggplot(grouped_data_decline, aes(x = PC1, y = decline, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
 ylab("Fraction of species showing decline  \n per network")+
  xlab("PC1")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  



subset_data_evorescue <- subset(species_data1, select = c(evorescue, PC1, PC2, width_competition, 
                                                competition_strength, Connectance, h2, 
                                                Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_evorescue <- aggregate(evorescue ~ h2 + PC1 + PC2 + Network_size + 
                            width_competition + `Nested\`ness` + Connectance + 
                            Temperature_shift + webname, 
                          data = subset_data_evorescue, 
                          FUN = mean, na.rm = TRUE)

s2<-ggplot(grouped_data_evorescue, aes(x = PC1, y = evorescue, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species showing rescue \n per network")+
  xlab("PC1")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  


## growth 

subset_data_growth <- subset(species_data1, select = c(growth, PC1, PC2, width_competition, 
                                                          competition_strength, Connectance, h2, 
                                                          Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_growth<- aggregate(growth ~ h2 + PC1 + PC2 + Network_size + 
                                      width_competition + `Nested\`ness` + Connectance + 
                                      Temperature_shift + webname, 
                                    data = subset_data_growth, 
                                    FUN = mean, na.rm = TRUE)

s3<-ggplot(grouped_data_growth, aes(x = PC1, y = growth, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species showing Growth \n per network")+
  xlab("PC1")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  



## decline and Growth

subset_data_stasis_decline<- subset(species_data1, select = c(decline_stasis, PC1, PC2, width_competition, 
                                                       competition_strength, Connectance, h2, 
                                                       Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_decline_stasis<- aggregate(decline_stasis ~ h2 + PC1 + PC2 + Network_size + 
                                  width_competition + `Nested\`ness` + Connectance + 
                                  Temperature_shift + webname, 
                                data = subset_data_stasis_decline, 
                                FUN = mean, na.rm = TRUE)

s4<-ggplot(grouped_data_decline_stasis, aes(x = PC1, y = decline_stasis, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing decline and stasis per network")+
  xlab("PC1")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  




## Stasis

subset_data_stasis<- subset(species_data1, select = c(stasis, PC1, PC2, width_competition, 
                                                              competition_strength, Connectance, h2, 
                                                              Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_stasis<- aggregate(stasis ~ h2 + PC1 + PC2 + Network_size + 
                                          width_competition + `Nested\`ness` + Connectance + 
                                          Temperature_shift + webname, 
                                        data = subset_data_stasis, 
                                        FUN = mean, na.rm = TRUE)

s5<-ggplot(grouped_data_stasis, aes(x = PC1, y = stasis, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing Stasis  per network")+
  xlab("PC1")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  

(s6<-species_data1 %>% 
    select(PC1,PC2,proportion_persisted,width_competition,competition_strength, Connectance, h2,Network_size, `Nested\`ness`, Temperature_shift) %>% 
    group_by(h2, PC1,PC2,Network_size,width_competition,`Nested\`ness`,Connectance,Temperature_shift) %>% 
    ggplot(aes(x=PC1, y=proportion_persisted,color=factor(h2)))+
    geom_point(size=3,alpha=0.8)+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.5 , alpha=0.1) +
    # geom_boxplot(alpha=0.5,aes(fill=fact1or(h2)),width = 0.5,size = 0.4)+
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    ylab("Fraction of species \n that survived per network")+
    xlab("PC1")+
    ylim(c(0,1))+
    labs(color="model")+
    theme_cowplot())


ggpubr::ggarrange(s1,s2,s3,s4, s5,s6, nrow=2,ncol = 3, labels = c("A","B","C","D","E","F"),
                  common.legend = T,legend = "bottom")



####################### FIgure S5 #####################
###### PC2- network size################################


subset_data_decline <- subset(species_data1, select = c(decline, PC1, PC2, width_competition, 
                                                        competition_strength, Connectance, h2, 
                                                        Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean decline grouped by multiple variables
grouped_data_decline <- aggregate(decline ~ h2 + PC1 + PC2 + Network_size + 
                                    width_competition + `Nested\`ness` + Connectance + 
                                    Temperature_shift + webname, 
                                  data = subset_data_decline, 
                                  FUN = mean, na.rm = TRUE)


s1<-ggplot(grouped_data_decline, aes(x = PC2, y = decline, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing decline per network")+
  xlab("PC2")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  



subset_data_evorescue <- subset(species_data1, select = c(evorescue, PC1, PC2, width_competition, 
                                                          competition_strength, Connectance, h2, 
                                                          Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_evorescue <- aggregate(evorescue ~ h2 + PC1 + PC2 + Network_size + 
                                      width_competition + `Nested\`ness` + Connectance + 
                                      Temperature_shift + webname, 
                                    data = subset_data_evorescue, 
                                    FUN = mean, na.rm = TRUE)

s2<-ggplot(grouped_data_evorescue, aes(x = PC2, y = evorescue, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing rescue per network")+
  xlab("PC2")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  


## growth 

subset_data_growth <- subset(species_data1, select = c(growth, PC1, PC2, width_competition, 
                                                       competition_strength, Connectance, h2, 
                                                       Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_growth<- aggregate(growth ~ h2 + PC1 + PC2 + Network_size + 
                                  width_competition + `Nested\`ness` + Connectance + 
                                  Temperature_shift + webname, 
                                data = subset_data_growth, 
                                FUN = mean, na.rm = TRUE)

s3<-ggplot(grouped_data_growth, aes(x = PC1, y = growth, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing Growth per network")+
  xlab("PC2")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  



## decline and Growth

subset_data_stasis_decline<- subset(species_data1, select = c(decline_stasis, PC1, PC2, width_competition, 
                                                              competition_strength, Connectance, h2, 
                                                              Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_decline_stasis<- aggregate(decline_stasis ~ h2 + PC1 + PC2 + Network_size + 
                                          width_competition + `Nested\`ness` + Connectance + 
                                          Temperature_shift + webname, 
                                        data = subset_data_stasis_decline, 
                                        FUN = mean, na.rm = TRUE)

s4<-ggplot(grouped_data_decline_stasis, aes(x = PC1, y = decline_stasis, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing decline and stasis per network")+
  xlab("PC2")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  




## Stasis

subset_data_stasis<- subset(species_data1, select = c(stasis, PC1, PC2, width_competition, 
                                                      competition_strength, Connectance, h2, 
                                                      Network_size, `Nested\`ness`, Temperature_shift, webname))

# Aggregate mean evorescue grouped by multiple variables
grouped_data_stasis<- aggregate(stasis ~ h2 + PC1 + PC2 + Network_size + 
                                  width_competition + `Nested\`ness` + Connectance + 
                                  Temperature_shift + webname, 
                                data = subset_data_stasis, 
                                FUN = mean, na.rm = TRUE)

s5<-ggplot(grouped_data_stasis, aes(x = PC1, y = stasis, color = h2)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.5 , alpha=0.1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab("Fraction of species \n showing Stasis   per network")+
  xlab("PC2")+
  labs(color="model")+
  ylim(c(0,1))+
  theme_cowplot()  

(s6<-species_data1 %>% 
    select(PC1,PC2,proportion_persisted,width_competition,competition_strength, Connectance, h2,Network_size, `Nested\`ness`, Temperature_shift) %>% 
    group_by(h2, PC1,PC2,Network_size,width_competition,`Nested\`ness`,Connectance,Temperature_shift) %>% 
    ggplot(aes(x=PC1, y=proportion_persisted,color=factor(h2)))+
    geom_point(size=3,alpha=0.8)+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.5 , alpha=0.1) +
    # geom_boxplot(alpha=0.5,aes(fill=fact1or(h2)),width = 0.5,size = 0.4)+
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    ylab("Fraction of species that survived \n per network")+
    xlab("PC2")+
    ylim(c(0,1))+
    labs(color="model")+
    theme_cowplot())


#figure S5
ggpubr::ggarrange(s1,s2,s3,s4, s5,s6, nrow=2,ncol = 3, labels = c("A","B","C","D","E","F"),
                  common.legend = T,legend = "bottom")













