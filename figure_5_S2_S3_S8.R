
rm(list=ls())

#required libraries
require(deSolve)
require(cowplot)
library(dplyr)
library(readr)
library(beepr)
library(ggplot2)
library(viridis)
library(igraph)
library(bipartite)
library(plyr)


#load simulated data
load("final_dat_evorescue.RData")



species_data1<- sp_data %>%  filter(width_competition == 0.2) #all_data_spdat2

#arranging the data for plotting
df_long <- species_data1 %>% select(Species,h2,density,width_competition,competition_strength,trait_lag,max_genetic_variance,
                                    change_in_mean_trait,change_in_mean_gvariance,
                                    Degree,decline,evorescue,stasis,decline_stasis,growth,minimum_path_generalist) %>% 
  tidyr::pivot_longer(cols=c(decline, evorescue,stasis,growth,decline_stasis),
               names_to="category", values_to = "Value") 

df_summary <- df_long %>%
  group_by(Degree,change_in_mean_trait,density,change_in_mean_gvariance,competition_strength,trait_lag,minimum_path_generalist, Value,
           width_competition, h2, category) %>%
  dplyr::summarise(Count = n()) %>%
  ungroup()

#renaming some columns to manuscript appropriate names
df_summary$h2<- revalue(df_summary$h2, c("no_evo"= "No evolution", "evo_trait" = "trait evolution", 
                                         "evo_trait_and_variance"="evolution trait and variance" ))


# remanimg some categories 
df_summary$category<- plyr::revalue(df_summary$category, c("decline"="decline",
                                                           "decline_stasis"="decline_stasis",
                                                           "evorescue"="rescue",
                                                           "growth"="growth",
                                                           "stasis"="stasis"))
# Competition strength vs. response categories
(a1<-df_summary %>% 
    ggplot( aes(x = competition_strength, y = Value, color = h2)) +
    geom_point(size=1.5, alpha=0.25) + 
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.2 , alpha=0.15) +
    xlim(c(0,0.39))+
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    # Size of dots corresponds to the count
    facet_wrap(.~category,ncol=5,nrow=1) +  # Facet grid by h2 and width_competition
    labs(x = "Overall competition strength faced by a species", title="B",
         y = "Chances of a metric", size = "Frequency", color = "model") +
    scale_y_continuous(breaks = c(0, 1)) +
    scale_x_continuous(breaks = c(0, 0.1,0.2, 0.3)) +# Ensure y-axis only shows 0 and 1
    theme_classic() +
    theme(legend.position = "bottom"))


# Trait lag vs response categories
(a2<- df_summary %>% 
    ggplot( aes(x = trait_lag, y = Value, color = h2)) +
    geom_point(size=1.5, alpha=0.25) + 
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.2 , alpha=0.15) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    # Size of dots corresponds to the count
    facet_wrap(.~category,ncol=5,nrow=1)+ # Facet grid by h2 and width_competition
    labs(x = "Species trait lag", y = "Chances of a metric", title="C", color = "model") +
    scale_y_continuous(breaks = c(0, 1)) +  # Ensure y-axis only shows 0 and 1
    theme_classic() +
    theme(legend.position = "bottom"))


# Degree vs response categories
(a3<-df_summary  %>% 
    ggplot( aes(x = Degree, y = Value, color = h2)) +
    geom_point(size=1.5, alpha=0.25) + 
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.2 , alpha=0.15) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    # Size of dots corresponds to the count
    facet_wrap(.~ category,ncol=5,nrow=1)+ # Facet grid by h2 and width_competition
    labs(x = "Species Degree", y = "Chance of a metric", title="D",size = "Frequency", color = "model") +
    scale_y_continuous(breaks = c(0, 1)) +  # Ensure y-axis only shows 0 and 1
    theme_classic() +
    theme(legend.position = "bottom"))


# minimum path generalist vs. response categories
(a4<-df_summary  %>%  
  ggplot( aes(x = minimum_path_generalist, y = Value, color = h2)) +
    geom_point(size=1.5, alpha=0.25) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.2 , alpha=0.15) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  # Size of dots corresponds to the count
  facet_wrap(. ~ category, ncol=5,nrow=1) +  # Facet grid by h2 and width_competition
  labs(x = "least path to super generalist", y = "Chance of a metric", 
       size = "Frequency", color = "model", title="E") +
  scale_y_continuous(breaks = c(0, 1)) +  # Ensure y-axis only shows 0 and 1
  theme_classic() +
  theme(legend.position = "bottom"))



# boxplot of proportions 

df_long2 <- species_data1 %>% select(Species,h2,width_competition,competition_strength,trait_lag, change_in_mean_trait,
                                     Degree,proportion_persisted,proportion_rescued,proportion_growth,
                                     proportion_stasis,proportion_declined, proportion_decline_stasis,minimum_path_generalist) %>% 
  tidyr::pivot_longer(cols=c(proportion_persisted,proportion_rescued,proportion_growth,
                      proportion_stasis,proportion_declined,proportion_decline_stasis),
               names_to="category", values_to = "Value") 

df_long2$Value[df_long2$Value > 1] <- 1

library(plyr)
df_long2$h2<- revalue(df_long2$h2, c("no_evo"= "No evolution", "evo_trait" = "trait evolution", "evo_trait_and_variance"="evolution trait and variance" ))
df_long2$category<- revalue(df_long2$category, 
                            c("proportion_persisted"=" Survived",
                              "proportion_declined"= "Declined", 
                              "proportion_stasis" = "Stasis", 
                              "proportion_rescued"="Rescued",
                              "proportion_growth"="Growth",
                              "proportion_decline_stasis" = "Decline-stasis"
                            ))



(bx<-df_long2  %>% 
  ggplot( aes(x = category, y = Value, color = factor(h2))) +
  geom_boxplot(alpha=0.1) + 
  geom_point(aes(colour = factor(h2)),alpha=1,size=1,
             position = position_jitterdodge(0.05))+
  ylim(c(0,1))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +  # Facet grid by h2 and width_competition
  labs(x = "", y = "Proportion", 
       size = "Frequency", title="A", color = "Model") +
  scale_y_continuous(breaks = c(0, 1)) +  # Ensure y-axis only shows 0 and 1
  theme_classic() +
  theme(legend.position = "bottom"))


a0<-ggpubr::ggarrange(a1,a2,a3,a4, nrow=2,ncol=2, common.legend = T,legend = "bottom")

############# Figure 4 ##############
library(ggpubr)
ggpubr::ggarrange(bx,a0,nrow=2,ncol=1)









###############---------Figure S8---------##################
# Create the dot plot with ggplot2


df_long$category<- plyr::revalue(df_long$category, c("decline"="decline",
                                                     "decline_stasis"="decline_stasis",
                                                     "evorescue"="rescue",
                                                     "growth"="growth",
                                                     "stasis"="stasis"))

n1<-df_long  %>% filter(h2 == "no_evo") %>% 
  ggplot( aes(x = trait_lag, y = Value, color = h2)) +
  geom_point(size=1.5, alpha=0.25) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.2 , alpha=0.15) +
  scale_color_manual(values = c( "#E7B800"))+
  # Size of dots corresponds to the count
  facet_wrap(.~ category,ncol=5,nrow=1)+ # Facet grid by h2 and width_competition
  labs(x = "Initial trait-lag", y = "Chance of a metric", title="A",size = "Frequency", color = "model") +
  scale_y_continuous(breaks = c(0, 1)) +  # Ensure y-axis only shows 0 and 1
  theme_classic() +
  theme(legend.position = "bottom")

n2<-df_long  %>% filter(h2 == "no_evo") %>% 
  ggplot( aes(x = Degree, y = Value, color = h2)) +
  geom_point(size=1.5, alpha=0.25) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = T, size=1.2 , alpha=0.15) +
  scale_color_manual(values = c( "#E7B800"))+
  # Size of dots corresponds to the count
  facet_wrap(.~ category,ncol=5,nrow=1)+ # Facet grid by h2 and width_competition
  labs(x = "Species degree", y = "Chance of a metric", title="B",size = "Frequency", color = "model") +
  scale_y_continuous(breaks = c(0, 1)) +  # Ensure y-axis only shows 0 and 1
  theme_classic() +
  theme(legend.position = "bottom")


#figure S8--------
ggpubr::ggarrange(n1,n2)


########################## Figure S2 ###################################
(s2<-df_long2 %>%
    ggplot( aes(x = Degree, y = abs(change_in_mean_trait), color = h2)) +
    geom_point(size=3.5, alpha=0.15) + 
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + # Facet grid by h2 and width_competition
    labs(x = "Degree", title="",
         y = "Absolute change in mean trait", size = "Frequency", color = "model") + # Ensure y-axis only shows 0 and 1
    theme_classic() +
    theme(legend.position = "bottom"))


############################## Figure S3 ##################################
(s3<-df_summary %>% 
    ggplot( aes(x = change_in_mean_trait, y = Value, color = h2)) +
    geom_point(size=1.5, alpha=0.25) + 
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.2 , alpha=0.15) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    # Size of dots corresponds to the count
    facet_wrap(.~category,ncol=5,nrow=1) +  # Facet grid by h2 and width_competition
    labs(x = "Change in mean trait", title="",
         y = "Chances of a metric", size = "Frequency", color = "model") +
    scale_y_continuous(breaks = c(0, 1)) +# Ensure y-axis only shows 0 and 1
    theme_classic() +
    theme(legend.position = "bottom"))

