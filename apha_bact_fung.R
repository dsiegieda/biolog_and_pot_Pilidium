library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(dplyr)
set.seed(1233323)
library(agricolae)
library(microViz)


############## BACTERIA
ps<-readRDS('16S/ps.16S.RDS')

#choose specific measures if all at once won't work
alpha.diversity <- estimate_richness(ps, measures =  c('Simpson', 'Shannon'))
alpha.diversity$ENS <- exp(alpha.diversity$Shannon)

#add those alpha diverstiy measure to the dataframe
data_b <- cbind(sample_data(ps), alpha.diversity)

labels <- c('one'='0 days', 'two'='9 days')


p_values <- data_b %>% 
  group_by(sampling,niche) %>%
  do({
    wilcox_health <- wilcox_test(ENS ~ health, data = .)
    
    # Print details for each niche for debugging
    print(wilcox_health)
    
    tibble(
      niche = unique(.$niche),
      wilcox_health_p = pvalue(wilcox_health) 
    )
  })



bact<-data_b %>%
  ggplot(aes(x=sampling, y=ENS, fill=health))+
  geom_boxplot(width=0.35)+
  #geom_jitter(width = 0.5)+
  facet_grid(~niche)+
  stat_compare_means(aes(group=sampling), label.y = 1700)+
  stat_compare_means(aes(group=health),label = "p.format", label.y = 1500)+
  labs(x=NULL,
       y='ENS')+
  scale_y_continuous()+
  #ggtitle('alpha diversity for 16S')+
  theme_bw()+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  scale_x_discrete(labels=c('one'= '0 days', 'two'= '9 days'))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

bact + 
  geom_text(data = p_values, 
            aes(x = 1.5, y = 1700, label = paste0("Health p = ", round(wilcox_health_p, 4))),
            inherit.aes = FALSE, 
            size = 4, 
            color = "black") +
  geom_text(data = p_values, 
            aes(x = 2.5, y = 1700, label = paste0("Sampling p = ", round(wilcox_sampling_p, 4))),
            inherit.aes = FALSE, 
            size = 4, 
            color = "black") 

############## FUNGI ###########
ps<-readRDS('ITS/ps.ITS.RDS')

ps<- ps_filter(ps, niche!='shoots')

#choose specific measures if all at once won't work
alpha.diversity <- estimate_richness(ps, measures =  c('Simpson', 'Shannon'))
alpha.diversity$ENS <- exp(alpha.diversity$Shannon)

#add those alpha diverstiy measure to the dataframe
data_f <- cbind(sample_data(ps), alpha.diversity)



fungi<-data_f %>%
  ggplot(aes(x=sampling, y=ENS, fill=health))+
  geom_boxplot(width=0.35)+
  #geom_jitter(width = 0.5)+
  facet_grid(~niche)+
  stat_compare_means(aes(group=sampling), label.y = 110)+
  stat_compare_means(aes(group=health),label = "p.format", label.y = 100)+
  labs(x=NULL,
       y='ENS')+
  scale_y_continuous()+
  #ggtitle('alpha diversity for ITS')+
  theme_bw()+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  scale_x_discrete(labels=c('one'= '0 days', 'two'= '9 days'))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

ggarrange(bact, fungi, common.legend = T, labels = c('A','B'))
ggsave('plots/alpha_bact_fung2.png', width = 13, height = 4.5, bg='white')




data_f$niche <- as.factor(data_f$niche)
data_f$sampling <- as.factor(data_f$sampling)
data_f$health   <- as.factor(data_f$health)


data_b$niche <- as.factor(data_b$niche)
data_b$sampling <- as.factor(data_b$sampling)
data_b$health   <- as.factor(data_b$health)
str(data_b)


### BETWEEN HEALTH 
data_f %>%
  group_by(sampling, niche) %>%
  do({
    wilcox_health <- wilcox_test(ENS ~ health, data = .)
    
    # Print details for each sampling for debugging
    print(wilcox_health)
    
    tibble(
      niche = unique(.$niche),
      wilcox_health_p = pvalue(wilcox_health)
          )
  })


data_b %>%
  group_by(sampling,niche) %>%
  do({
    wilcox_health <- wilcox_test(ENS ~ health, data = .)
    
    # Print details for each niche for debugging
    print(wilcox_health)
    
    tibble(
      niche = unique(.$niche),
      wilcox_health_p = pvalue(wilcox_health) 
    )
  })

##### NOPE ############

### BETWEEN SAMPLINGS

data_b %>%
  group_by(niche) %>%
  summarise(
    wilcox_health_p = wilcox.test(ENS ~ health)$p.value,
    .groups = 'drop'
  )


data_f %>%
  group_by(niche) %>%
  summarise(
    wilcox_health_p = wilcox.test(ENS ~ health)$p.value,
    .groups = 'drop'
  )
