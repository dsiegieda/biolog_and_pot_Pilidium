library(RColorBrewer)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(lme4)
library(lmerTest) 
library(emmeans)
library(broom)


### FUNG  

"""
ps<-readRDS('ITS/ps.ITS.RDS')

ps<- ps_filter(ps, niche!='shoots')

#choose specific measures if all at once won't work
alpha.diversity <- estimate_richness(ps, measures =  c('Simpson', 'Shannon'))
alpha.diversity$ENS <- exp(alpha.diversity$Shannon)

#add those alpha diverstiy measure to the dataframe
data_f <- cbind(sample_data(ps), alpha.diversity)%>% 
  rownames_to_column()

writexl::write_xlsx(data_f, 'alpha_f.xlsx')
"""

data_f<- readxl::read_xlsx('alpha_f.xlsx') %>% 
  column_to_rownames()


results <- data_f %>%
  group_by(niche, sampling) %>%
  group_modify(~ {
    model <- lmer(ENS ~ health + (1 | pot), data = .x)
    emm <- emmeans(model, pairwise ~ health)
    tidy(emm$contrasts)  # extract comparison results
  })


p_labels <- results %>%
  mutate(
    p_label = sprintf("p = %.3f", p.value), 
    health = "healthy", 
    ENS = max(data_f$ENS, na.rm = TRUE)*1.15)


labels <- c('one'='0 days', 'two'='9 days')



plot_f<-data_f %>%
  ggplot(aes(x=health, y=ENS, fill=health))+
  geom_boxplot(width=0.35)+
  facet_grid(sampling ~niche,
             labeller = labeller(sampling = labels))+
  geom_text(data=p_labels, aes(label=p_label))+
  scale_y_continuous()+
  #ggtitle('alpha diversity for ITS')+
  theme_bw()+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  scale_x_discrete(labels=c('one'= '0 days', 'two'= '9 days'))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))




############## BACT ###########
"""
ps<-readRDS('16S/ps.16S.RDS')

#choose specific measures if all at once won't work
alpha.diversity <- estimate_richness(ps, measures =  c('Simpson', 'Shannon'))
alpha.diversity$ENS <- exp(alpha.diversity$Shannon)

#add those alpha diverstiy measure to the dataframe
data_b <- cbind(sample_data(ps), alpha.diversity) %>% 
  rownames_to_column()


writexl::write_xlsx(data_b, 'alpha_b.xlsx')
"""
data_b<- readxl::read_xlsx('alpha_b.xlsx') %>% 
  column_to_rownames()


results <- data_b %>%
  group_by(niche, sampling) %>%
  group_modify(~ {
    model <- lmer(ENS ~ health + (1 | pot), data = .x)
    emm <- emmeans(model, pairwise ~ health)
    tidy(emm$contrasts)  # extract comparison results
  })



facet_max <- data_b %>%
  group_by(sampling) %>%
  summarize(ENS = max(ENS, na.rm = TRUE))

# Join this with your p_labels to position labels correctly
p_labels_corrected <- results %>%
  left_join(facet_max, by = c( "sampling")) %>%
  mutate(
    p_label = ifelse(p.value < 0.001, "p < 0.001", sprintf("p = %.3f", p.value)),
    health = "healthy",
    # Position label at 90% of each facet's max value
    y_pos = ENS * 0.9
  )


plot_b<-data_b %>%
  ggplot(aes(x=health, y=ENS, fill=health))+
  geom_boxplot(width=0.35)+
  facet_grid(sampling ~niche,
             labeller = labeller(sampling = labels),
             scales='free')+
  geom_text(data=p_labels_corrected, aes(label=p_label))+
  scale_y_continuous()+
  #ggtitle('alpha diversity for ITS')+
  theme_bw()+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  scale_x_discrete(labels=c('one'= '0 days', 'two'= '9 days'))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


ggarrange(plot_b, plot_f, common.legend = T, labels='AUTO')
ggsave('plots/alpha_16s_ITS_lmer.png', height = 5, width = 12, bg='white')
