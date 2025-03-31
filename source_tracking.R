library(tidyverse)
library(ggpubr)

meta_bact <- read.table(file = "16S/metadata_16S.csv", header=T, sep=';', row.names = 1) %>% 
  rownames_to_column() %>% 
  select(rowname, niche) %>% 
  rename('sink'='niche') %>% 
  rename('SampleID'='rowname')

bact_healthy <- read.table('source_tracking/16S_healthy_mixing_proportions.txt', header = T,
                           sep = '\t')

df_healthy <- bact_healthy %>% 
  merge(., meta_bact) %>% 
  select(!SampleID) %>% 
  pivot_longer(cols = c('bulk.soil':'Unknown'), names_to = 'source', values_to = 'percent') %>% 
  mutate(health = 'healthy')

df_healthy$source <- gsub('bulk.soil', 'bulk soil', df_healthy$source)

bact_unhealthy <- read.table('source_tracking/16S_unhealthy_mixing_proportions.txt', header = T,
                           sep = '\t')

df_unhealthy <- bact_unhealthy %>% 
  merge(., meta_bact) %>% 
  select(!SampleID) %>% 
  pivot_longer(cols = c('bulk.soil':'Unknown'), names_to = 'source', values_to = 'percent') %>% 
  mutate(health = 'unhealthy')

df_unhealthy$source <- gsub('bulk.soil', 'bulk soil', df_unhealthy$source)

df_bact <- bind_rows(df_healthy, df_unhealthy) %>% 
  mutate(micro='bacteria')



##### FUNGI #####
meta_f <- read.table(file = "ITS/metadata_ITS.csv", header=T, sep=';', row.names = 1) %>% 
  rownames_to_column() %>% 
  select(rowname, niche) %>% 
  rename('sink'='niche') %>% 
  rename('SampleID'='rowname')

f_healthy <- read.table('source_tracking/its_healthy_mixing_proportions.txt', header = T,
                           sep = '\t')

f_df_healthy <- f_healthy %>% 
  merge(., meta_f) %>% 
  select(!SampleID) %>% 
  pivot_longer(cols = c('bulk.soil':'Unknown'), names_to = 'source', values_to = 'percent') %>% 
  mutate(health = 'healthy')

f_df_healthy$source <- gsub('bulk.soil', 'bulk soil', f_df_healthy$source)

f_unhealthy <- read.table('source_tracking/its_unhealthy_mixing_proportions.txt', header = T,
                             sep = '\t')

f_df_unhealthy <- f_unhealthy %>% 
  merge(., meta_f) %>% 
  select(!SampleID) %>% 
  pivot_longer(cols = c('bulk.soil':'Unknown'), names_to = 'source', values_to = 'percent') %>% 
  mutate(health = 'unhealthy')

f_df_unhealthy$source <- gsub('bulk.soil', 'bulk soil', f_df_unhealthy$source)

df_f <- bind_rows(f_df_healthy, f_df_unhealthy) %>% 
  mutate(micro='fungi')

df_big <- bind_rows(df_bact, df_f)

prop<-df_big %>% 
  ggplot(aes(x=source, y= percent, fill = health))+
  stat_summary(position = position_dodge(width = 0.25),
               fun.data='mean_cl_boot', linewidth=1.1, size=1, shape=21, color='black')+
  facet_grid(micro~sink)+
  ggpubr::stat_compare_means(label = "p.signif",hide.ns = TRUE, label.y = 0.75, size=8)+
  scale_y_continuous(labels=scales::percent_format(scale=100))+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size=14),
        legend.position = 'none')
ggsave('plots/source_tracking.png', height = 5, width =12)


"""
df_big %>% 
  ggplot(aes(x=source, y= percent, fill = micro))+
  stat_summary()+
  facet_grid(health~sink)+
  ggpubr::stat_compare_means(label = 'p.signif'')+
  scale_y_continuous(labels=scales::percent_format(scale=100))+
  scale_fill_manual(values=c(bacteria='darkblue', fungi = 'darkred'))+
  theme_bw()
"""
