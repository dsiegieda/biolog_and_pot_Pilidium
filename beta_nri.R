"""The code was developed by Dominika Siegieda, PhD, Institute of Agrophysics, Polish Academy of Sciences. This research was funded in whole by National Science Centre, Poland, contract number:2022/45/N/NZ9/02089"""

library(microeco)
library(file2meco)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(jcolors)

runs = 2000

########## SOIL ##########
ps<-readRDS('16S/ps.16S_rar.RDS')
ps <- phyloseq::subset_samples(ps, niche=="bulk soil")
ps <- phyloseq::subset_samples(ps, sampling=="two")
meco <- phyloseq2meco(ps)

meco$cal_betadiv(wei_unifrac = TRUE)

t1 <- trans_nullmodel$new(meco, filter_thres = 0.0005)

t1$cal_ses_betampd(runs = runs, abundance.weighted = TRUE)


# add betaNRI matrix to beta_diversity list
t1$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd


# create trans_beta class, use measure "betaNRI"
t2 <- trans_beta$new(dataset = t1, group = "health", measure = "betaNRI")

# transform the distance for each group
t2$cal_group_distance()

# see the help document for more methods, e.g. "anova" and "KW_dunn"
t2$cal_group_distance_diff(method = "wilcox")

# plot the results
g1 <- t2$plot_group_distance()
g1 + 
  geom_hline(yintercept = -2, linetype = 2) + 
  geom_hline(yintercept = 2, linetype = 2)+
  scale_color_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  theme()


distances_bs<-t2$res_group_distance %>% 
  mutate(niche='bulk soil')


##### RHIZO ####
ps<-readRDS('16S/ps.16S_rar.RDS')
ps <- phyloseq::subset_samples(ps, niche=="rhizosphere")
ps <- phyloseq::subset_samples(ps, sampling=="two")
meco <- phyloseq2meco(ps)

meco$cal_betadiv(wei_unifrac = TRUE)

t1 <- trans_nullmodel$new(meco, filter_thres = 0.0005)

t1$cal_ses_betampd(runs = runs, abundance.weighted = TRUE)


# add betaNRI matrix to beta_diversity list
t1$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd


# create trans_beta class, use measure "betaNRI"
t2 <- trans_beta$new(dataset = t1, group = "health", measure = "betaNRI")

# transform the distance for each group
t2$cal_group_distance()

# see the help document for more methods, e.g. "anova" and "KW_dunn"
t2$cal_group_distance_diff(method = "wilcox")

# plot the results
g1 <- t2$plot_group_distance()
g1 + 
  geom_hline(yintercept = -2, linetype = 2) + 
  geom_hline(yintercept = 2, linetype = 2)+
  scale_color_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  theme()


distances_rhizo<-t2$res_group_distance %>% 
  mutate(niche='rhizosphere')



##### ROOTS ####
ps<-readRDS('16S/ps.16S_rar.RDS')
ps <- phyloseq::subset_samples(ps, niche=="roots")
ps <- phyloseq::subset_samples(ps, sampling=="two")
meco <- phyloseq2meco(ps)

meco$cal_betadiv(wei_unifrac = TRUE)

t1 <- trans_nullmodel$new(meco, filter_thres = 0.0005)

t1$cal_ses_betampd(runs = runs, abundance.weighted = TRUE)


# add betaNRI matrix to beta_diversity list
t1$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd


# create trans_beta class, use measure "betaNRI"
t2 <- trans_beta$new(dataset = t1, group = "health", measure = "betaNRI")

# transform the distance for each group
t2$cal_group_distance()

# see the help document for more methods, e.g. "anova" and "KW_dunn"
t2$cal_group_distance_diff(method = "wilcox")

# plot the results
g1 <- t2$plot_group_distance()
g1 + 
  geom_hline(yintercept = -2, linetype = 2) + 
  geom_hline(yintercept = 2, linetype = 2)+
  scale_color_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  theme()


distances_roots<-t2$res_group_distance %>% 
  mutate(niche='roots')


##### ALL ####

df <- distances_bs %>% 
  rbind(distances_rhizo, distances_roots)


nri<-df %>% 
  ggplot(aes(y=Value, x= niche, fill=health))+
  geom_violin(scale = "count", adjust = .75, position = position_dodge(width = 0.5))+
  stat_compare_means(method='t.test', label = "p.signif", hide.ns=T, label.y=0.8)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))+
  geom_hline(yintercept = -2, linetype = 2) + 
  geom_hline(yintercept = 2, linetype = 2) +
  theme_bw()+
  labs(y='BetaNRI',
       x='')+
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=16),
        plot.margin = unit(c(1,1,0,2), "lines"))

ggarrange(prop, nri, ncol = 2,labels = c("A", "B"), font.label = list(size=18), widths = c(3,1))
ggsave('plots/source_tracking_and_beta.png', height = 5, width =14)

