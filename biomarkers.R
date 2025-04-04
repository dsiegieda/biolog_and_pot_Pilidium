"""The code was developed by Dominika Siegieda, PhD, Institute of Agrophysics, Polish Academy of Sciences. This research was funded in whole by National Science Centre, Poland, contract number:2022/45/N/NZ9/02089"""


library(phyloseq)
library(microbiomeMarker)
library(ggpubr)

ps<-readRDS('16S/ps.16S.RDS')

ps_1<- subset_samples(ps, sampling == 'two')
soil_1 <-  subset_samples(ps_1, niche == 'bulk soil')
roots_1 <-  subset_samples(ps_1, niche == 'roots')
rhizo_1 <-  subset_samples(ps_1, niche == 'rhizosphere')

ml_soil_1 <- run_sl(soil_1,
                    taxa_rank = "Genus",
                group = "health",
                nfolds = 10,
                nrepeats = 10,
                top_n = 10,
                norm = "TSS",
                method = "SVM")


p_soil<-plot_ef_bar(ml_soil_1)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))

ml_rhizo_1 <- run_sl(rhizo_1,
                     taxa_rank = "Genus",
                group = "health",
                nfolds = 10,
                nrepeats = 10,
                top_n = 10,
                norm = "none",
                method = "SVM")

p_rhizo<-plot_ef_bar(ml_rhizo_1)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))


ml_roots_1<- run_sl(roots_1,
                   group = "health",
                   taxa_rank = "Genus",
                   nfolds = 10,
                   nrepeats = 10,
                   top_n = 10,
                   norm = "none",
                   method = "SVM")

p_roots<-plot_ef_bar(ml_roots_1)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))

bact<-ggarrange(p_soil, p_rhizo, p_roots, common.legend = T, ncol=3, 
          labels = c('bulk soil', 'rhizosphere', 'roots'), hjust = -0.75,
          vjust = 0.25, legend = 'none', font.label='plain')+
  theme(plot.margin = margin(1,0.1,0.1,0.1, "cm"))


##### FUNGI #####
ps<-readRDS('ITS/ps.ITS.RDS')

ps_1<- subset_samples(ps, sampling == 'two')
soil_1 <-  subset_samples(ps_1, niche == 'bulk soil')
roots_1 <-  subset_samples(ps_1, niche == 'roots')
rhizo_1 <-  subset_samples(ps_1, niche == 'rhizosphere')

ml_soil_1 <- run_sl(soil_1,
                    taxa_rank = "Genus",
                    group = "health",
                    nfolds = 10,
                    nrepeats = 10,
                    top_n = 15,
                    norm = "TSS",
                    method = "SVM")


p_soil<-plot_ef_bar(ml_soil_1)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))

ml_rhizo_1 <- run_sl(rhizo_1,
                     taxa_rank = "Genus",
                     group = "health",
                     nfolds = 10,
                     nrepeats = 10,
                     top_n = 15,
                     norm = "none",
                     method = "SVM")

p_rhizo<-plot_ef_bar(ml_rhizo_1)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))


ml_roots_1<- run_sl(roots_1,
                    group = "health",
                    taxa_rank = "Genus",
                    nfolds = 10,
                    nrepeats = 10,
                    top_n = 15,
                    norm = "none",
                    method = "SVM")

p_roots<-plot_ef_bar(ml_roots_1)+
  scale_fill_manual(values=c(healthy='lightblue', unhealthy = 'orange'))

fung<-ggarrange(p_soil, p_rhizo, p_roots, common.legend = T, ncol=3, 
                labels = c('A) bulk soil', 'B) rhizosphere', 'C) roots'), hjust = -0.15,
                vjust = 0.05, legend = 'none', font.label='plain')+
  theme(plot.margin = margin(1,0.1,0.1,0.1, "cm"))


#ggarrange(bact, fung, ncol = 1, labels = 'AUTO')
#ggsave('plots/10biomarkers.png', height = 6, width = 12, bg='white')
ggsave('plots/10biomarkers_fungi_bigger_genus.png', height = 4, width = 10, bg='white')
