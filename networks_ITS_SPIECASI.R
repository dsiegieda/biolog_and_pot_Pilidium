"""The code was developed by Dominika Siegieda, PhD, Institute of Agrophysics, Polish Academy of Sciences. This research was funded in whole by National Science Centre, Poland, contract number:2022/45/N/NZ9/02089"""


### zwiększenie limitu pamięci dla LINUXA    ulimit::memory_limit(200000)



library(phyloseq)
library(NetCoMi)
library(ulimit)


############ ITS ######
ps<-readRDS('ITS/ps.ITS_rar.RDS')
ps <- phyloseq::subset_samples(ps, sampling=="two")
ps <- phyloseq::subset_samples(ps, niche=="rhizosphere")
ps <- tax_glom(ps, taxrank = 'Genus')


# Taxonomic table
taxtab <- as(tax_table(ps), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
ps <- renameTaxa(ps, 
                 pat = "<name>", 
                 substPat = "<name>_<subst_name>(<subst_R>)",
                 numDupli = "Genus")

ps_c <- phyloseq::subset_samples(ps, health=="healthy")
ps_un <- phyloseq::subset_samples(ps, health=="unhealthy")

n_yes <- phyloseq::nsamples(ps_c)


reps <- 30 # make it higher


net_season <- netConstruct(data = ps_c, 
                           data2 = ps_un,
                           taxRank = "Genus",
                           measure = "spieceasi",
                           measurePar = list(method = "glasso",
                                             pulsar.params = list(rep.num = reps),
                                             nlambda=20,
                                             lambda.min.ratio=1e-4,
                                             symBetaMode = "ave"),
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 200),
                           sparsMethod = "threshold",
                           zeroMethod = "pseudo",
                           verbose = 2,
                           seed=12)


props_season <- netAnalyze(net_season, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubQuant = 0.95,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)

summary(props_season)

plotHeat(net_season$assoMat1, textUpp = "none", textLow = "none")
edges<-net_season$edgelist1

png('plots/its.png', width = 2800, height = 1100, units = "px", pointsize = 12)

plot<-plot(props_season, 
           sameLayout = TRUE, 
           repulsion = 0.95,
           layoutGroup = 'union',
           rmSingles = "none", 
           nodeTransp = 0,
           nodeColor = "cluster",
           nodeSize = "clr", 
           labelScale = FALSE,
           cexNodes = 0.5, 
           cexLabels = 0,
           cexHubLabels = 0,
           cexTitle = 0.5,
           edgeTranspLow = 0,
           edgeTranspHigh = 0,
           edgeInvisFilter = 'threshold',
           edgeInvisPar = 0.4,
           groupNames = c("", ""),
           hubBorderCol  = "gray40")

dev.off()

png('plots/its_labels.png', width = 2800, height = 1100, units = "px", pointsize = 12)

plot<-plot(props_season, 
           sameLayout = TRUE, 
           repulsion = 0.95,
           layoutGroup = 'union',
           rmSingles = "none", 
           nodeTransp = 0,
           nodeColor = "cluster",
           nodeSize = "clr", 
           labelScale = FALSE,
           cexNodes = 0.5, 
           cexLabels = 2,
           cexHubLabels = 0,
           cexTitle = 0.5,
           edgeTranspLow = 0,
           edgeTranspHigh = 0,
           edgeInvisFilter = 'threshold',
           edgeInvisPar = 0.4,
           groupNames = c("", ""),
           hubBorderCol  = "gray40")

dev.off()


comp_season <- netCompare(props_season, 
                          permTest = FALSE, 
                          verbose = FALSE,
                          seed = 123456)

summary(comp_season, 
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 10)




