require(phylotate)
require(phangorn)
require(rwty)
require(phytools)
require(tidyverse)
require(geiger)
require(beanplot)
require(FossilSim)
require(ggtree)
library(devtools)
require(RevGadgets)

windwaker <- c("#19647E", "#28AFB0", "#DDCECD", "#37392E")


##############################
# Fossil Sampling Experiment #
##############################

setwd("~/Desktop/cloud_run")

#treatments <- as.factor(treatments)
#true_trees <- c("output/all.tree", "output/rnd.tree", "output/brw.tree", "output/pd.tree")
run_names <- seq(1, 50, 1)
run_names <- paste0("Run_", run_names)

sym.dif <- vector()
#bra.dif <- vector()
#quad.dif <- vector()
sym.dif.extinct <- vector()
#bra.dif.extinct <- vector()
#quad.dif.extinct <- vector()
sym.dif.extant <- vector()



for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatments <- c("all", "rnd", "brw", "pd")
  for(j in 1:length(treatments)){
    setwd(treatments[j])
    mcc <- read.mrbayes(file = "output/sim.mcc.tre")
    true_tree <- read.tree(paste0("output/", treatments[j], ".tree"))
    fossils <- getExtinct(true_tree)
    
    sym.dif <- c(sym.dif, RF.dist(mcc@phylo, true_tree, normalize = T)[[1]])
    #bra.dif <- c(bra.dif, treedist(mcc@phylo, true_tree)[[2]])
    #quad.dif <- c(quad.dif, treedist(mcc@phylo, true_tree)[[4]])
    fossils <- getExtinct(true_tree)
    
    the_living <- getExtant(true_tree)
    mcc_extinct <- drop.tip(mcc@phylo, the_living)
    true_tree_extinct <- drop.tip(true_tree, the_living)
    
    sym.dif.extinct <- c(sym.dif.extinct, RF.dist(mcc_extinct, true_tree_extinct, normalize = T)[[1]])
    #bra.dif.extinct <- c(bra.dif.extinct, treedist(mcc_extinct, true_tree_extinct)[[2]])
    #quad.dif.extinct <- c(quad.dif.extinct, treedist(mcc_extinct, true_tree_extinct)[[4]])
    
    the_dead <- getExtinct(true_tree)
    mcc_extant <- drop.tip(mcc@phylo, the_dead)
    true_tree_extant <- drop.tip(true_tree, the_dead)
    
    sym.dif.extant <- c(sym.dif.extant, RF.dist(mcc_extant, true_tree_extant, normalize = T)[[1]])
    
    print(getwd())
    setwd("..")
    }
setwd("..")
}

treatments <- rep(c("Cmp", "Rnd", "Trt", "Div"), 50)

run_names_full <- vector()

distances <- data_frame(sym.dif, sym.dif.extant, sym.dif.extinct, treatments)

save(distances, file = "distances.Rdata")
setwd("~/Desktop/cloud_run")
load("distances.Rdata")
distances$treatments <- as.factor(distances$treatments)


kw.sym.dif <- kruskal.test(sym.dif ~ treatments, data = distances)
kw.sym.dif
require(FSA)
Dunn <- dunnTest(sym.dif ~ treatments, data = distances, method = "bonferroni")
Dunn

kw.sym.dif <- kruskal.test(sym.dif.extant ~ treatments, data = distances)
kw.sym.dif
require(FSA)
Dunn <- dunnTest(sym.dif.extant ~ treatments, data = distances, method = "bonferroni")
Dunn

kw.sym.dif <- kruskal.test(sym.dif.extinct ~ treatments, data = distances)
kw.sym.dif
require(FSA)
Dunn <- dunnTest(sym.dif.extinct ~ treatments, data = distances, method = "bonferroni")
Dunn


full_house <- ggplot(data = distances, aes(x = sym.dif, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div")) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Robinson-Foulds distance", limits = c(-0.1, 0.8))+
  ggtitle("a")

the_living <- ggplot(data = distances, aes(x = sym.dif.extant, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div")) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Robinson-Foulds distance", limits = c(-0.1, 0.8))+
  ggtitle("b")

the_dead <- ggplot(data = distances, aes(x = sym.dif.extinct, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div")) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Robinson-Foulds distance", limits = c(-0.1, 0.8))+
  ggtitle("c")


tree_distance <- grid.arrange(full_house, the_living, the_dead, nrow = 3)
ggsave("distance.png", plot = tree_distance, device = "png", width = 9.5, height = 12.5)


cmp <- subset(distances, distances$treatments == "Cmp")
CI95(cmp$sym.dif)
CI95(cmp$sym.dif.extant)
CI95(cmp$sym.dif.extinct)

rnd <- subset(distances, distances$treatments == "Rnd")
CI95(rnd$sym.dif)
CI95(rnd$sym.dif.extant)
CI95(rnd$sym.dif.extinct)

trt <- subset(distances, distances$treatments == "Trt")
CI95(trt$sym.dif)
CI95(trt$sym.dif.extant)
CI95(trt$sym.dif.extinct)

div <- subset(distances, distances$treatments == "Div")
CI95(div$sym.dif)
CI95(div$sym.dif.extant)
CI95(div$sym.dif.extinct)


median <- subset(parameters, parameters$treatments == "Median")
CI95(median$speciation)
CI95(median$extinction)

high <- subset(parameters, parameters$treatments == "High")
CI95(high$speciation)
CI95(high$extinction)


beanplot(sym.dif ~ treatments, data = distances,
         ylab = "symmetric distance", what = c(1,1,1,0),
         col = "darkslategray3")

dist.aov <- aov(sym.dif ~ treatments, data = distances)
summary(dist.aov)

beanplot(sym.dif.extinct ~ treatments, data = distances,
         ylab = "symmetric distance", what = c(1,1,1,0),
         col = "darkslategray3")

dist.aov.extinct <- aov(sym.dif.extinct ~ treatments, data = distances)
summary(dist.aov.extinct)

########################
# Imbalance Experiment #
########################


setwd("~/Desktop/comb_trees_100")

run_namesH <- seq(1, 20, 1)
run_namesH <- paste0("High_", run_namesH)
run_namesM <- seq(1, 20, 1)
run_namesM <- paste0("Median_", run_namesM)
run_namesL <- seq(1, 20, 1)
run_namesL <- paste0("Low_", run_namesL)
run_names <- c(run_namesL, run_namesM, run_namesH)
treatments <- c(replicate(20, "1.Low"), replicate(20, "2.Median"), replicate(20, "3.High"))

sym.dif <- vector()
bra.dif <- vector()
sym.dif.extinct <- vector()
bra.dif.extinct <- vector()
quad.dif.extinct <- vector()
quad.dif <- vector()
sym.dif.extant <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatment_name <- treatments[i]
  mcc <- read.mrbayes(file = "output/sim.mcc.tre")
  
  true_tree <- read.tree("output/true.tree")
  fossils <- getExtinct(true_tree)
  
  sym.dif <- c(sym.dif, RF.dist(mcc@phylo, true_tree, normalize = T)[[1]])
  #bra.dif <- c(bra.dif, treedist(mcc@phylo, true_tree)[[2]])
  #quad.dif <- c(quad.dif, treedist(mcc@phylo, true_tree)[[4]])
  
  fossils <- getExtinct(true_tree)

  the_living <- getExtant(true_tree)
  mcc_extinct <- drop.tip(mcc@phylo, the_living)
  true_tree_extinct <- drop.tip(true_tree, the_living)
  
  sym.dif.extinct <- c(sym.dif.extinct, RF.dist(mcc_extinct, true_tree_extinct, normalize = T)[[1]])
  #bra.dif.extinct <- c(bra.dif.extinct, treedist(mcc_extinct, true_tree_extinct)[[2]])
  #quad.dif.extinct <- c(quad.dif.extinct, treedist(mcc@phylo, true_tree)[[4]])
  
  the_dead <- getExtinct(true_tree)
  mcc_extant <- drop.tip(mcc@phylo, the_dead)
  true_tree_extant <- drop.tip(true_tree, the_dead)
  
  sym.dif.extant <- c(sym.dif.extant, RF.dist(mcc_extant, true_tree_extant, normalize = T)[[1]])
  
  
  print(getwd())
  setwd("..")
}

treatments <- c(replicate(20, "Low"), replicate(20, "Median"), replicate(20, "High"))

distances <- data_frame(run_names, sym.dif, sym.dif.extant, sym.dif.extinct, treatments)
save(distances, file = "distances.Rdata")
setwd("~/Desktop/comb_trees_100")
load("distances.Rdata")
distances$treatments <- as.factor(distances$treatments)

kw.sym.dif <- kruskal.test(sym.dif ~ treatments, data = distances)
kw.sym.dif
require(FSA)
Dunn <- dunnTest(sym.dif ~ treatments, data = distances, method = "bonferroni")
Dunn

kw.sym.dif <- kruskal.test(sym.dif.extant ~ treatments, data = distances)
kw.sym.dif
require(FSA)
Dunn <- dunnTest(sym.dif.extant ~ treatments, data = distances, method = "bonferroni")
Dunn

kw.sym.dif <- kruskal.test(sym.dif.extinct ~ treatments, data = distances)
kw.sym.dif
require(FSA)
Dunn <- dunnTest(sym.dif.extinct ~ treatments, data = distances, method = "bonferroni")
Dunn


low <- subset(distances, distances$treatments == "Low")
CI95(low$sym.dif)
CI95(low$sym.dif.extant)
CI95(low$sym.dif.extinct)

median <- subset(distances, distances$treatments == "Median")
CI95(median$sym.dif)
CI95(median$sym.dif.extant)
CI95(median$sym.dif.extinct)

high <- subset(distances, distances$treatments == "High")
CI95(high$sym.dif)
CI95(high$sym.dif.extant)
CI95(high$sym.dif.extinct)



full_house <- ggplot(data = distances, aes(x = sym.dif, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Low", "Median", "High")) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Robinson-Foulds distance", limits = c(-0.15, 0.6))+
  ggtitle("a")

the_living <- ggplot(data = distances, aes(x = sym.dif.extant, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Low", "Median", "High")) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Robinson-Foulds distance", limits = c(-0.15, 0.6))+
  ggtitle("b")

the_dead <- ggplot(data = distances, aes(x = sym.dif.extinct, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Low", "Median", "High")) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Robinson-Foulds distance", limits = c(-0.15, 0.6))+
  ggtitle("c")


tree_distance <- grid.arrange(full_house, the_living, the_dead, nrow = 3)
ggsave("distance.png", plot = tree_distance, device = "png", width = 9.5, height = 12.5)




beanplot(quad.dif ~ treatments, data = distances,
         ylab = "symmetric distance", what = c(1,1,1,0),
         col = "darkslategray3")

dist.aov <- aov(sym.dif ~ treatments, data = distances)
summary(dist.aov)

beanplot(sym.dif.extinct ~ treatments, data = distances,
         ylab = "symmetric distance", what = c(1,1,1,0),
         col = "darkslategray3")

dist.aov.extinct <- aov(sym.dif.extinct ~ treatments, data = distances)
summary(dist.aov.extinct)

#####################################
# Clade Level Extinction Experiment #
#####################################

setwd("~/Desktop/clade_extinction")

run_names <- seq(1,44,1)
run_names <- paste0("clade_extinction_", run_names)

sym.dif <- vector()
bra.dif <- vector()
quad.dif <- vector()
sym.dif.extinct <- vector()
bra.dif.extinct <- vector()
quad.dif.extinct <- vector()
clade_size <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])

  mcc <- read.mrbayes(file = "output/sim.mcc.tre")
  true_tree <- read.tree("output/true.tree")
  

  sym.dif <- c(sym.dif, RF.dist(mcc@phylo, true_tree, normalize = T)[[1]])
  bra.dif <- c(bra.dif, treedist(mcc@phylo, true_tree)[[2]])
  quad.dif <- c(quad.dif, treedist(mcc@phylo, true_tree)[[4]])
  
  fossils <- getExtinct(true_tree)
  clade_size <- c(clade_size, length(fossils))
  
  the_living <- getExtant(true_tree)
  mcc_extinct <- drop.tip(mcc@phylo, the_living)
  true_tree_extinct <- drop.tip(true_tree, the_living)
  
  sym.dif.extinct <- c(sym.dif.extinct, RF.dist(mcc_extinct, true_tree_extinct, normalize = T)[[1]])
  bra.dif.extinct <- c(bra.dif.extinct, treedist(mcc_extinct, true_tree_extinct)[[2]])
  quad.dif.extinct <- c(quad.dif.extinct, treedist(mcc_extinct, true_tree_extinct)[[4]])
  
  
  print(getwd())
  setwd("..")
}

plot(clade_size, sym.dif)
cor.test(clade_size, sym.dif)
plot(clade_size, bra.dif)
cor.test(clade_size, bra.dif)

plot(clade_size, sym.dif.extinct)
cor.test(clade_size, sym.dif.extinct)
plot(clade_size, bra.dif.extinct)
cor.test(clade_size, bra.dif.extinct)

plot(cophylo(true_tree_extinct, mcc_extinct))
