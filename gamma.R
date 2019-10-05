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
require(tidytree)
require(ggridges)
#and a function to determine the length of the tree
tree.length <- function(x) {
  max(diag(vcv(x)))
}

CI95 <- function(x){
  a <- mean(x)
  n <- length(x)
  s <- sd(x)
  error <- qnorm(0.975)*s/sqrt(n)
  left <- a-error
  right <- a+error
  c(a, left,right)
}

cynthwave <- c("#FFA9E7", "#FF84E8", "#7F2CCB", "#414361", "#2A2D43")

##############################
# Fossil Sampling Experiment #
##############################

setwd("~/Desktop/cloud_run")

run_names <- seq(1, 50, 1)
run_names <- paste0("Run_", run_names)

all_gammas <- vector()
all_ages <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatments <- c("all", "rnd", "brw", "pd")
  for(j in 1:length(treatments)){
    setwd(treatments[j])
    #we load in all the trees, however we have to drop the extinct tips because gamma
    #assumes an ultrametric tree
    mcc <- read.mrbayes(file = "output/sim.mcc.tre")
    true_tree <- read.tree(paste0("output/", treatments[j], ".tree"))
    fossils <- getExtinct(true_tree)
    mcc <- drop.tip(mcc@phylo, fossils)
    mcc <- force.ultrametric(mcc)
    true_tree <- drop.tip(true_tree, fossils)
    
    
    trees <- read.tree(file = "output/sim.trees")
    trees <- trees[2502:10001] #drop 25% as burnin
    
    gamma.true <- gammaStat(true_tree)
    age.true <- tree.length(true_tree)
    gamma.mcc <- gammaStat(mcc)
    age.mcc <- tree.length(mcc)
    gamma.post <- vector()
    age.post <- vector()
    
    
    for(i in 1:length(trees)){
      x <- trees[[i]]
      if(typeof(x) == "list") {
      x <- drop.tip(x, fossils)
      x <- force.ultrametric(x)
      y <- gammaStat(x)
      z <- tree.length(x)
      gamma.post <- c(gamma.post, y)
      age.post <- c(age.post, z)
      print(i)
      }
      else{
        print(i)
      }
    }
    #the if/else statement is there because some of the samples are corrupted - just a bit of errori handling
    
    gamma.post <- CI95(gamma.post)
    all_gammas <- rbind(all_gammas, c(gamma.true, gamma.mcc, gamma.post))
    age.post <- CI95(age.post)
    all_ages <- rbind(all_ages, c(age.true, age.mcc, age.post))
    print(getwd())
    setwd("..")
  }
  setwd("..")
}

treatments <- rep(c("Cmp", "Rnd", "Trt", "Div"), 50)

####################
#get fossil numbers#
####################

fossil_number <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatments <- c("all", "rnd", "brw", "pd")
  for(j in 1:length(treatments)){
    setwd(treatments[j])
    true_tree <- read.tree(paste0("output/", treatments[j], ".tree"))
    fossil_number <- c(fossil_number, length(getExtinct(true_tree)))
    print(getwd())
    setwd("..")
  }
  setwd("..")
}

all_gammas <- as.data.frame(all_gammas)
colnames(all_gammas) <- c("gammas_true", "gammas_mcc", "mean_post_gamma", "CI_low", "CI_high")
all_gammas$precision <- all_gammas$CI_high - all_gammas$CI_low
all_gammas$bias <- all_gammas$gammas_true - all_gammas$gammas_mcc
all_gammas$treatments <- treatments
all_gammas$treatments <- as.factor(all_gammas$treatments)
all_gammas$fossil_number <- as.integer(fossil_number)
write.csv(all_gammas, file = "gamma_scores.csv")

all_ages <- as.data.frame(all_ages)
colnames(all_ages) <- c("ages_true", "ages_mcc", "mean_post_ages", "CI_low", "CI_high")
all_ages$precision <- all_ages$CI_high - all_ages$CI_low
all_ages$bias <- all_ages$ages_true - all_ages$ages_mcc
all_ages$treatments <- treatments
all_ages$treatments <- as.factor(all_ages$treatments)
all_ages$fossil_number <- as.factor(fossil_number)
write.csv(all_ages, file = "origin_ages.csv")

save(all_gammas, file = "gamma_scores.Rdata")
save(all_ages, file = "origin_ages.Rdata")

### ANLYSES ###

setwd("~/Desktop/cloud_run")
load("gamma_scores.Rdata")
load("origin_ages.Rdata")
all_ages$treatments <- as.factor(all_ages$treatments)
all_gammas$treatments <- as.factor(all_gammas$treatments)


kw.age <- kruskal.test(bias ~ treatments, data = all_ages)
kw.age

kw.gammas <- kruskal.test(bias ~ treatments, data = all_gammas)
kw.gammas

age.b <-ggplot(data = all_ages, aes(x = bias, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), expand = c(0.01, 0)) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Origin age differences (Ma)")+
  ggtitle("a")

age.p <-ggplot(data = all_ages, aes(x = precision, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), expand = c(0.01, 0)) +
  scale_fill_manual(values = cynthwave, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  scale_x_continuous(name = "Origin age precision (Ma)")+
  theme(axis.title.y = element_blank())

gam.b <-ggplot(data = all_gammas, aes(x = bias, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), expand = c(0.01, 0)) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = expression(paste(gamma, " score differences")))+
  ggtitle("b") 

gam.p <-ggplot(data = all_gammas, aes(x = precision, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 3, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), expand = c(0.01, 0)) +
  scale_fill_manual(values = cynthwave, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  scale_x_continuous(name = expression(paste(gamma, " score precision")))+
  theme(axis.title.y = element_blank())

age_and_gamma <- grid.arrange(age.b, gam.b, nrow = 2)
ggsave("age_and_gamma_ridge.png", plot = age_and_gamma, device = "png", width = 9, height = 10)


cmp <- subset(all_gammas, all_gammas$treatments == "Cmp")
rnd <- subset(all_gammas, all_gammas$treatments == "Rnd")
trt <- subset(all_gammas, all_gammas$treatments == "Trt")
div <- subset(all_gammas, all_gammas$treatments == "Div")
CI95(cmp$bias)
CI95(rnd$bias)
CI95(trt$bias)
CI95(div$bias)


cmp <- subset(all_ages, all_ages$treatments == "Cmp")
rnd <- subset(all_ages, all_ages$treatments == "Rnd")
trt <- subset(all_ages, all_ages$treatments == "Trt")
div <- subset(all_ages, all_ages$treatments == "Div")
CI95(cmp$bias)
CI95(rnd$bias)
CI95(trt$bias)
CI95(div$bias)


aov.ages.b <- aov(bias ~ treatments, data = all_ages)
summary(aov.ages.b)

aov.ages.p <- aov(precision ~ treatments, data = all_ages)
summary(aov.ages.p)

aov.gamma.b <- aov(bias ~ treatments, data = all_gammas)
summary(aov.gamma.b)

aov.gamma.p <- aov(bias ~ treatments, data = all_gammas)
summary(aov.gamma.p)

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
treatments <- c(replicate(20, "Low"), replicate(20, "Median"), replicate(20, "High"))
#treatments <- as.factor(treatments)

all_gammas <- vector()
all_ages <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  mcc <- read.mrbayes(file = "output/sim.mcc.tre")
  true_tree <- read.tree(file = "output/true.tree")
  fossils <- getExtinct(true_tree)
  mcc <- drop.tip(mcc@phylo, fossils)
  mcc <- force.ultrametric(mcc)
  true_tree <- drop.tip(true_tree, fossils)
  
  
  trees <- read.tree(file = "output/sim.trees")
  trees <- trees[2502:10001] #drop 25% as burnin
  
  gamma.true <- gammaStat(true_tree)
  age.true <- tree.length(true_tree)
  gamma.mcc <- gammaStat(mcc)
  age.mcc <- tree.length(mcc)
  gamma.post <- vector()
  age.post <- vector()
  
  
  for(i in 1:length(trees)){
    x <- trees[[i]]
    if(typeof(x) == "list") {
      x <- drop.tip(x, fossils)
      x <- force.ultrametric(x)
      y <- gammaStat(x)
      z <- tree.length(x)
      gamma.post <- c(gamma.post, y)
      age.post <- c(age.post, z)
      print(i)
    }
    else{
      print(i)
    }
  }
  #the if/else statement is there because some of the samples are corrupted - just a bit of errori handling
  
  gamma.post <- CI95(gamma.post)
  all_gammas <- rbind(all_gammas, c(gamma.true, gamma.mcc, gamma.post))
  age.post <- CI95(age.post)
  all_ages <- rbind(all_ages, c(age.true, age.mcc, age.post))
  print(getwd())
  setwd("..")
}

####################
#get fossil numbers#
####################

fossil_number <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  true_tree <- read.tree(file = "output/true.tree")
  fossil_number <- c(fossil_number, length(getExtinct(true_tree)))
  print(getwd())
  setwd("..")
}



all_gammas <- as.data.frame(all_gammas)
colnames(all_gammas) <- c("gammas_true", "gammas_mcc", "mean_post_gamma", "CI_low", "CI_high")
all_gammas$precision <- all_gammas$CI_high - all_gammas$CI_low
all_gammas$bias <- all_gammas$gammas_true - all_gammas$gammas_mcc
all_gammas$treatments <- treatments
all_gammas$treatments <- as.factor(all_gammas$treatments)
all_gammas$fossil_number <- as.integer(fossil_number)
write.csv(all_gammas, file = "gamma_scores.csv")

all_ages <- as.data.frame(all_ages)
colnames(all_ages) <- c("ages_true", "ages_mcc", "mean_post_ages", "CI_low", "CI_high")
all_ages$precision <- all_ages$CI_high - all_ages$CI_low
all_ages$bias <- all_ages$ages_true - all_ages$ages_mcc
all_ages$treatments <- treatments
all_ages$treatments <- as.factor(all_ages$treatments)
all_ages$fossil_number <- as.factor(fossil_number)
write.csv(all_ages, file = "origin_ages.csv")

save(all_gammas, file = "gamma_scores.Rdata")
save(all_ages, file = "origin_ages.Rdata")

############
# analyses #
############


setwd("~/Desktop/comb_trees_100")
load("gamma_scores.Rdata")
load("origin_ages.Rdata")

kw.age <- kruskal.test(bias ~ treatments, data = all_ages)
kw.age

kw.gammas <- kruskal.test(bias ~ treatments, data = all_gammas)
kw.gammas


low <- subset(all_ages, all_ages$treatments == "Low")
median <- subset(all_ages, all_ages$treatments == "Median")
high <- subset(all_ages, all_ages$treatments == "High")
CI95(low$bias)
CI95(median$bias)
CI95(high$bias)

low <- subset(all_gammas, all_gammas$treatments == "Low")
median <- subset(all_gammas, all_gammas$treatments == "Median")
high <- subset(all_gammas, all_gammas$treatments == "High")
CI95(low$bias)
CI95(median$bias)
CI95(high$bias)


age.b <-ggplot(data = all_ages, aes(x = bias, y = treatments, fill = treatments)) + 
        geom_density_ridges(scale = 4, alpha = .7) +
        scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
        scale_fill_manual(values = windwaker, guide = FALSE)+
        theme_ridges(center_axis_labels = T)+
        theme(axis.title.y = element_blank())+
        scale_x_continuous(name = "Origin age differences (Ma)")+
        ggtitle("a")

age.p <-ggplot(data = all_ages, aes(x = precision, y = treatments, fill = treatments)) + 
        geom_density_ridges(scale = 4, alpha = .7) +
        scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
         scale_fill_manual(values = cynthwave, guide = FALSE)+
       theme_ridges(center_axis_labels = T)+
        theme(axis.title.y = element_blank())+
        scale_x_continuous(name = "Origin age precision (Ma)")
  

gam.b <-ggplot(data = all_gammas, aes(x = bias, y = treatments, fill = treatments)) + 
        geom_density_ridges(scale = 4, alpha = .7) +
        scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
        theme(axis.title.y = element_blank())+
        scale_x_continuous(name = expression(paste(gamma, " score differences")))+
        ggtitle("b") 

gam.p <-ggplot(data = all_gammas, aes(x = precision, y = treatments, fill = treatments)) + 
        geom_density_ridges(scale = 4, alpha = .7) +
        scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
  scale_fill_manual(values = cynthwave, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
        scale_x_continuous(name = expression(paste(gamma, " score precision")))+
        theme(axis.title.y = element_blank())

age_and_gamma <- grid.arrange(age.b, gam.b, nrow = 2)
ggsave("age_and_gamma_ridge.png", plot = age_and_gamma, device = "png", width = 9, height = 10)

aov.ages.b <- aov(bias ~ treatments, data = all_ages)
summary(aov.ages.b)

aov.gamma.b <- aov(bias ~ treatments, data = all_gammas)
summary(aov.gamma.b)

aov.ages.p <- aov(precision ~ treatments, data = all_ages)
summary(aov.ages.p)
TukeyHSD(aov.ages.p)

aov.gamma.p <- aov(precision ~ treatments, data = all_gammas)
summary(aov.gamma.p)
