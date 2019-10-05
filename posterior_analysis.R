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
#and a function to determine the length of the tree
tree.length <- function(x) {
  max(diag(vcv(x)))
}

windwaker <- c("#19647E", "#28AFB0", "#DDCECD", "#37392E")


##############################
# Fossil Sampling Experiment #
##############################

setwd("~/Desktop/cloud_run")

#treatments <- as.factor(treatments)
#true_trees <- c("output/all.tree", "output/rnd.tree", "output/brw.tree", "output/pd.tree")
run_names <- seq(1, 50, 1)
run_names <- paste0("Run_", run_names)

all_post <- vector()
mean_post <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatments <- c("all", "rnd", "brw", "pd")
  for(j in 1:length(treatments)){
  setwd(treatments[j])
  mcc <- read.mrbayes(file = "output/sim.mcc.tre")
  true <- read.tree(paste0("output/", treatments[j], ".tree"))
  
  posteriors <- data.frame(mcc@data$node, mcc@data$index, mcc@data$posterior)
  confidence <- mcc@data$age_0.95_HPD
  high <- vector()
  low <- vector()
  age <- vector()
  for(i in 1:length(confidence)){
    x <- confidence[[i]]
    low <- c(low, x[1])
    high <- c(high, x[2])
    age <- c(age, mean(c(x[1], x[2])))
  }
  
  posteriors$age <- age
  posteriors$low95 <- low
  posteriors$high95 <- high
  
  nodes <- mcc@phylo$edge[,1]
  nodes <- unique(nodes)
  
  posteriors <- posteriors[which(posteriors$mcc.data.node %in% nodes),]
  
  fossils <- getExtinct(true_tree)
  
  fossils_per_node <- function(tree, nodes, fossils){
    desc <- Descendants(tree, nodes, type = "tips")
    desc <- desc[[1]]
    desc <- true_tree$tip.label[desc]
    desc <- desc %in% fossils == TRUE
    return(length(which(desc == TRUE)))
  }
  
  posteriors$mcc.data.node <- as.numeric(posteriors$mcc.data.node)
  posteriors$mcc.data.posterior <- as.numeric(posteriors$mcc.data.posterior)
  
  nodes <- posteriors$mcc.data.node
  fossil_descendents <- vector()
  
  for(i in 1:length(nodes)){
    fossil_descendents <- c(fossil_descendents, fossils_per_node(true_tree, nodes = nodes[i], fossils = fossils))
  }
  
  
  posteriors$fossil_descendents <- fossil_descendents
  
  tree <- mcc@phylo
  n <- length(tree$tip.label)
  ee<-setNames(tree$edge[sapply(1:n,function(x,y)   
    which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ff <- which(names(ee) %in% fossils)
  
  gg <- ee[ff]
  names(gg) <- NULL
  
  posteriors$fossil_associated <-posteriors$mcc.data.node %in% gg
  
  treatment_name <- treatments[j]
  
  posteriors$treatments <- treatment_name
  #posteriors$run_names <- run_names[j]
  
  all_post <- rbind(all_post, posteriors)
  
  mean_post <- c(mean_post, mean(posteriors$mcc.data.posterior))
  
  write.csv(posteriors, file = "posteriors.csv")
  print(getwd())
  setwd("..")
  }
  setwd("..")
}

all_post$treatments <- as.factor(all_post$treatments)
levels(all_post$treatments) <- c("Cmp","Trt","Div", "Rnd")

save(all_post, file = "posterior.Rdata")

### ANLYSES ###

setwd("~/Desktop/cloud_run")
load("posterior.Rdata")

posterior_support_fos <- ggplot(data = all_post, aes(x = mcc.data.posterior, y = treatments, fill = treatments)) + 
        geom_density_ridges(scale = 4, alpha = .7) +
        scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), expand = c(0.01, 0)) +
        scale_fill_manual(values = windwaker, guide = FALSE)+
        theme_ridges(center_axis_labels = T)+
        theme(axis.title.y = element_blank())+
        scale_x_continuous(name = "Posterior node support", limits = c(0,1))

ggsave("posterior_support.png", plot = posterior_support_fos, device = "png", width = 9.5, height = 6)


cmp <- subset(all_post, all_post$treatments == "Cmp")
rnd <- subset(all_post, all_post$treatments == "Rnd")
trt <- subset(all_post, all_post$treatments == "Trt")
div <- subset(all_post, all_post$treatments == "Div")

CI95(cmp$mcc.data.posterior)
CI95(rnd$mcc.data.posterior)
CI95(trt$mcc.data.posterior)
CI95(div$mcc.data.posterior)

all_post$fossil_associated <- as.factor(all_post$fossil_associated)

#with all treatment
post.fossil.aov <- aov(mcc.data.posterior ~ treatments + fossil_associated, 
                       data = all_post)
summary(post.fossil.aov)

post.fossil.kw <- kruskal.test(mcc.data.posterior ~ treatments, 
                               data = all_post)
post.fossil.kw

#post.fossil.kw <- kruskal.test(mcc.data.posterior ~ fossil_associated, 
                               data = all_post)

#post.fossil.kw

plot(all_post$fossil_descendents, all_post$mcc.data.posterior)
cor.test(all_post$fossil_descendents, all_post$mcc.data.posterior)

TukeyHSD(post.fossil.aov)
#without all treatment

all_post <- all_post[which(all_post$treatments != "all"),]
post.fossil.aov <- aov(mcc.data.posterior ~ treatments + fossil_associated + age, 
                       data = all_post)
summary(post.fossil.aov)


beanplot(mcc.data.posterior ~ treatments, data = all_post,
         ylab = "posterior support", what = c(1,1,1,0),
         col = "darkslategray3")

#without posts = 1 
all_post <- all_post[which(all_post$mcc.data.posterior != 1),]
post.fossil.aov <- aov(mcc.data.posterior ~ treatments + fossil_associated + age, 
                       data = all_post)
summary(post.fossil.aov)


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

all_post <- vector()
mean_post <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatment_name <- treatments[i]
  mcc <- read.mrbayes(file = "output/sim.mcc.tre")

  posteriors <- data.frame(mcc@data$node, mcc@data$index, mcc@data$posterior)
  confidence <- mcc@data$age_0.95_HPD
  high <- vector()
  low <- vector()
  age <- vector()
    for(i in 1:length(confidence)){
      x <- confidence[[i]]
      low <- c(low, x[1])
      high <- c(high, x[2])
      age <- c(age, mean(c(x[1], x[2]))) #this is wrong
      
      ######################
      #this is the solution#
      ######################
      
      #tree_text <- mcc@treetext
      #tree_text <- read.tree(text = tree_text)
      #heights <- nodeHeights(tree_text)
      #heights <- unique(heights[,1])
      #ages <- tree.length(tree_text) - heights
      
      #mcc <- as_data_frame(mcc)
      #mcc <- mcc %>% arrange(branch.length)
      #mcc[(unique(mcc$branch.length)),]
      #tree_text_df <- as_data_frame(tree_text)
      #tree_text_df <- tree_text_df %>% arrange(branch.length)
      
      #you need to rearrange so branch lengths in same order as tree_text
      #then you can give it the nodes of tree_text
      
    }
  
  posteriors$age <- age
  posteriors$low95 <- low
  posteriors$high95 <- high

  nodes <- mcc@phylo$edge[,1]
  nodes <- unique(nodes)

  posteriors <- posteriors[which(posteriors$mcc.data.node %in% nodes),]

  true_tree <- read.tree(file = "output/true.tree")
  fossils <- getExtinct(true_tree)

  fossils_per_node <- function(tree, nodes, fossils){
    desc <- Descendants(tree, nodes, type = "tips")
    desc <- desc[[1]]
    desc <- true_tree$tip.label[desc]
    desc <- desc %in% fossils == TRUE
    return(length(which(desc == TRUE)))
    }

  posteriors$mcc.data.node <- as.numeric(posteriors$mcc.data.node)
  posteriors$mcc.data.posterior <- as.numeric(posteriors$mcc.data.posterior)

  nodes <- posteriors$mcc.data.node
  fossil_descendents <- vector()

  for(i in 1:length(nodes)){
    fossil_descendents <- c(fossil_descendents, fossils_per_node(true_tree, nodes = nodes[i], fossils = fossils))
  }


  posteriors$fossil_descendents <- fossil_descendents

  #you should look at the posterior of the clades containig fossil species
  
  #parents <- mcc@phylo$edge[1:length(mcc@phylo$tip.label),1]
  #tips <- (mcc@phylo$tip.label)
  #tips  <- data.frame(tips, parents)
  #tips_word <- as.character(tips$tips)
  #tips <- tips[which(tips_word %in% fossils),]
  #tips$parents <- as.numeric(tips$parents)
  #fossil_descendent <- posteriors$mcc.data.index %in% tips$parents
  
  tree <- mcc@phylo
  n <- length(tree$tip.label)
  ee<-setNames(tree$edge[sapply(1:n,function(x,y)   
    which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ff <- which(names(ee) %in% fossils)
  
  gg <- ee[ff]
  names(gg) <- NULL
  
  posteriors$fossil_associated <-posteriors$mcc.data.node %in% gg
  
  
#this is all very annoying
#true_heights <- unique(nodeHeights(true_tree)[,1])
#true_ages <- tree.length(true_tree) - true_heights
#true_edges <- unique(true_tree$edge[,1])
#true_ages <- cbind(true_edges, true_ages)

  posteriors$treatments <- treatment_name
  #posteriors$run_names <- run_names[i]
  
  all_post <- rbind(all_post, posteriors)
  
  mean_post <- c(mean_post, mean(posteriors$mcc.data.posterior))

  write.csv(posteriors, file = "posteriors.csv")
  print(getwd())
  setwd("..")
}

mean_post <- cbind(run_names, mean_post, treatments)



all_post$treatments <- as.factor(all_post$treatments)

save(all_post, file = "posterior.Rdata")


setwd("~/Desktop/comb_trees_100")
load(file = "posterior.Rdata")

### ANLYSES ###
  
posterior_support_im <- ggplot(data = all_post, aes(x = mcc.data.posterior, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  scale_x_continuous(name = "Posterior node support", limits = c(0,1))
ggsave("posterior_support.png", plot = posterior_support_im, device = "png", width = 9.5, height = 6)

low <- subset(all_post, all_post$treatments == "Low")
median <- subset(all_post, all_post$treatments == "Median")
high <- subset(all_post, all_post$treatments == "High")

CI95(low$mcc.data.posterior)
CI95(median$mcc.data.posterior)
CI95(high$mcc.data.posterior)


#all_post$treatments <- as.factor(all_post$treatments)

boxplot(mcc.data.posterior ~ treatments, data = all_post)

post.aov <- aov(mcc.data.posterior ~ treatments + fossil_descendents, data = all_post)

post.aov <- aov(mcc.data.posterior ~ treatments + fossil_descendents + age, data = all_post)

post.fossil.aov <- aov(mcc.data.posterior ~ treatments + fossil_associated + age, 
                       data = all_post)


summary(post.fossil.aov)

TukeyHSD(post.aov)

post.kru <- kruskal.test(mcc.data.posterior ~ treatments, data = all_post)


beanplot(mcc.data.posterior ~ treatments, data = all_post,
         ylab = "posterior support", what = c(1,1,1,0),
         col = "darkslategray3")


beanplot(mcc.data.posterior ~ treatments, data = all_post,
         ylab = "posterior support", what = c(1,1,1,0),
         col = "darkslategray3")

beanplot(mcc.data.posterior ~ fossil_associated + treatments, data = all_post,
         ylab = "posterior support", what = c(1,1,1,0),
         col = "darkslategray3")

write.csv(all_post, file = "this.csv")

#mcc <- read_annotated(file = "sim.mcc.tre")

#we weant to get ages with 95CI

#we also want to get posterior for every clade
#then we can ask is there a difference between the average posteriors:
# - between treatments
# - across time
# - between fossils and extant taxa - is their an effect of having fossil taxa

#####################################
# Clade Level Extinction Experiment #
#####################################

setwd("~/Desktop/clade_extinction")

run_names <- seq(1,44,1)
run_names <- paste0("clade_extinction_", run_names)

all_post <- vector()
mean_post <- vector()

sym.dif <- vector()
bra.dif <- vector()
clade_size <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  this_run <- run_names[i]

  mcc <- read.mrbayes(file = "output/sim.mcc.tre")
  
  posteriors <- data.frame(mcc@data$node, mcc@data$index, mcc@data$posterior)
  confidence <- mcc@data$age_0.95_HPD
  high <- vector()
  low <- vector()
  age <- vector()
  for(i in 1:length(confidence)){
    x <- confidence[[i]]
    low <- c(low, x[1])
    high <- c(high, x[2])
    age <- c(age, mean(c(x[1], x[2])))
  }
  
  posteriors$age <- age
  posteriors$low95 <- low
  posteriors$high95 <- high
  
  nodes <- mcc@phylo$edge[,1]
  nodes <- unique(nodes)
  
  posteriors <- posteriors[which(posteriors$mcc.data.node %in% nodes),]
  
  true_tree <- read.tree("output/true.tree")
  fossils <- getExtinct(true_tree)
  
  fossils_per_node <- function(tree, nodes, fossils){
    desc <- Descendants(tree, nodes, type = "tips")
    desc <- desc[[1]]
    desc <- true_tree$tip.label[desc]
    desc <- desc %in% fossils == TRUE
    return(length(which(desc == TRUE)))
  }
  
  posteriors$mcc.data.node <- as.numeric(posteriors$mcc.data.node)
  posteriors$mcc.data.posterior <- as.numeric(posteriors$mcc.data.posterior)
  
  nodes <- posteriors$mcc.data.node
  fossil_descendents <- vector()
  
  for(i in 1:length(nodes)){
    fossil_descendents <- c(fossil_descendents, fossils_per_node(true_tree, nodes = nodes[i], fossils = fossils))
  }
  
  
  posteriors$fossil_descendents <- fossil_descendents
  
  #you should look at the posterior of the clades containig fossil species
  
  #parents <- mcc@phylo$edge[1:length(mcc@phylo$tip.label),1]
  #tips <- (mcc@phylo$tip.label)
  #tips  <- data.frame(tips, parents)
  #tips_word <- as.character(tips$tips)
  #tips <- tips[which(tips_word %in% fossils),]
  #tips$parents <- as.numeric(tips$parents)
  #fossil_descendent <- posteriors$mcc.data.index %in% tips$parents
  
  tree <- mcc@phylo
  n <- length(tree$tip.label)
  ee<-setNames(tree$edge[sapply(1:n,function(x,y)   
    which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ff <- which(names(ee) %in% fossils)
  
  gg <- ee[ff]
  names(gg) <- NULL
  
  posteriors$fossil_associated <-posteriors$mcc.data.node %in% gg
  
  
  #this is all very annoying
  #true_heights <- unique(nodeHeights(true_tree)[,1])
  #true_ages <- tree.length(true_tree) - true_heights
  #true_edges <- unique(true_tree$edge[,1])
  #true_ages <- cbind(true_edges, true_ages)
  
  posteriors$clade_size <- length(fossils)
  #this_run <- rep(run_names[i], length(posteriors$mcc.data.node))
  posteriors$run_names <- this_run
  
  all_post <- rbind(all_post, posteriors)
  
  mean_post <- c(mean_post, mean(posteriors$mcc.data.posterior))
  
  write.csv(posteriors, file = "posteriors.csv")
  
  #lets do something with tree distance as well
  
  sym.dif <- c(sym.dif, treedist(mcc@phylo, true_tree)[[1]])
  bra.dif <- c(bra.dif, treedist(mcc@phylo, true_tree)[[2]])
  clade_size <- c(clade_size, length(fossils))
  
  print(getwd())
  setwd("..")
}
post.aov <- aov(mcc.data.posterior ~ clade_size + run_names, data = all_post)
summary(post.aov)



