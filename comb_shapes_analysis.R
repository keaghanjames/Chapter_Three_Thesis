
require(CollessLike)
require(phytools)
require(parallel)
source("~/Desktop/cloud_run/functions.R") #or wherever


cores <- detectCores()
tt = log(50)-log(2)



##########################
# Calculatre Tree Length #
##########################

#calculate the distance of a tip to the current root

tree.length <- function(x) {
  max(diag(vcv(x)))
}

##################
# Root and Ultra #
##################

#a messy bit of code based on force.ultrametric() which makes all extant tips equal to the maximum age
#and to add a root to the tree
root.and.ultra <- function(tree, root_length = 1, scale = tt) {

  extant <- getExtant(tree) #extract extant tips
  h<-diag(vcv(tree))
  d<-max(h)-h
  d <- d[extant]
  ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
           y=tree$edge[,2])
  names(ii) <- tree$tip.label
  ii <- ii[extant] 
  names(ii) <- NULL
  tree$edge.length[ii]<-tree$edge.length[ii]+d
  tree.length <- tree.length(tree)#calculate the distance of a tip to the current root
  #add a new root edge
  tree$root.edge <- 0
  # give it length equal to the tree length and add 1 length between it and the tree root
  root<-list(edge=matrix(c(3,1,3,2),2,2,byrow=TRUE),
             edge.length=c((tree.length+root_length),1),tip.label=c("Z","NA"),
             Nnode=1)
  class(root)<-"phylo"
  tree<-paste.tree(root,tree)
  #plotTree(tree, type = "fan")
  tree$edge.length<- tree$edge.length/max(nodeHeights(tree)[,2])*scale
  tree <- tree
}

#generate some birth-death trees
trees <- pbtree(n = 30, b = 1.7, d = 0.7, t = tt, scale = log(50)-log(2), nsim = 1000)
#to calculate colless_scores for the full tree
trees_igraph <- lapply(trees, FUN = as.igraph.phylo)
colless_scores <- mclapply(trees_igraph, colless.like.index, mc.cores = cores)
colless_scores <- unlist(colless_scores)

#the problem with this approach is that there is a clear correlation between 
#coless scores and the number of extinct lineages on a tree
extinct <- vector()
for(i in 1:length(trees)){
  x <- length(getExtinct(trees[[i]]))
  extinct <- c(extinct, x)
}

plot(extinct, colless_scores)
cor.test(extinct, colless_scores)

#this means we have a bias in our data, with the different
#balanace treatments have biased numbers of fossil species
#is this simply because more imbalanced trees have more tips 

size <- vector()
for(i in 1:length(trees)){
  x <- length(trees[[i]]$tip.label)
  size <- c(size, x)
}

plot(size, colless_scores)
cor.test(size, colless_scores)

#yes that seems to be driving the correlation - let's see what happens if
#we standerdise them
colless_scores_standerdised <- colless_scores/size

#now let's look if that clears up the correlation between extinction and
#imbalance

plot(extinct, colless_scores_standerdised)
cor.test(extinct, colless_scores_standerdised)

#it certainly improves it - we loose about 40 per cent of the correlation
#however, there is still a corelation - thus having more extinct tips drives
#up imbalance


#therefor may be better off taking the colless score after
#forcing the tree to be ultramtric
#lineages

#to calculate colless scores for the ultrametic tree
ultrametric_trees <- mclapply(trees, force.ultrametric)
ultrametric_trees_igraph <- mclapply(ultrametric_trees, FUN = as.igraph.phylo)


colless_scores <- mclapply(ultrametric_trees_igraph, colless.like.index, mc.cores = cores)
colless_scores <- unlist(colless_scores)
colless_scores_standerdised <- colless_scores/size
names(colless_scores_standerdised) <-  seq(1,1000,1)
#colless_scores_standerdised <- sort(colless_scores, decreasing = T)

#does this account for the bias
plot(extinct, colless_scores_standerdised)
cor.test(extinct, colless_scores_standerdised)
#no - there is still a correlation, although it is worth noting that the strength of that correlation is very low

#let's try an alternative approach where we measure the imablance of just the extant tips on our tree
#function to trim extint tips from tree
trim_dead_ends <- function(x){
  drop.tip(x, getExtinct(x))
}

extant_trees <- mclapply(trees, trim_dead_ends)
extant_trees_igraph <- mclapply(extant_trees, FUN = as.igraph.phylo)


colless_scores <- mclapply(extant_trees_igraph, colless.like.index, mc.cores = cores)
colless_scores <- unlist(colless_scores)
colless_scores_standerdised <- colless_scores
names(colless_scores_standerdised) <-  seq(1,1000,1)
#does this account for the bias
plot(extinct, colless_scores_standerdised)
cor.test(extinct, colless_scores_standerdised)
#yes it does!
colless_scores <- sort(colless_scores, decreasing = T)



#these latter three approaches represent different contributions of extinction and speciation
#to imbalance of the trees - in the first extinction creates imbalanced trees, in the second speciation 
#and extinction are contributing to imbalance - meaning the diversification as a whole is creating imbalance
#and in the last is only speciation that is creating an imabalnced tree...

#force ultrametric and give all the trees root tip Z
trees <- mclapply(trees, root.and.ultra)
trees <- mclapply(trees, ladderize)
#and this is a loop to ensure that the multiphylo object is saved properly
aa <- trees[[1]]
for(i in 2:length(trees)){
  aa <- c(aa, trees[[i]])
}
trees <- aa

#we still extract from trees - so we are still getting the trees with extinct lineages
high_colless <- colless_scores[1:20]
high_colless_trees <- trees[as.numeric(names(high_colless))]
median_colless <- colless_scores[491:510]
median_colless_trees <- trees[as.numeric(names(median_colless))]
low_colless <- colless_scores[981:1000]
low_colless_trees <- trees[as.numeric(names(low_colless))]

#discard the trees we didn't sample
trees <- c(high_colless_trees, median_colless_trees, low_colless_trees)
#and create a little vector to remember who is who
trees_colless <- c(high_colless, median_colless, low_colless)
names(trees_colless[1:20]) <- "high"
names(trees_colless[21:40]) <- "median"
names(trees_colless[41:60]) <- "low"

#generate trait, alignment and tipdate data
traits <- mclapply(trees, traits.simulate)
alignments <- mclapply(trees, alignment.simulate)
tip_dates <- mclapply(trees, tip.dates)

#we might also want to look at the number of extinct species on our trees to ensure that 
#we don't have any that are far too small
extinct <- vector()
for(i in 1:length(trees)){
  x <- length(getExtinct(trees[[i]]))
  extinct <- c(extinct, x)
}
hist(extinct, main = "Fossil sp. similated trees")
quantile(extinct)

#also we might want to know something about the number of tips
tips_n <- vector()
for(i in 1:length(trees)){
  tips_n <- c(tips_n, length(trees[[i]]$tip.label))
}

#and why not just get the number of extant too
extant_num <- vector()
for(i in 1:length(trees)){
  x <- length(getExtant(trees[[i]]))
  extant_num <- c(extant_num, x)
}


#now we want to start putting together the files to run our experiment
setwd("~/Desktop/comb_trees")
run_names_high <- seq(1, 20, 1)
run_names_high <- paste0("High_", run_names_high)
run_names_median <- seq(1, 20, 1)
run_names_median <- paste0("Median_", run_names_median)
run_names_low <- seq(1, 20, 1)
run_names_low <- paste0("Low_", run_names_low)
run_names <- c(run_names_high, run_names_median, run_names_low)



for(i in 1:length(run_names)){
  system(paste("mkdir", run_names[i]))
  setwd(paste(run_names[i]))
  system(paste("mkdir data"))
  system(paste("mkdir output"))
  write.tree(trees[[i]], file = "output/true.tree")
  write.nexus.data(traits[[i]], file = "data/states.nex", format = "standard",interleaved = F)
  write.nexus.data(alignments[[i]], file = "data/alignment.nex")
  write.table(tip_dates[[i]], file = "data/taxa.tsv", quote = F, row.names = F, sep = "\t")
  setwd("..")
}

run_list <- file("run_list_c6.txt")
file.names <- vector()
for(i in 51:60){
  file.name <- paste0("cd ~/comb_trees/", run_names[i])
  file.name <- paste0(c(file.name, "rb ~/comb_trees/scripts/initialise.Rev"))
  file.names <- c(file.names, file.name)
}
writeLines(file.names, run_list)
close(run_list)

############################################

#okay - now we need something that is going to read in our old trees, randomly sample the fossils,
#alter the traits and tip dates accordingly

setwd("~/Desktop/comb_trees_50")
sample.prop <- 0.5

for(i in 1:length(run_names)){
setwd(run_names[i])

#sample our fossils and create a new tree with just the sampled
tree <- read.tree("output/true.tree")
extinct <- getExtinct(tree)
sample_fossils <- sample(extinct, round(length(extinct)*sample.prop), replace = F)
sample.tree <- drop.tip(tree, tree$tip.label[-match(c(getExtant(tree), sample_fossils), tree$tip.label)])
write.tree(sample.tree, file = "output/true.tree")

traits <- read.nexus.data("data/states.nex")
traits <- traits[c(getExtant(sample.tree), sample_fossils)]
write.nexus.data(traits, file = "data/states.nex")

taxa <- read.table("data/taxa.tsv")
rownames(taxa) <- taxa$V1
tips <- sample.tree$tip.label
taxa <- taxa[tips,]
rownames(taxa) <- NULL
write.table(taxa, file = "data/taxa.tsv", quote = F, row.names = F, col.names = c("taxon", "min", "max"), sep = "\t")

setwd("..")
}
