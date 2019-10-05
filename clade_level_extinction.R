require(phytools)
require(phangorn)
require(ape)
require(Rfast)
require(gtools)
require(caper)
require(knitr)
require(rwty)
require(parallel)
source("~/Desktop/cloud_run/functions.R") #or wherever

setwd("~/Desktop/clade_extinction")

###############################
# Clade-Level Extinction Tree #
###############################

#you could rewrite this so that it it was a transformation applied to extisting trees...

clade.extinction.tree <- function(tip_stop = 50, time_stop = log(50)-log(2), root_length = 1, 
                                 birth_rate = 1.2, extinction_rate = 0, scale = time_stop,
                                 min_fossils = tip_stop*0.1, max_fossils = tip_stop*0.5){
mono <- FALSE
while(mono == FALSE) {
tree <- pbtree(b = birth_rate, extinction_rate = 0, n = tip_stop, t = time_stop, scale = scale)
extant <- getExtant(tree)
bg_extinctions <- getExtinct(tree) 

target <- Descendants(tree, type = "tips")
target <- lapply(target, length)
target <- which(target >= 10 & target <= 20)
target <- sample(target, 1)

attach <- getParent(tree, target)

edge <- which(tree$edge[,1] == attach & tree$edge[,2] == target)
edge_length <- tree$edge.length[edge]

last <- length(tree$tip.label)+1

desc <- Descendants(tree, target, type = "tips")
desc <- desc[[1]]
desc <- tree$tip.label[desc]
is.monophyletic(tree, desc)

#try binding a tip to preserve the structure
#tree <- bind.tip(tree, "to_die", edge.length = tree.length(extinct_tree), where = attach)


extinct_tree <- drop.tip(tree,tree$tip.label[-match(desc, tree$tip.label)])
extinct_tree$tip.label <- paste0("c", extinct_tree$tip.label)

shorten <- function(x){ 
  x*sample(seq(0.1,0.8,.05),1)
  }

extinct_tree$edge.length<- lapply(extinct_tree$edge.length, shorten)
edge_length <- shorten(edge_length)

clade_extinct <- extinct_tree$tip.label

tree.length <- function(x) {
  max(diag(vcv(x)))
}

#tree <- bind.tree(tree, extinct_tree, where = attach)
tree <- bind.tree(tree, extinct_tree, where = attach)
tree <- multi2di(tree)
poly <- which(tree$edge.length == 0)
tree$edge.length <- unlist(tree$edge.length)

tree.length <- function(x) {
  max(diag(vcv(x)))
}

tree$edge.length[poly] <- (edge_length*0.5) #do this smarter
tree <- drop.tip(tree, desc)

#traits <- traits.simulate(bind.tree)
#fancyTree(tree, type = "extinction")
a <- extant %in% desc #which members of the murdered clade are in extant
b <- which(a == FALSE) #get positions of extant tips who are not murdered tips
extant <- extant[b] #select only those extant tips who are not murdered tips 

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
tree$edge.length<- tree$edge.length/max(nodeHeights(tree)[,2])*1
mono <- is.monophyletic(tree, extinct_tree$tip.label)
}
#plot(ladderize(tree))
tree <- tree
}

#and here is where we can actually start to build our sample
trees <- pbtree(n = 20)
for(i in 1:2000){
  x <- clade.extinction.tree()
  trees <- c(trees, x)
  print(i)
}
trees <- trees[2:length(trees)]
save(trees, file = "1000trees.Rdata")

#we probably also want to hold onto just the tree with the extinct clade

extinct_trees <- pbtree(n = 20)
for(i in 1:length(trees)){
  tree <- trees[[i]]
  extinct <- drop.tip(tree,tree$tip.label[-match(getExtinct(tree), tree$tip.label)])
  extinct_trees <- c(extinct_trees, extinct)
}
extinct_trees <- extinct_trees[2:length(extinct_trees)]
save(extinct_trees, file = "1000extincttrees.Rdata")

n.extinct <- vector()
for(i in 1:length(extinct_trees)){
  x <- extinct_trees[[i]]
  x <- length(x$tip.label)
  n.extinct <- c(n.extinct, x)
}

positions <- vector()    
samples <- sort(unique(n.extinct))
samples #actually check what these are, weird 1s turn up all the time
length(samples)
samples <- samples[2:12]

for(i in 1:length(samples)){
  x <- which(n.extinct == samples[i])
  x <- sample(x,4, replace = F)
  positions <- c(positions, x)
}

trees <- trees[positions]
extinct_trees <- extinct_trees[positions]

save(trees, file = "trees.Rdata")
save(extinct_trees, file = "extinct_trees.Rdata")

#generate trait, alignment and tipdate data
traits <- mclapply(trees, traits.simulate)
alignments <- mclapply(trees, alignment.simulate)
tip_dates <- mclapply(trees, tip.dates)


#creating our directories for the analysis

run_names <- seq(1, length(trees), 1)
run_names <- paste0("clade_extinction_", run_names)

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


run_list <- file("run_list_c8.txt")
file.names <- vector()
for(i in 40:44){
  file.name <- paste0("cd ~/clade_extinction/", run_names[i])
  file.name <- paste0(c(file.name, "rb ~/clade_extinction/scripts/initialise.Rev"))
  file.names <- c(file.names, file.name)
}
writeLines(file.names, run_list)
close(run_list)


