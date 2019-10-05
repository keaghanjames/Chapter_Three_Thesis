require(phytools)
require(phangorn)
require(ape)
require(Rfast)
require(gtools)
require(caper)
require(knitr)
require(rwty)

##########################
# Calculatre Tree Length #
##########################

#calculate the distance of a tip to the current root

tree.length <- function(x) {
  max(diag(vcv(x)))
  }

##########################################
# Fossilised Birth-Death Tree Simulation #
##########################################

fbd.simulate <- function(tip_stop = 30, time_stop = log(50)-log(2), root_length = 1, birth_rate = 1.2, extinction_rate = 0.9, scale = time_stop, nsim = 1) {
  #simulate a tree
  num.fossils <- vector()
  #for(i in 1:100){
  tree<-pbtree(b = birth_rate, d = extinction_rate, t = time_stop, n = tip_stop, scale = scale, nsim = nsim)
  extant <- getExtant(tree) #extract extant tips
  
  #a messy bit of code based on force.ultrametric() which makes all extant tips equal to the maximum age
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
  
#####################
#Simulate Trait Data #
#####################


traits.simulate <- function(tree, tips = tree$tip.label, n_traits = 100, n_states = 2, ...){
  states <- vector()
    for(i in 1:n_traits){
      nchanges <- round(rlnorm(1, meanlog = 1.7, sdlog = .5)) #built stochasticity into the number of changes
      q <- nchanges/sum(tree$edge.length)
      Q <- matrix(q,n_states,n_states)
      diag(Q) <- -q
      rownames(Q)<-colnames(Q)<-as.character(c(1:n_states)-1)
      trait.tree <- sim.history(tree,Q=Q, ...)
    states <- cbind(states, trait.tree$states)
    print(nchanges)
  }
  states <- states
}

#write.nexus.data(states, file = "all/data/states.nex", format = "standard",interleaved = F)


###########################################
# Simulate Alignment Data for Extant Taxa #
###########################################

alignment.simulate <- function(tree, seq_length = 1000, rate = .015, extant.only = T, format = "matrix", ...){
  alignment <- genSeq(tree, l = seq_length, format = format, rate = rate, ...) #simulate sequence alignment
  if(extant.only == T){
  extant <- getExtant(tree)
  extant.alignment <- alignment[extant,] #extract sequence data just for extant tips
  extant.alignment <- as.DNAbin(extant.alignment) #convert to DNAbin format
  } else {
  alignment <- as.DNAbin(alignment)  
  }
}

#write.nexus.data(extant.alignment, file = "all/data/alignment.nex", format = "dna", interleaved = F)



#n = the number of fossils to be sampled

###########################
# Random Sampling Functon #
###########################


random.fossils <- function(tree, n) {
  extinct <- getExtinct(tree)
  extant.tree <- drop.tip(tree, extinct)
  sample <- sample(extinct, n)
  sample <- c(sample, extant.tree$tip.label)
  rnd.fossils.tree <- drop.tip(tree, tree$tip.label[-match(sample, tree$tip.label)])
  rnd.fossils.tree
}

#write.tree(rnd.fossils.tree, file = "random.tree")
#sample.states <- states[tips,]



##############################
# Brownian Sampling Function #
##############################

brownian.fossils <- function(tree, n, ...) {
  extinct <- getExtinct(tree)
  extant.tree <- drop.tip(tree, extinct)
  BMtrait <- fastBM(tree) #returns a matrix of tip lables and bm value
  BMtrait <- BMtrait[extinct]
  BMtrait <- sort(BMtrait, decreasing = T) #organises tips from largest to smallest
  sample <- names(BMtrait[1:n]) #samples the half with the highest BM
  #add the extant taxa to the sample and create to sample tree
  sample <- c(sample, extant.tree$tip.label)
  brw.fossils.tree <- drop.tip(tree,tree$tip.label[-match(sample, tree$tip.label)])
  brw.fossils.tree
}

#write.tree(brw.fossils.tree, file = "brownian.tree")
#sample.states <- states[tips,]


########################
# PD Sampling Function #
########################

#requires the program pda to be installed in working directory
#will also write an ultrametric extinct tree and a txt file output of pda

pd.fossils <- function (tree, n, force_ultrametric = F) {
  #set up file for runing pda
  extant <- getExtant(tree)
    if(force_ultrametric == T) {
      extinct.tree <- drop.tip(tree, extant)
      extinct.tree <- force.ultrametric(extinct.tree)
    } else {
      extinct.tree <- drop.tip(tree, extant)
    }
    write.tree(extinct.tree, file = "extinct.tree")
    execute <- combine_words(c("./pda extinct.tree -k", n, "-g pda.extinct.txt"), sep = " ")
    system(execute)
    
  #read in output of pda and convert to list of selected tips
  sample <- read.delim(file = "pda.extinct.txt")
  sample <- sample[[1]]
  sample
  sample <- as.character(sample[13:(n+12)])
  
  #create tree with sampled fossil tips
  extinct <- getExtinct(tree)
  extant.tree <- drop.tip(tree, extinct)
  sample <- c(sample, extant.tree$tip.label)
  pd.fossils.tree <- drop.tip(tree,tree$tip.label[-match(sample, tree$tip.label)])
  pd.fossils.tree
}

#write.tree(pd.fossils.tree, file = "pd.tree")
#sample.states <- states[tips,]


#####################
# Extract tip dates #
#####################

#returs a data frame of first and last appearance dates for every tip on the tree

tip.dates <- function(tree, age.uncertainty = T , age.uncertainty.ammount = sample(seq(0.1,1,.1),1)){
  
  tips <- tree$tip.label
  extinct <- getExtinct(tree)
  LADs <- vector()
  FADs <- vector()
  nodes<-sapply(tips,function(x,y) which(y==x),y=tree$tip.label)
  edge.lengths<-setNames(tree$edge.length[sapply(nodes,
        function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))
  
  if(age.uncertainty == T){

    for(i in 1:length(tips)){
      lineage.age <- pd.calc(tree, tips[i], root.edge = T)
      lineage.age <- lineage.age[1]
      LADs <- c(LADs, lineage.age)
      first.appearance <- lineage.age - (age.uncertainty.ammount*edge.lengths[i])
      FADs <- c(FADs, first.appearance)
    }
    
  } else {
    for(i in 1:length(tips)){
      lineage.age <- pd.calc(tree, tips[i], root.edge = T)
      lineage.age <- lineage.age[1]
      LADs <- c(LADs, lineage.age)
      FADs <- c(FADs, lineage.age)
    }
  }

  tree.length <- tree.length(tree)#calculate the distance of a tip to the current root
  #and we use tree.length to get tip dates for the extinct taxa
  LADs <- tree.length - LADs
  FADs <- tree.length - FADs

  #this messy bit of codes fixes the ages of all of the extant tips to zero
  names(LADs) <- tips
  last.extcintion <- min(LADs[extinct])
  LADs[LADs < last.extcintion] <- 0
  names(LADs) <- NULL

  #create the data frame of our sample
  sample.taxa <- data.frame(taxon =  tips, min = LADs, max = FADs)

}

###############################
# Get traits/dates for sampled tips #
###############################

#where 'data' is a the data frame produced by traits.simulate() or tip.dates()
#and 'tree' is the tree of the subset of tips to be sampled

sample.data <- function(data,sample_tree){
  data[sample_tree$tip.label,]
}


############################################
# Create Set of Directories for Simulation #
############################################

#run_names <- read.table("run_names.csv", head = T, sep = ",", stringsAsFactors = F)
#run_names <- run_names$runName

write_files <- function(x) {
  file.name <- (x)
  system(paste("mkdir", file.name))
  setwd(file.name)
  
  system("mkdir all")
  system("mkdir all/data")
  system("mkdir all/output")
  
  system("mkdir pd")
  system("mkdir pd/data")
  system("mkdir pd/output")
  system("cp ~/Desktop/cloud_run/essentials/pda pda")
  
  system("mkdir rnd")
  system("mkdir rnd/data")
  system("mkdir rnd/output")
  
  
  system("mkdir brw")
  system("mkdir brw/data")
  system("mkdir brw/output")
  
  setwd("..")
}

#########################################
# Execute Simulation Across Directories #
#########################################

simple.run <- function(x, ...) {
  
  #load in the necessary functions
  source("~/Desktop/cloud_run/functions.R")
  #setwd("~/Desktop/SimStudy/simple")
  setwd(x)  
  #code.path <- "~/Desktop/simple/essentials/codepath/"
  
  #first we simulate a tree using the fbd.simulate function
  all.tree <- fbd.simulate(tip_stop = 30, extinction_rate = 0.9, ...)
  write.tree(all.tree, file = "all/output/all.tree")
  
  #then we generate a set of traits using traits.simulate()
  traits <- traits.simulate(all.tree, n_traits = 100)
  write.nexus.data(traits, file = "all/data/states.nex", format = "standard",interleaved = F)
  
  
  #next we generate an alignment for the tree 
  #but extract only the extant taxa's sequences
  alignment <- alignment.simulate(all.tree, extant.only = T)
  #and we save the alignment to working directory
  write.nexus.data(alignment, file = "all/data/alignment.nex", format = "dna", interleaved = F)
  write.nexus.data(alignment, file = "rnd/data/alignment.nex", format = "dna", interleaved = F)
  write.nexus.data(alignment, file = "brw/data/alignment.nex", format = "dna", interleaved = F)
  write.nexus.data(alignment, file = "pd/data/alignment.nex", format = "dna", interleaved = F)
  
  #and then we get the tip dates for our all tree
  all.taxa <- tip.dates(all.tree)
  write.table(all.taxa, file = "all/data/taxa.tsv", quote = F, row.names = F, sep = "\t")
  
  
  #next we get a sample of extinct taxa using random.fossils()
  rnd.tree <- random.fossils(all.tree, n = 15)
  write.tree(rnd.tree, file = "rnd/output/rnd.tree")
  rnd.traits <- sample.data(traits, rnd.tree)
  write.nexus.data(rnd.traits, file = "rnd/data/states.nex", format = "standard",interleaved = F)
  rnd.taxa <- sample.data(all.taxa, rnd.tree)
  write.table(rnd.taxa, file = "rnd/data/taxa.tsv", quote = F, row.names = F, sep = "\t")
  
  
  #next we get a sample of extinct taxa using brownian.fossils()
  brw.tree <- brownian.fossils(all.tree, n = 15)
  write.tree(brw.tree, file = "brw/output/brw.tree")
  brw.traits <- sample.data(traits, brw.tree)
  write.nexus.data(brw.traits, file = "brw/data/states.nex", format = "standard",interleaved = F)
  brw.taxa <- sample.data(all.taxa, brw.tree)
  write.table(brw.taxa, file = "brw/data/taxa.tsv", quote = F, row.names = F, sep = "\t")
  
  
  #next we get a sample of extinct taxa using pd.fossils()
  #remember this requires that pda be in the working directory
  
  pd.tree <- pd.fossils(all.tree, n = 15, force_ultrametric = F)
  write.tree(pd.tree, file = "pd/output/pd.tree")
  pd.traits <- sample.data(traits, pd.tree)
  write.nexus.data(pd.traits, file = "pd/data/states.nex", format = "standard",interleaved = F)
  pd.taxa <- sample.data(all.taxa, pd.tree)
  write.table(pd.taxa, file = "pd/data/taxa.tsv", quote = F, row.names = F, sep = "\t")
  
  setwd("..")
  
}