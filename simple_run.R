require(phytools)
require(phangorn)
require(ape)
require(Rfast)
require(gtools)
require(caper)
require(knitr)
require(rwty)

#load in the necessary functions
source("~/Desktop/cloud_run/functions.R")
#setwd("~/Desktop/SimStudy/simple")
#code.path <- "~/Desktop/simple/essentials/codepath/"

#first we simulate a tree using the fbd.simulate function
all.tree <- fbd.simulate(tip_stop = 30, extinction_rate = 0.9)
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


#WILL NEED TO EDIT REVBAYES SOURCE FILE

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


