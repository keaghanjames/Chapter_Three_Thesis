require(phytools)
require(phangorn)
require(ape)
require(Rfast)
require(gtools)
require(caper)
require(knitr)
require(rwty)
require(parallel)


setwd("~/Desktop/cloud_run") #or wherever you want to put it
source("~/Desktop/cloud_run/functions.R") #or wherever
run_names <- seq(1, 20, 1)
run_names <- paste0("Run_", run_names)

run_list <- file("run_list.txt")

file.names <- vector()
for(i in (1:length(run_names))){
  file.name <- paste0("cd ~/cloud_run/", run_names[i])
  file.name <- paste0(c(file.name, "rb ~/cloud_run/scripts/initialise.Rev"))
  file.names <- c(file.names, file.name)
}
writeLines(file.names, run_list)
close(run_list)

numCores <- detectCores()
mclapply(run_names, write_files, mc.cores = numCores)
mclapply(run_names, simple.run, mc.cores = numCores)
#lapply(run_names, simple.run)
