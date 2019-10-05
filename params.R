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

cynthwave <- c("#FFA9E7", "#FF84E8", "#7F2CCB", "#414361", "#2A2D43")
windwaker <- c("#19647E", "#28AFB0", "#DDCECD", "#37392E", "#CF5C36")

CI95 <- function(x){
  a <- mean(x)
  n <- length(x)
  s <- sd(x)
  error <- qnorm(0.975)*s/sqrt(n)
  left <- a-error
  right <- a+error
  c(a, left,right)
}


##############################
# Fossil Sampling Experiment #
##############################

setwd("~/Desktop/cloud_run")

run_names <- seq(1, 50, 1)
run_names <- paste0("Run_", run_names)
speciation <- vector()
extinction <- vector()
#treatments <- vector()
rec.speciation <- vector()
rec.extinction <- vector()

for(i in 1:length(run_names)){
  setwd(run_names[i])
  treatments <- c("all", "rnd", "brw", "pd")
  for(j in 1:length(treatments)){
    setwd(treatments[j])
      log <- read.table("output/sim.log", header = T)
      log <- (log[2501:10001,])
      speciation <- c(speciation, median(log$speciation_rate))
      extinction <- c(extinction, median(log$extinction_rate))
      
      # does the paramter fall in the C1 #
      
      sp.up <- CI95(log$speciation_rate)[2]
      sp.low <-  CI95(log$speciation_rate)[3]
      
      rec.speciation <- c(rec.speciation, (1.7 <= sp.up & 1.7 >= sp.low))
      
      ex.up <- CI95(log$extinction_rate)[2]
      ex.low <- CI95(log$extinction_rate)[3]  
      rec.extinction <- c(rec.extinction, (ex.up >= 0.9 & ex.low <= 0.9))
      
      print(getwd())
      setwd("..")
  }
  setwd("..")
}

treatments <- rep(c("Cmp","Rnd", "Trt", "Div"), 50)

#speciation <- speciation - 1.7 #the true value?
#extinction <- extinction - 0.7 #the true value?


parameters <- tibble(speciation, extinction, treatments)
parameters <- parameters[-which(is.na(speciation)),]
parameters$treatments <- as.factor(parameters$treatments)

save(parameters, file = "parameters.Rdata")
setwd("~/Desktop/cloud_run")
load(file = "parameters.Rdata")

shapiro.test(parameters$speciation[which(parameters$treatments == "Cmp")])
shapiro.test(parameters$speciation[which(parameters$treatments == "Rnd")])
shapiro.test(parameters$speciation[which(parameters$treatments == "Trt")])
shapiro.test(parameters$speciation[which(parameters$treatments == "Div")])

cmp <- subset(parameters, parameters$treatments == "Cmp")
CI95(cmp$speciation)
CI95(cmp$extinction)

rnd <- subset(parameters, parameters$treatments == "Rnd")
CI95(rnd$speciation)
CI95(rnd$extinction)

trt <- subset(parameters, parameters$treatments == "Trt")
CI95(trt$speciation)
CI95(trt$extinction)

div <- subset(parameters, parameters$treatments == "Div")
CI95(div$speciation)
CI95(div$extinction)

aov.speciation <- aov(speciation ~ treatments, data = parameters)
summary(aov.speciation)


kw.speciation <- kruskal.test(speciation ~ treatments, data = parameters)
kw.speciation

require(FSA)
Dunn <- dunnTest(speciation ~ treatments, data = parameters, method = "bonferroni")
Dunn

kw.extinction <- kruskal.test(extinction ~ treatments, data = parameters)
kw.extinction
Dunn <- dunnTest(extinction ~ treatments, data = parameters, method = "bonferroni")
Dunn


Dunn <- dunnTest(speciation ~ treatments, data = parameters, method = "bh")

cldList(comparison = Dunn$Comparison, p.value = Dunn$P.adj, threshold  = 0.05)

ggplot(data = parameters, aes(x = speciation, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), (expand = c(0.01, 0))) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  geom_vline(xintercept = 1.7, color = windwaker[5], size = 1.3) +
  scale_x_continuous(name = "Speciation rate", limits = c(0.5, 2.5)) +
  annotate("text", x = 2, y = 6, label = expression(paste(lambda, "= 1.7")), color = windwaker[5], size = 6)+
  ggtitle("a")
#ggsave(filename = "speciation.png", plot = last_plot(), device = "png", width = 10, height = 5.5)
  
aov.extinction <- aov(extinction ~ treatments, data = parameters)
summary(aov.extinction)

kw.extinction <- kruskal.test(extinction ~ treatments, data = parameters)
kw.extinction

require(FSA)
Dunn <- dunnTest(extinction ~ treatments, data = parameters, method = "bh")
Dunn

require(rcompanion)
cldList(comparison = Dunn$Comparison, p.value = Dunn$P.adj, threshold  = 0.05)

Mann.Whitney <- pairwise.wilcox.test(parameters$extinction, parameters$treatments,
                                     p.adjust.method = "bonferroni")
Mann.Whitney
Mann.Whitney <- Mann.Whitney$p.value
fullPTable(Mann.Whitney)

require(multcompView)

multcompLetters(Mann.Whitney, compare="<", threshold=0.05, Letters=letters, reversed = FALSE)

extinction.plot <- ggplot(data = parameters, aes(x = extinction, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Cmp", "Rnd", "Trt", "Div"), (expand = c(0.01, 0))) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  geom_vline(xintercept = 0.9, color = windwaker[5], size = 1.3) +
  annotate("text", x = 1.2, y = 6, label = expression(paste(mu, "= 0.9")), color = windwaker[5], size = 6) +
  scale_x_continuous(name = "Extinction rate", limits = c(0.5, 2.5)) +
  ggtitle("b")

parameters.plot <- grid.arrange(speciation.plot, extinction.plot, nrow = 2)
ggsave(filename = "parameters.png", plot = parameters.plot, device = "png", width = 9, height = 8)


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

speciation <- vector()
extinction <- vector()
#treatments <- vector()
rec.speciation <- vector()
rec.extinction <- vector()


for(i in 1:length(run_names)){
  setwd(run_names[i])
  log <- read.table("output/sim.log", header = T)
  log <- (log[2501:10001,])
  speciation <- c(speciation, median(log$speciation_rate))
  extinction <- c(extinction, median(log$extinction_rate))
  
  # does the paramter fall in the C1 #
  
  sp.up <- CI95(log$speciation_rate)[2]
  sp.low <-  CI95(log$speciation_rate)[3]
  
  rec.speciation <- c(rec.speciation, (1.7 <= sp.up & 1.7 >= sp.low))
  
  ex.up <- CI95(log$extinction_rate)[2]
  ex.low <- CI95(log$extinction_rate)[3]  
  rec.extinction <- c(rec.extinction, (ex.up >= 0.9 & ex.low <= 0.9))
  
  print(getwd())
  setwd("..")
}

#treatments <- c(rep("Low", 150020), rep("Median", 150020), rep("High", 150020))

parameters <- tibble(speciation, extinction, treatments)
parameters$treatments <- as.factor(parameters$treatments)

save(parameters, file = "parameters.Rdata")
setwd("~/Desktop/comb_trees_100")
load(file = "parameters.Rdata")

kw.speciation <- kruskal.test(speciation ~ treatments, data = parameters)
kw.speciation

require(FSA)
Dunn <- dunnTest(speciation ~ treatments, data = parameters, method = "bonferroni")
Dunn

kw.extinction <- kruskal.test(extinction ~ treatments, data = parameters)
kw.extinction
Dunn <- dunnTest(extinction ~ treatments, data = parameters, method = "bonferroni")
Dunn


low <- subset(parameters, parameters$treatments == "Low")
CI95(low$speciation)
CI95(low$extinction)

median <- subset(parameters, parameters$treatments == "Median")
CI95(median$speciation)
CI95(median$extinction)

high <- subset(parameters, parameters$treatments == "High")
CI95(high$speciation)
CI95(high$extinction)

aov.speciation <- aov(speciation ~ treatments, data = parameters)
summary(aov.speciation)

speciation.plot <- ggplot(data = parameters, aes(x = speciation, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  geom_vline(xintercept = 1.7, color = windwaker[5], size = 1.3) +
  annotate("text", x = 1.3, y = 5, label = expression(paste(lambda, " = 1.7")), color = windwaker[5], size = 6)+
  scale_x_continuous(name = "Speciation rate", limits = c(0.05, 2.5))+
  ggtitle("c")



aov.extinction <- aov(extinction ~ treatments, data = parameters)
summary(aov.extinction)
TukeyHSD(aov.extinction)

extinction.plot <- ggplot(data = parameters, aes(x = extinction, y = treatments, fill = treatments)) + 
  geom_density_ridges(scale = 4, alpha = .7) +
  scale_y_discrete(limits = c("Low", "Median", "High"), expand = c(0.01, 0)) +
  scale_fill_manual(values = windwaker, guide = FALSE)+
  theme_ridges(center_axis_labels = T)+
  theme(axis.title.y = element_blank())+
  geom_vline(xintercept = 0.9, color = windwaker[5], size = 1.3) +
  annotate("text", x = 0.6, y = 5, label = expression(paste(mu, " = 0.9")), color = windwaker[5], size = 6)+
  scale_x_continuous(name = "Extinction rate", limits = c(0.05, 2.5))+
  ggtitle("d")

parameters.plot <- grid.arrange(speciation.plot, extinction.plot, nrow = 2)
ggsave(filename = "parameters.png", plot = parameters.plot, device = "png", width = 9, height = 8)


kw.extinction <- kruskal.test(extinction ~ treatments, data = parameters)


