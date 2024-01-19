######Code for Cardoso et al. (2024)
######Calculating functional diversity metrics using neighbor-joining trees
######code authors: Cardoso, Guillerme, Mammola, Matthews, Carvalho
######Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

library(ape)
library(BAT)

######analyses for Fig 1 and demonstration of how the method works

mat = matrix(c(0, 1, 7, 7, 1, 0, 6, 6, 7, 6, 0, 4, 7, 6, 4, 0), nrow = 4)
colnames(mat) = rownames(mat) = c("A", "B", "C", "D")
mat = as.dist(mat)
tree = nj(mat)
plot(tree)
plot(tree, "u")

comm1 = c(1,1,1,0)
comm2 = c(1,1,0,1)
comm = rbind(comm1, comm2)
colnames(comm) = c("A", "B", "C", "D")

alpha(comm, tree)
originality(comm, tree, relative = FALSE)
uniqueness(comm, tree, relative = FALSE)
contribution(comm, tree, F, relative = FALSE)
dispersion(comm, tree, relative = FALSE)
evenness(comm, tree)
beta(comm, tree)

# R code to generate the simulation and explore correlations

## ------------
## Preparatory steps
## ------------

## Set seed

set.seed(42)

## Set working directory

setwd("/Users/stefanomammola/Desktop/PAPERS IN CORSO/Cardoso et al. HYPERTREE")

## Loading R packages for simulations

library("dads") #devtools::install_github("TGuillerme/dads", force = TRUE)
library("BAT")

# Custom function
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

## ------------
## Running multiple simulations
## ------------

# Set scenarios
n.trait <- c(1,2,4,8)
n.sp    <- 100
n.sim   <- 10

# Set simul parameters
bd_params <- make.bd.params(speciation = runif,
                            extinction = runif,
                            joint      = TRUE)

# Generate simulations of traits
Sp_Traits_BM <- list()

for(i in 1:length(n.trait)) { 
  
  # the Brownian motion (BM) where the variance at t+1 is independent of the variance at t-1 (which effectively results in a constant increase in trait variance through time).
  trait_processes_BM <- make.traits(process = BM.process, n = n.trait[i])
  
  Sp_Traits_BM[[i]] <- replicate(n.sim, 
                                 dads(bd.params      = bd_params,
                                      stop.rule      = list(max.living = n.sp),
                                      traits         = trait_processes_BM,
                                      null.error     = 100), 
                                 simplify = FALSE)
  
}

names(Sp_Traits_BM) <- paste0(rep("Trait_", length (n.trait)), n.trait)

saveRDS(Sp_Traits_BM, file = "tree_traits_BM.rds")

# Analyse brownian motion simuls ------------------------------------------------

Sp_Traits_BM <- readRDS("tree_traits_BM.rds")

n.trait <- c(1,2,4,8)
comm_list <- list()

# Generate simuls (first loops over number of traits, next on species)
for(i in 1 : length(Sp_Traits_BM)) {  
  
  message(paste("Selecting trait n° ", n.trait[i]))
  
  Sp_Traits_i <- Sp_Traits_BM[[i]]
  
  #Loop ver traits
  for(j in 1 : length(Sp_Traits_i)) { 
    
    traits_j <- drop.fossil.dads(Sp_Traits_i[[j]])
    
    # Generate the random communities for the run
    if(j == 1){ 
      
      comm <- sample(c(rep(0,90), rep(1,10))) #Community with 10 species
      
      for(k in seq(from = 20, to = 100, by = 10) ) 
        comm <- rbind(comm, sample(c(rep(0, 100 - k), rep(1, k))) ) #Communities with 20 to 100 species
      
      rownames(comm) <- paste0(rep("Comm",10), 1:10)
      colnames(comm) <- rownames(traits_j)
      
      comm_list[[i]] <- comm # store
    } 
    
    # select traits
    traits_j <- traits_j$data[rownames(traits_j$data) %in% traits_j$tree$tip.label, , drop = FALSE]
    #nrow(traits_j) #Should be 100!
    
    message(paste0("-------- Simulation ", j , " out of ", length(Sp_Traits_i)))
    
    message("------ Estimating nj trees ------")
    
    #Calculate stats with nj Trees
    tree_nj <- tree.build(traits_j, distance = "gower", func = "nj")
    
    r_tree_nj <- alpha(comm = comm, tree = tree_nj) #richness
    d_tree_nj <- dispersion(comm = comm, tree = tree_nj) #divergence
    e_tree_nj <- evenness(comm = comm, tree = tree_nj, method = "expected", func = "camargo") #regularity
    
    message("------ Estimating UPGMA trees ------")
    
    #Calculate stats with nj Trees
    tree_upgma <- tree.build(traits_j, distance = "gower", func = "upgma")
    
    r_tree_upgma <- alpha(comm = comm, tree = tree_upgma) #richness
    d_tree_upgma <- dispersion(comm = comm, tree = tree_upgma) #divergence
    e_tree_upgma <- evenness(comm = comm, tree = tree_upgma) #regularity
    
    message("------ Estimating hypervolumes ------")
    
    #Calculate stats with hypervolumes
    hv <- quiet(kernel.build(comm = comm, trait = traits_j, method = "box", cores = 10))
    
    r_hv <- kernel.alpha(hv) ; names(r_hv) <- NULL #richness 
    d_hv <- kernel.dispersion(hv) ; names(d_hv) <- NULL #divergence
    e_hv <- quiet(kernel.evenness(hv)) ; names(e_hv) <- NULL #regularity
    
    #Store the results
    
    if(j > 1) {
      
      n_trait         <- c(n_trait, rep(ncol(traits_j),10))
      n_species       <- c(n_species, seq(from = 10, to = 100, by = 10))
      simulation      <- c(simulation, rep(j, nrow(comm)))
      ric_tree_nj     <- c(ric_tree_nj, r_tree_nj)
      div_tree_nj     <- c(div_tree_nj, d_tree_nj)
      reg_tree_nj     <- c(reg_tree_nj, e_tree_nj)
      ric_tree_upgma  <- c(ric_tree_upgma, r_tree_upgma)
      div_tree_upgma  <- c(div_tree_upgma, d_tree_upgma) 
      reg_tree_upgma  <- c(reg_tree_upgma, e_tree_upgma) 
      ric_hv          <- c(ric_hv, r_hv)
      div_hv          <- c(div_hv, d_hv)
      reg_hv          <- c(reg_hv, e_hv)
      
    } else { 
      
      n_trait         <- rep(ncol(traits_j),10)
      n_species       <- seq(from = 10, to = 100, by = 10)
      simulation      <- rep(j, nrow(comm))
      ric_tree_nj     <- r_tree_nj 
      div_tree_nj     <- d_tree_nj 
      reg_tree_nj     <- e_tree_nj 
      ric_tree_upgma  <- r_tree_upgma
      div_tree_upgma  <- d_tree_upgma 
      reg_tree_upgma  <- e_tree_upgma 
      ric_hv          <- r_hv
      div_hv          <- d_hv
      reg_hv          <- e_hv
      
    }
  }
  
  if(i > 1) {
    results <- rbind(results, data.frame(
      process = rep("BM", length(simulation)),
      n_trait,
      n_species,
      simulation,
      ric_tree_nj,
      div_tree_nj,
      reg_tree_nj,
      ric_tree_upgma,
      div_tree_upgma, 
      reg_tree_upgma,  
      ric_hv,   
      div_hv,  
      reg_hv))
  } else {
    
    results <- data.frame(process = rep("BM", length(simulation)),
                          n_trait,
                          n_species,
                          simulation,
                          ric_tree_nj,
                          div_tree_nj,
                          reg_tree_nj,
                          ric_tree_upgma,
                          div_tree_upgma, 
                          reg_tree_upgma,  
                          ric_hv,   
                          div_hv,  
                          reg_hv)
  }
}

## Save
saveRDS(results, file = "results_BM.rds")
saveRDS(comm_list, file = "random_communities.rds")

## ------------
## Checking results
## ------------

## Loading R packages for data analysis

library("dplyr")
library("tidyverse")
library("GGally")
library("gridExtra")
library("cowplot")

# Brownian motion ---------------------------------------------------------

## Checking correlation among Richness, Divergence and Regularity in nj-trees

results_BM <- readRDS(file = "results_BM.rds")

results_BM <- results_BM %>% mutate_at(vars(n_trait), as.factor)
levels(results_BM$n_trait) <- c("1 Trait", "2 Traits", "4 Traits", "8 Traits")

colors_manual = c("blue", "orange", "grey10", "purple")

# To enhance visualization...
results_BM$ric_hv <- log(results_BM$ric_hv+1)
results_BM$reg_hv <- sqrt(asin(results_BM$reg_hv))

# Check correlations (as scatterplots), variable distribution, and print correlation coefficient
(plot1 <- results_BM %>% 
    dplyr::select(n_trait, n_species, ric_tree_nj, div_tree_nj, reg_tree_nj) %>% 
    GGally::ggpairs(aes(colour = n_trait, fill = n_trait, alpha = 0.6),
                    columns = 3:5,
                    upper = list(continuous = GGally::wrap(ggally_cor, method = "spearman", stars = F)),
                    title = "Neighbour-joining tree",
                    columnLabels = c(
                      "Richness",
                      "Divergence",
                      "Regularity")) + 
    scale_color_manual(values = rev(colors_manual))+
    scale_fill_manual(values = rev(colors_manual))+
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), units = "cm"),
      strip.text=element_text(color = "grey10", face = "bold", size = 10),
      strip.background = element_rect(colour=NA, fill="white")))

# Is there cor with upgma trees and hypervolumes?
strip_label_custom <- c("nj tree", "upgma tree", "hypervolume")

(plot2 <- results_BM %>% 
    dplyr::select(n_trait, ric_tree_nj, ric_tree_upgma, ric_hv) %>% 
    GGally::ggpairs(aes(colour = n_trait, fill = n_trait, alpha = 0.6),
                    columns = 2:4,
                    upper = list(continuous = GGally::wrap(ggally_cor, method = "spearman", stars = F)),
                    title = "Richness",
                    columnLabels = strip_label_custom) + 
    scale_color_manual(values = rev(colors_manual))+
    scale_fill_manual(values = rev(colors_manual))+
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), units = "cm"),
      strip.text=element_text(color = "grey10", face = "bold", size = 10),
      strip.background = element_rect(colour=NA, fill="white")))

(plot3 <- results_BM %>% 
    dplyr::select(n_trait, div_tree_nj, div_tree_upgma, div_hv) %>% 
    GGally::ggpairs(aes(colour = n_trait, fill = n_trait, alpha = 0.6),
                    columns = 2:4,
                    upper = list(continuous = GGally::wrap(ggally_cor,method = "spearman", stars = F)),
                    title = "Divergence",
                    columnLabels = strip_label_custom) + 
    scale_color_manual(values = rev(colors_manual))+
    scale_fill_manual(values = rev(colors_manual))+
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), units = "cm"),
      strip.text=element_text(color = "grey10", face = "bold", size = 10),
      strip.background = element_rect(colour=NA, fill="white")))

(plot4 <- results_BM %>% 
    dplyr::select(n_trait, reg_tree_nj, reg_tree_upgma, reg_hv) %>% 
    GGally::ggpairs(aes(colour = n_trait, fill = n_trait, alpha = 0.6),
                    columns = 2:4,
                    upper = list(continuous = GGally::wrap(ggally_cor, method = "spearman", stars = F)),
                    title = "Regularity",
                    columnLabels = strip_label_custom) + 
    scale_color_manual(values = rev(colors_manual))+
    scale_fill_manual(values = rev(colors_manual))+
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), units = "cm"),
      strip.text=element_text(color = "grey10", face = "bold", size = 10),
      strip.background = element_rect(colour=NA, fill="white")))

# Saving plots
ggsave(plot = cowplot::plot_grid(
  ggmatrix_gtable(plot1),
  ggmatrix_gtable(plot2),
  ggmatrix_gtable(plot3),
  ggmatrix_gtable(plot4),
  nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")), filename = "Figure_2.pdf",width = 38, height = 32, units = "cm")


########################### Quality of the functional space (Appendix S1) ##################################################################


########################### Preliminary notes #############################################
## Note: Data - 10 replicates of 2 Processes (BM or OU) x n_traits (1, 2, 4 and 8 traits)
## 10 communities (10, 20, ..., 100 species) for each combination

## Note: Quality - A tree was rebuilt for each community, using the traits of species present in the community.
# This tree (NOT the initial tree with all the species) was used to evaluate the quality
# of the functional space

## BAT quality functions were adapted to accept several communities at once.
# tree.quality.com
# pcoa.quality.com

########################### Tree quality with community data function #####################
tree.quality.com <- function (comm, traits, tree.builder = "nj") {
  comm <- comm # community data
  traits <- traits # trait data
  tree.builder <- tree.builder # either nj or upgma
  
  ## cycle through all communities
  res <- NULL
  for (i in 1: nrow(comm)) {
    present <- which (comm[i,] > 0)
    tr_present <- as.data.frame(traits [present,])
    distance <- gower (tr_present)
    tree <- tree.build (tr_present, distance = "gower", func = paste(tree.builder))
    quality <- msd(distance, as.dist(cophenetic (tree)))
    res <- c(res, quality)
    
  }
  res
  
}

# tree.quality.com(comm, run1_traits_extant, "upgma")

########################### PCOA quality function #########################
pcoa.quality.com <- function (comm, traits, k = 2) {
  comm <- comm # community data
  traits <- traits # trait data
  k = k # number of axes
  
  ## cycle through all communities
  res = NULL
  for (i in 1: nrow(comm)) {
    present <- which (comm[i,] > 0)
    tr_present <- as.data.frame(traits [present,])
    distance <- gower (tr_present)
    pcoa <- cmdscale(distance, k= k , add = T)
    quality <- msd(distance, dist(pcoa$points))
    res <- c(res, quality)
    
  }
  res
  
}

# pcoa.quality.com(comm, run1_traits_extant, 2)

########################### libraries ######################
# Libraries
library(ape)
library(dads)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(BAT)

########################### Data import ##################################
# to read
tree_traits_OU <- readRDS("C:/path/tree_traits_OU.rds")
tree_traits_BM <- readRDS("C:/path/tree_traits_BM.rds")
random_communities <- readRDS("C:/path/random_communities.rds")

# to read
Sp_Traits_BM <- readRDS("tree_traits_BM.rds")
Sp_Traits_OU <- readRDS("tree_traits_OU.rds")



########################### Results :: NJ vs UPGMA #####################################################
########################### BM process :: 1 trait  #######################
# a list of 10 replicates of 1 trait
trait1 <- Sp_Traits_BM[[1]] 

# 10 random communities
comm1 <- random_communities[[1]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa <- NULL
for (i in 1:length(trait1)){
  
  # cycle through the 10 replicates
  run1 <- trait1[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  pcoa <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 1)
  res_pcoa <- c(res_pcoa, pcoa)
  
}

resBM_1trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa", 100)), n_traits = rep(1, 300), n_species = rep(seq(10,100, 10), 30), quality = c(res_nj, res_upgma, res_pcoa))

########################### BM process :: 2 traits #######################
# a list of 10 replicates of 2 traits
trait2 <- Sp_Traits_BM[[2]] 

# 10 random communities
comm1 <- random_communities[[2]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa <- NULL
for (i in 1:length(trait2)){
  
  # cycle through the 10 replicates
  run1 <- trait2[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  pcoa <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 1)
  res_pcoa <- c(res_pcoa, pcoa)
  
  
}

resBM_2trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa-1axis", 100)), n_traits = rep(2, 300), n_species = rep(seq(10,100, 10), 30), quality = c(res_nj, res_upgma, res_pcoa))


########################### BM process :: 4 traits #######################
# a list of 10 replicates of 4 trait
trait4 <- Sp_Traits_BM[[3]] 

# 10 random communities
comm1 <- random_communities[[3]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa <- NULL
for (i in 1:length(trait4)){
  
  # cycle through the 10 replicates
  run1 <- trait4[[i]]
  
  ## drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  pcoa <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 2)
  res_pcoa <- c(res_pcoa, pcoa)
  
  
}

resBM_4trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa-2 axes", 100)), n_traits = rep(4, 300), n_species = rep(seq(10,100, 10), 30), quality = c(res_nj, res_upgma, res_pcoa))

########################### BM process :: 8 traits #######################
# a list of 10 replicates of 8 trait
trait8 <- Sp_Traits_BM[[4]] 

# 10 random communities
comm1 <- random_communities[[4]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa2 <- NULL
res_pcoa4 <- NULL
for (i in 1:length(trait8)){
  
  # cycle through the 10 replicates
  run1 <- trait8[[i]]
  
  ## drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  ## PCOA with 2 axes
  pcoa2 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 2)
  res_pcoa2 <- c(res_pcoa2, pcoa2)
  
  ## PCOA with 4 axes
  pcoa4 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 4)
  res_pcoa4 <- c(res_pcoa4, pcoa4)
  
  
}

resBM_8trait <- data.frame (process = rep("BM", 400), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa-2 axes", 100), rep("pcoa-4 axes", 100)), n_traits = rep(8, 400), n_species = rep(seq(10,100, 10), 40), quality = c(res_nj, res_upgma, res_pcoa2, res_pcoa4))


########################### OU process :: 1 trait  #######################
# a list of 10 replicates of 1 trait
trait1 <- Sp_Traits_OU[[1]] 

# 10 random communities
comm1 <- random_communities[[1]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa <- NULL
for (i in 1:length(trait1)){
  
  # cycle through the 10 replicates
  run1 <- trait1[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  pcoa <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 1)
  res_pcoa <- c(res_pcoa, pcoa)
  
  
}

resOU_1trait = data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa", 100)), n_traits = rep(1, 300), n_species = rep(seq(10,100, 10), 30), quality = c(res_nj, res_upgma, res_pcoa))

########################### OU process :: 2 traits #######################
# a list of 10 replicates of 2 traits
trait2 <- Sp_Traits_OU[[2]] 

# 10 random communities
comm1 <- random_communities[[2]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa <- NULL
for (i in 1:length(trait2)){
  
  # cycle through the 10 replicates
  run1 <- trait2[[i]]
  
  ## drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  pcoa <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 1)
  res_pcoa <- c(res_pcoa, pcoa)
  
  
}

resOU_2trait <- data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa", 100)), n_traits = rep(2, 300), n_species = rep(seq(10,100, 10), 30), quality = c(res_nj, res_upgma, res_pcoa))


########################### OU process :: 4 traits #######################
# a list of 10 replicates of 4 traits
trait4 <- Sp_Traits_OU[[3]] 

# 10 random communities
comm1 <- random_communities[[3]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa <- NULL
for (i in 1:length(trait4)){
  
  # cycle through the 10 replicates
  run1 <- trait4[[i]]
  
  ## drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  pcoa <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 2)
  res_pcoa <- c(res_pcoa, pcoa)
  
  
}

resOU_4trait <- data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa-2 axes", 100)), n_traits = rep(4, 300), n_species = rep(seq(10,100, 10), 30), quality = c(res_nj, res_upgma, res_pcoa))


########################### OU process :: 8 traits #######################
# a list of 10 replicates of 8 traits
trait8 <- Sp_Traits_OU[[4]] 

# 10 random communities
comm1 <- random_communities[[4]]

res_nj <- NULL
res_upgma <- NULL
res_pcoa2 <- NULL
res_pcoa4 <- NULL
for (i in 1:length(trait8)){
  
  # cycle through the 10 replicates
  run1 <- trait8[[i]]
  
  ## Nota: adicionei isto pra dar as 100 espécies (está no script Simulation V5)
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- tree.quality.com (comm = comm1, traits = traits_j, tree.builder = "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  ## PCOA with 2 axes
  pcoa2 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 2)
  res_pcoa2 <- c(res_pcoa2, pcoa2)
  
  ## PCOA with 4 axes
  pcoa4 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 4)
  res_pcoa4 <- c(res_pcoa4, pcoa4)
  
  
}

resOU_8trait <- data.frame (process = rep("OU", 400), method = c(rep("nj", 100), rep("upgma", 100), rep("pcoa-2 axes", 100), rep("pcoa-4 axes", 100)), n_traits = rep(8, 400), n_species = rep(seq(10,100, 10), 40), quality = c(res_nj, res_upgma, res_pcoa2, res_pcoa4))


########################### Plots NJ vs UPGMA (Appendix S1: Fig. S1.1) ############################################
## OU process
data <- rbind (subset (resOU_1trait, method %in% c("nj", "upgma")), 
               subset (resOU_2trait, method %in% c("nj", "upgma")),
               subset (resOU_4trait, method %in% c("nj", "upgma")), 
               subset (resOU_8trait, method %in% c("nj", "upgma"))
)


# Grouped
png ("All_OU2.png", width = 900, height = 600)


data %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=quality, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent") +
  facet_wrap(~ n_traits, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("OU process") +
  xlab("Number of species") +
  ylab("Quality of the functional space") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,1)

dev.off()


## BM process
data <- rbind (subset (resBM_1trait, method %in% c("nj", "upgma")), 
               subset (resBM_2trait, method %in% c("nj", "upgma")),
               subset (resBM_4trait, method %in% c("nj", "upgma")), 
               subset (resBM_8trait, method %in% c("nj", "upgma"))
)

# Grouped
png ("All_BM2.png", width = 900, height = 600)

data %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=quality, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent") +
  facet_wrap(~ n_traits, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("BM process") +
  xlab("Number of species") +
  ylab("Quality of the functional space") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,1)

dev.off()





########################### Results :: PCOA 2 to 8 axes (8 traits) ######################################
#### Quality results for a sequence of 2 to 8 axes of a PCOA for datasets with 8 traits


########################### BM process :: 8 traits #####################################
# a list of 10 replicates of 8 trait
trait8 <- Sp_Traits_BM[[4]] 

# 10 random communities
comm1 <- random_communities[[4]]

res_pcoa2 <- NULL
res_pcoa3 <- NULL
res_pcoa4 <- NULL
res_pcoa5 <- NULL
res_pcoa6 <- NULL
res_pcoa7 <- NULL
res_pcoa8 <- NULL
for (i in 1:length(trait8)){
  
  # cycle through the 10 replicates
  run1 <- trait8[[i]]
  
  ## drop fosil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  
  ## PCOA with 2 axes
  pcoa2 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 2)
  res_pcoa2 <- c(res_pcoa2, pcoa2)
  
  ## PCOA with 3 axes
  pcoa3 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 3)
  res_pcoa3 <- c(res_pcoa3, pcoa3)
  
  ## PCOA with 4 axes
  pcoa4 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 4)
  res_pcoa4 <- c(res_pcoa4, pcoa4)
  
  ## PCOA with 5 axes
  pcoa5 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 5)
  res_pcoa5 <- c(res_pcoa5, pcoa5)
  
  ## PCOA with 6 axes
  pcoa6 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 6)
  res_pcoa6 <- c(res_pcoa6, pcoa6)
  
  ## PCOA with 7 axes
  pcoa7 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 7)
  res_pcoa7 <- c(res_pcoa7, pcoa7)
  
  ## PCOA with 8 axes
  pcoa8 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 8)
  res_pcoa8 <- c(res_pcoa8, pcoa8)
  
  
}

resBM_8trait_PCOA <- data.frame (process = rep("BM", 700), method = c(rep("pcoa-2 axes", 100), rep("pcoa-3 axes", 100), rep("pcoa-4 axes", 100), rep("pcoa-5 axes", 100), rep("pcoa-6 axes", 100), rep("pcoa-7 axes", 100), rep("pcoa-8 axes", 100)), n_traits = rep(8, 700), n_species = rep(seq(10, 100, 10), 70), quality = c(res_pcoa2, res_pcoa3, res_pcoa4, res_pcoa5, res_pcoa6, res_pcoa6, res_pcoa8))

data <- resBM_8trait_PCOA

# Grouped
png ("Plot_8T_BM_PCOA.png", width = 900, height = 1200)

data %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=quality, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  # facet_grid(~fct_relevel(method, "pcoa-2 axes", "pcoa-3 axes", "pcoa-4 axes", "pcoa-5 axes", "pcoa-6 axes", "pcoa-7 axes", "pcoa-8 axes")) +
  facet_wrap(~fct_relevel(method, "pcoa-2 axes", "pcoa-3 axes", "pcoa-4 axes", "pcoa-5 axes", "pcoa-6 axes", "pcoa-7 axes", "pcoa-8 axes"), scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("BM process: 8 traits") +
  xlab("Number of species") +
  ylab("Quality of the functional space") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,1) 


dev.off()

########################### OU process :: 8 traits #############################
# a list of 10 replicates of 8 trait
trait8 <- Sp_Traits_OU[[4]] 

# 10 random communities
comm1 <- random_communities[[4]]

res_pcoa2 <- NULL
res_pcoa3 <- NULL
res_pcoa4 <- NULL
res_pcoa5 <- NULL
res_pcoa6 <- NULL
res_pcoa7 <- NULL
res_pcoa8 <- NULL
for (i in 1:length(trait8)){
  
  # cycle through the 10 replicates
  run1 <- trait8[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  
  
  ## PCOA with 2 axes
  pcoa2 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 2)
  res_pcoa2 <- c(res_pcoa2, pcoa2)
  
  ## PCOA with 3 axes
  pcoa3 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 3)
  res_pcoa3 <- c(res_pcoa3, pcoa3)
  
  ## PCOA with 4 axes
  pcoa4 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 4)
  res_pcoa4 <- c(res_pcoa4, pcoa4)
  
  ## PCOA with 5 axes
  pcoa5 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 5)
  res_pcoa5 <- c(res_pcoa5, pcoa5)
  
  ## PCOA with 6 axes
  pcoa6 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 6)
  res_pcoa6 <- c(res_pcoa6, pcoa6)
  
  ## PCOA with 7 axes
  pcoa7 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 7)
  res_pcoa7 <- c(res_pcoa7, pcoa7)
  
  ## PCOA with 8 axes
  pcoa8 <- pcoa.quality.com(comm = comm1, traits = traits_j, k = 8)
  res_pcoa8 <- c(res_pcoa8, pcoa8)
  
  
}

resOU_8trait_PCOA <- data.frame (process = rep("OU", 700), method = c(rep("pcoa-2 axes", 100), rep("pcoa-3 axes", 100), rep("pcoa-4 axes", 100), rep("pcoa-5 axes", 100), rep("pcoa-6 axes", 100), rep("pcoa-7 axes", 100), rep("pcoa-8 axes", 100)), n_traits = rep(8, 700), n_species = rep(seq(10, 100, 10), 70), quality = c(res_pcoa2, res_pcoa3, res_pcoa4, res_pcoa5, res_pcoa6, res_pcoa6, res_pcoa8))

data <- resOU_8trait_PCOA

# Grouped
png ("Plot_8T_OU_PCOA.png", width = 900, height = 1200)

data %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=quality, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  # facet_grid(~fct_relevel(method, "pcoa-2 axes", "pcoa-3 axes", "pcoa-4 axes", "pcoa-5 axes", "pcoa-6 axes", "pcoa-7 axes", "pcoa-8 axes")) +
  facet_wrap(~fct_relevel(method, "pcoa-2 axes", "pcoa-3 axes", "pcoa-4 axes", "pcoa-5 axes", "pcoa-6 axes", "pcoa-7 axes", "pcoa-8 axes"), scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("OU process: 8 traits") +
  xlab("Number of species") +
  ylab("Quality of the functional space") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,1) 


dev.off()


########################### Sensitivity to outliers (Appendix S2) #######################

########################### Preliminary notes ##################
# Note: Data - 10 replicates of 2 Processes (BM or OU) x n_traits (1, 2, 4 and 8 traits)
## 10 communities (10, 20, ..., 100 species) for each combination

# Outlier: species with maximum average distance to all other species 

## Note: 
# - A tree was built for each community and used to measure FD.
# - The same tree was used to measure FD without the outlier
# - % CHANGE = (FD_all - FD_w/outlier) / FD_all


########################### Function to calculate % change in FD after removing outliers ####
FD.outlier = function (comm, traits, method = "nj") {
  comm = comm
  traits = traits
  method = method
  
  ## cycle through all communities
  res = NULL
  for (i in 1: nrow(comm)) {
    present <- which (comm[i,] > 0)
    tr_present <- as.data.frame(traits [present,])
    
    # distance matrix with all species
    distance <- gower (tr_present)
    
    # Maximum value of the average distance from one species to all other
    dd <- rowMeans(as.matrix(distance))
    outl <- which (dd == max(dd)) ## outlier
    
    # method switch
    method <- match.arg(method, c("nj","upgma", "hyper"))
    
    switch (method, nj = {
      
      ## all species
      tree_all <- tree.build (tr_present, distance = "gower", func = "nj")
      
      # matrix with comm i (first row) and comm i without the outlier (second row)
      mat <- matrix (1, nrow = 2, ncol=length(present))
      mat[2, outl] <- 0 # set the outlier with a zero (absence)
      
      # FD
      FD = alpha (mat, tree_all)
      
      ## sensitivity to outlier
      sst <- (FD[1,]-FD[2,]) / FD[1,] * 100
      names (sst) <- NULL
      
      
    }, upgma = {
      
      ## all species
      tree_all <- tree.build (tr_present, distance = "gower", func = "upgma")
      
      # matrix with comm i (first row) and comm i without the outlier (second row)
      mat <- matrix (1, nrow = 2, ncol=length(present))
      mat[2, outl] <- 0 # set the outlier with a zero (absence)
      
      # FD
      FD = alpha (mat, tree_all)
      
      ## sensitivity to outlier
      sst <- (FD[1,]-FD[2,]) / FD[1,] * 100
      names (sst) <- NULL
      
    }, hyper = {
      
      # matrix with comm i (first row) and comm i without the outlier (second row)
      mat <- matrix (1, nrow = 2, ncol=length(present))
      mat[2, outl] <- 0 # set the outlier with a zero (absence)
      rownames (mat) <- c("all", "out")
      
      # standardize traits into the range 0-1
      # kernel.build is not doing the standardization
      library(vegan)
      trait_range <- decostand(tr_present, "range")
      
      hv <- kernel.build(comm = mat, trait = trait_range, method = "box", abun = F, cores = 10)
      FD <- kernel.alpha(hv) ; names(hv) <- NULL 
      
      ## sensitivity to outlier
      sst <- (FD[1]-FD[2]) / FD[1] * 100
      names (sst) <- NULL
      
    })
    
    res <- c(res, sst)
    
  }
  res
  
}

# FD.outlier (comm1, traits_j, "nj")
# FD.outlier (comm1, traits_j, "upgma")
# FD.outlier (comm1, traits_j, "hyper")

########################### libraries ######################
# Libraries
library(parallel)
library(hypervolume)
library(ape)
library(dads)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(parallel)
library(BAT)

########################### Data import ##################################
# to read
tree_traits_OU <- readRDS("C:/path/tree_traits_OU.rds")
tree_traits_BM <- readRDS("C:/path/tree_traits_BM.rds")
random_communities <- readRDS("C:/path/random_communities.rds")


# to read
Sp_Traits_BM <- readRDS("tree_traits_BM.rds")
Sp_Traits_OU <- readRDS("tree_traits_OU.rds")

########################### BM process :: 1 trait  #######################
# a list of 10 replicates of 1 trait
trait1 <- Sp_Traits_BM[[1]] 

# 10 random communities
comm1 <- random_communities[[1]]


res_nj <- NULL
res_upgma <- NULL
# res_hyper <- NULL
for (i in 1:length(trait1)){
  
  # cycle through the 10 replicates
  run1 <- trait1[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  # hyper <- FD.outlier (comm1, traits_j, "hyper") # Hypervolumes are not possible when nº of axes = 1 
  # res_hyper <- c(res_hyper, hyper)
  
}

# outBM_1trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(1, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))
outBM_1trait <- data.frame (process = rep("BM", 200), method = c(rep("nj", 100), rep("upgma", 100)), n_traits = rep(1, 200), n_species = rep(seq(10,100, 10), 20), change = c(res_nj, res_upgma))


data = outBM_1trait

# Grouped
png ("out_1T_BM.png", width = 600, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("BM process: 1 trait") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### BM process :: 2 traits  #######################
# a list of 10 replicates of 2 trait
trait2 <- Sp_Traits_BM[[2]] 

# 10 random communities
comm1 <- random_communities[[2]]


res_nj <- NULL
res_upgma <- NULL
res_hyper <- NULL
for (i in 1:length(trait2)){
  
  # cycle through the 10 replicates
  run1 <- trait2[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  hyper <- FD.outlier (comm1, traits_j, "hyper")
  res_hyper <- c(res_hyper, hyper)
  
}

outBM_2trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(2, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))

data = outBM_2trait

# Grouped
png ("out_2T_BM.png", width = 900, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("BM process: 2 traits") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### BM process :: 4 traits  #######################
# a list of 10 replicates of 4 traits
trait4 <- Sp_Traits_BM[[3]] 

# 10 random communities
comm1 <- random_communities[[3]]


res_nj <- NULL
res_upgma <- NULL
res_hyper <- NULL
for (i in 1:length(trait4)){
  
  # cycle through the 10 replicates
  run1 <- trait4[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  hyper <- FD.outlier (comm1, traits_j, "hyper")
  res_hyper <- c(res_hyper, hyper)
  
}

outBM_4trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(4, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))

data <- outBM_4trait

# Grouped
png ("out_4T_BM.png", width = 900, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("BM process: 4 traits") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### BM process :: 8 traits  #######################
# a list of 10 replicates of 8 traits
trait8 <- Sp_Traits_BM[[4]] 

# 10 random communities
comm1 <- random_communities[[4]]


res_nj <- NULL
res_upgma <- NULL
res_hyper <- NULL
for (i in 1:length(trait8)){
  
  # cycle through the 10 replicates
  run1 <- trait8[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  hyper <- FD.outlier (comm1, traits_j, "hyper")
  res_hyper <- c(res_hyper, hyper)
  
}

outBM_8trait <- data.frame (process = rep("BM", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(8, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))

data <- outBM_8trait

# Grouped
png ("out_8T_BM.png", width = 900, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("BM process: 8 traits") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### OU process :: 1 trait  #######################
# a list of 10 replicates of 1 trait
trait1 <- Sp_Traits_OU[[1]] 

# 10 random communities
comm1 <- random_communities[[1]]


res_nj <- NULL
res_upgma <- NULL
# res_hyper <- NULL
for (i in 1:length(trait1)){
  
  # cycle through the 10 replicates
  run1 <- trait1[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  # hyper <- FD.outlier (comm1, traits_j, "hyper") # hypervolumes are not possible when nº of axes = 1
  # res_hyper <- c(res_hyper, hyper)
  
}

# outOU_1trait <- data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(1, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))
outOU_1trait <- data.frame (process = rep("OU", 200), method = c(rep("nj", 100), rep("upgma", 100)), n_traits = rep(1, 200), n_species = rep(seq(10,100, 10), 20), change = c(res_nj, res_upgma))


data <- outOU_1trait

# Grouped
png ("out_1T_OU.png", width = 600, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("OU process: 1 trait") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### OU process :: 2 traits  #######################
# a list of 10 replicates of 2 traits
trait2 <- Sp_Traits_OU[[2]] 

# 10 random communities
comm1 <- random_communities[[2]]


res_nj <- NULL
res_upgma <- NULL
res_hyper <- NULL
for (i in 1:length(trait2)){
  
  # cycle through the 10 replicates
  run1 <- trait2[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj = FD.outlier (comm1, traits_j, "nj")
  res_nj = c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  hyper <- FD.outlier (comm1, traits_j, "hyper")
  res_hyper <- c(res_hyper, hyper)
  
}

outOU_2trait <- data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(2, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))

data <- outOU_2trait

# Grouped
png ("out_2T_OU.png", width = 900, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("OU process: 2 traits") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### OU process :: 4 traits  #######################
# a list of 10 replicates of 4 traits
trait4 <- Sp_Traits_OU[[3]] 

# 10 random communities
comm1 <- random_communities[[3]]


res_nj <- NULL
res_upgma <- NULL
res_hyper <- NULL
for (i in 1:length(trait4)){
  
  # cycle through the 10 replicates
  run1 <- trait4[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  hyper <- FD.outlier (comm1, traits_j, "hyper")
  res_hyper <- c(res_hyper, hyper)
  
}

outOU_4trait <- data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(4, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))

data <- outOU_4trait

# Grouped
png ("out_4T_OU.png", width = 900, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("OU process: 4 traits") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()

########################### OU process :: 8 traits  #######################
# a list of 10 replicates of 8 traits
trait8 <- Sp_Traits_OU[[4]] 

# 10 random communities
comm1 <- random_communities[[4]]


res_nj <- NULL
res_upgma <- NULL
res_hyper <- NULL
for (i in 1:length(trait8)){
  
  # cycle through the 10 replicates
  run1 <- trait8[[i]]
  
  # drop fossil species
  run1 <- drop.fossil.dads(run1)
  
  # Note that the number of species is higher than 100, because it also includes the fossil species used to generate the current community.
  # To get the 100 extant species, run this:
  traits_j <- run1$data[rownames(run1$data) %in% run1$tree$tip.label, , drop = FALSE]
  
  nj <- FD.outlier (comm1, traits_j, "nj")
  res_nj <- c(res_nj, nj)
  
  upgma <- FD.outlier (comm1, traits_j, "upgma")
  res_upgma <- c(res_upgma, upgma)
  
  hyper = FD.outlier (comm1, traits_j, "hyper")
  res_hyper <- c(res_hyper, hyper)
  
}

outOU_8trait <- data.frame (process = rep("OU", 300), method = c(rep("nj", 100), rep("upgma", 100), rep("hyper", 100)), n_traits = rep(8, 300), n_species = rep(seq(10,100, 10), 30), change = c(res_nj, res_upgma, res_hyper))

data <- outOU_8trait

# Grouped
png ("out_8T_OU.png", width = 900, height = 600)

data_new <- data                              # Replicate data
data_new$method <- factor(data_new$method,      # Reordering group factor levels
                          levels = c("nj", "upgma", "hyper"))

data_new %>%
  mutate(n_species = factor(n_species)) %>%
  ggplot(aes(fill=method, y=change, x=n_species)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", show.legend = FALSE) +
  facet_wrap(~method, scale="free") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  ggtitle("OU process: 8 traits") +
  xlab("Number of species") +
  ylab("% Change") +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank()) +
  ylim(0,100)

dev.off()


#########################################################
####Kiwi NJ Analyses###################################
#############################################################

library(dplyr)

########################################
#########Load and format data############
#########################################

dat4 <- readRDS("kiwi_data.RDS")

#vector of five kiwi species names
kiwi <- c("Apteryx australis",
          "Apteryx haastii",
          "Apteryx mantelli",
          "Apteryx owenii",
          "Apteryx rowi")

##select and format traits (logged then scaled)
datJR <- dplyr::select(dat4, -species)
datJR2 <- apply(datJR, 2, function(x) scale(log(x))) %>%
  as.data.frame()
rownames(datJR2) <- dat4$species

#Generate the PCoA
ee_dist <- dist(datJR2)
trai.pcoa <- ape::pcoa(ee_dist)
TRAITS <- trai.pcoa$vectors[,1:5]#first 5 axes
TRAITS_ALL <- trai.pcoa$vectors#all axes

####################################################
##calculate distance and build the two dendrograms
#####################################################
eu_dist <- dist(TRAITS)#five PCoA axes
eu_dist_all <- dist(TRAITS_ALL)#all PCoA axes

#upgma
up_tree <- tree.build(TRAITS, distance = "euclidean", 
                      func = "upgma")
up_tree <- ape:::as.phylo(up_tree)

#nj
nj_tree <- tree.build(TRAITS, distance = "euclidean",
                      func = "nj")

nj_tree <- ape:::as.phylo(nj_tree)

######################################################
####BUILD HyperVolume & convex hull############
#################################################

##Create a presence-absence matrices with two assemblages:
#with and without kiwis (rows = sites, cols = species)
hyper_PA <- matrix(1, ncol = nrow(TRAITS), nrow = 2)
colnames(hyper_PA) <- rownames(TRAITS)
hwm <- which(colnames(hyper_PA) %in% kiwi)
hyper_PA[2, hwm] <- 0 #set kiwis to zero for one row
rownames(hyper_PA) <- c("withKiwi", "withoutKiwi")

###Build convex hull and hypervolume
hyper0 <- hyper.build(TRAITS, distance = "euclidean", axes = 0)

hyper <- kernel.build(comm = hyper_PA,
                      trait = hyper0, distance = "euclidean",
                      abund = FALSE, axes = 0)

hull <- hull.build(comm = hyper_PA,
                   trait = hyper0, distance = "euclidean",
                   axes = 0)

##Work out % change in FD with and without the kiwis 
#(Alpha (with kiwis) - Alpha (excluding kiwis) ) / alpha(with kiwis)
#UPGMA
up_alpha <- BAT::alpha(comm = hyper_PA, tree = up_tree)
((up_alpha[1] - up_alpha[2]) / up_alpha[1]) * 100
#NJ
nj_alpha <- BAT::alpha(comm = hyper_PA, tree = nj_tree)
((nj_alpha[1] - nj_alpha[2]) / nj_alpha[1]) * 100
#hyper
hyper_alpha <- kernel.alpha(hyper)
((hyper_alpha[1] - hyper_alpha[2]) / hyper_alpha[1]) * 100
#hull
hull_alpha <- hull.alpha(hull)
((hull_alpha[1] - hull_alpha[2]) / hull_alpha[1]) * 100

#######################################
###Tree and hyperspace Quality##########
##################################################

##tree quality (using the 5 PCoA axes distance matrix)
tree.quality(eu_dist, up_tree)#upgma
tree.quality(eu_dist, nj_tree)#nj
#Note that you get the same results if you use the
#full set of PCoA axes (eu_dist_all)

##hyper quality (from 1 through to 8 axes)
#use the all PCoA axes distance matrix
res_v <- vector(length = 8)
for (i in 1:8){
  TRAITS_v <- trai.pcoa$x[,1:i, drop = FALSE]
  hyper_v <- hyper.build(TRAITS_v, distance = "euclidean", axes = 0)
  res_v[i] <- hyper.quality(eu_dist_all, hyper_v)
}
res_v
#####################################
##plot hyper (Figure 3)
##############################################

#Note, hypervolume construction is stochastic, so the figure
#will not match exactly with the one in the paper
jpeg(file = "hyper.jpeg", width = 20, height = 20,
     units = "cm", res = 300)
plot(hyper, colors = c("gray42", "deepskyblue4"))
dev.off()
