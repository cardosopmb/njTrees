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
