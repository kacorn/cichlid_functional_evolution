################################################################################
# Scripts for data analyses for Martinez, Corn et al. 2024 American Naturalist #
################################################################################
#Statistical analyses and code by Katherine Corn

#load packages
require(tidyverse); require(geomorph); require(ape); require(phytools); 
require(geiger); require(reshape2)

#setwd(#YOUR REPO HERE#)

#assuming you set the working directory to the folder with the data files

Mcgee_2020_tree <- read.tree("Mcgee_2020_tree.tre")

species_table <- read_csv("classifiers.csv") %>%
  rename(Genus_species = genus_species, tree_name = phylo_species, 
         locality = region)

kinematic_components <- 
  read_csv("Avg Kinematic Components (logged _ scaled).csv") %>%
  rename("Genus_species" = `...1`, 
         "premax_protrusion" = "premax protrusion", 
         "maxillary_rotation" = "maxillary rotation", 
         "lower_jaw_rotation" = "lower jaw rotation", 
         "cranial_elevation" = "cranial elevation", 
         "hyoid_depression" = "hyoid depression", 
         "maximum_gape" = "gape")

kinesis_traits <- read_csv("kinesis traits.csv") %>% 
  rename("Genus_species" = `...1`)

all_kinematic_traits <- full_join(species_table, kinematic_components, 
                                  by = "Genus_species") %>%
  full_join(., kinesis_traits, by = "Genus_species")

#import PC axes
start_head_shapes_4pcs <- read_csv("PCs 1-4, START head shapes.csv") %>%
  rename(Genus_species = `...1`, PC1 = PC_1, PC2 = PC_2, PC3 = PC_3, PC4 = PC_4) %>%
  full_join(., species_table, by = "Genus_species") %>%
  dplyr::select(Genus_species, tree_name, locality, PC1, PC2, PC3, PC4)

trajectories_4pcs <- read_csv("PCs 1-4, TRAJECTORY Shape, ALIGNED.csv") %>%
  rename(Genus_species = `...1`, PC1 = PC_1, PC2 = PC_2, PC3 = PC_3, PC4 = PC_4) %>%
  full_join(., species_table, by = "Genus_species") %>%
  dplyr::select(Genus_species, tree_name, locality, PC1, PC2, PC3, PC4)

kinematics_4pcs <- read_csv("PCs 1-6, COMPONENTS.csv") %>%
  rename(Genus_species = `...1`, PC1 = Comp1, PC2 = Comp2, PC3 = Comp3, 
         PC4 = Comp4) %>%
  full_join(., species_table, by = "Genus_species") %>%
  dplyr::select(Genus_species, tree_name, locality, PC1, PC2, PC3, PC4)

#alignment of head shapes with Genus_species rownames
start_head_shapes_alignment <- read.csv("Avg. START head shapes, ALIGNED.csv", 
                                        row.names = 1)

#alignment of trajectories
trajectories_alignment <- read.csv("Avg. Trajectory Shapes, ALIGNED.csv", 
                                   row.names = 1)

#SET UP TREE THAT IS SUBSET OF MCGEE et al 2020 TREE
cichlid_tree <- ladderize(drop.tip(Mcgee_2020_tree, setdiff(Mcgee_2020_tree$tip.label, 
                                              species_table$tree_name)))

#match pruned tree and data to make sure data is in the correct order
species_table <- species_table %>%
  arrange(match(tree_name, cichlid_tree$tip.label))

#checks to make sure everything is in order and we're not missing anything
stopifnot(species_table$tree_name == cichlid_tree$tip.label)

#################################################
# RUN UNIVARIATE DISPARITY FOR KINEMATIC TRAITS #
#################################################


#set up data here
continuous_traits <- all_kinematic_traits %>% 
  dplyr::select(tree_name, premax_protrusion:log_kinesis_coeff_head)
discrete_regimes <- all_kinematic_traits %>% dplyr::select(tree_name, locality)


discrete <- as.vector(t(discrete_regimes[,2]))
names(discrete) <- t(discrete_regimes[,1])
#set up continuous data as a named data_frame
first_title <- names(continuous_traits[,1]) #set of vector to exclude species names column
continuous <- as.data.frame(continuous_traits %>%
                              dplyr::select(everything(), -one_of(first_title)))
rownames(continuous) <- t(continuous_traits[,1])

groups <- levels(as.factor(discrete))
trait_vec <- names(continuous)
#great! our data is all set up and ready to go

#set up data frame to catch disparity results
disparity_results <- data.frame(matrix(ncol = ncol(continuous), nrow = length(groups)))
colnames(disparity_results) <- trait_vec
rownames(disparity_results) <- groups
disparity_results$discrete <- groups

#set up data frame to catch pairwise comparisons significance results. number of pairwise comparisons formula = k*(k-1)/2
num_comparisons <- length(groups)*(length(groups)-1)/2
disparity_sig <- data.frame(matrix(ncol = length(groups) + 2,nrow = 0))
colnames(disparity_sig) <- c("comparison_group", groups, "trait")

sig_res <- data.frame(matrix(ncol = length(groups) + 2, nrow = length(groups)*ncol(continuous)))
colnames(sig_res) <- c("comparison_group", groups, "trait")
sig_res$comparison_group <- rep(groups)

i = 1 #iterates through continuous traits
row = 1 #iterate through rows of pairwise comparisons p-values

#sequential_disparity
stopifnot(!is.na(continuous[,i]))
stopifnot(!is.na(discrete))

for(i in 1:ncol(continuous)){
  
  j = 1 #iterates through groups within disparity analysis
  k = 1 #iterates through pairwise comparisons for significance - ROW INDICATOR
  comp = 1 #iterates rows of pairwise comparisons within each df of p-values
  
  print(paste("Now running disparity for", trait_vec[i]))
  
  run_disparity <- morphol.disparity(f1 = continuous[,i] ~ discrete, groups = ~ discrete, partial = FALSE, 
                                     iter = 10000, seed = NULL, print.progress = F)
  
  sig_results <- tibble(as.data.frame(run_disparity$PV.dist.Pval) %>%
                          rownames_to_column(var = "comparison_group")) %>%
    mutate(trait = rep(trait_vec[i]))
  
  disparity_sig <- rbind(disparity_sig, sig_results)
  
  
  #get pairwise p-values on disparities
  for(comp in 1:length(groups)){
    
    while (k <= length(groups)){
      sig_res[row, k + 1] <- run_disparity$PV.dist.Pval[k,comp]
      sig_res$trait[row] <- trait_vec[i]
      k <- k + 1
    }
    k = 1
    row <- row + 1
    comp <- comp + 1
  }
  
  #get pairwise disparities
  for (j in 1:length(groups)){
    disparity_results[j,i] <- run_disparity$Procrustes.var[j]
    j <- j + 1
  }
  
  i <- i + 1
  
}

#returns 3 sets of results
#1) disparity_results, with disparities
#2) sig_res, which has pairwise p-values
#3) disparity_sig, which also has pairwise p-values but formatted differently in case you want that

#clean up disparity_results a little bit before putting it in the list
clean_disparities <- disparity_results %>%
  rename(locality = discrete) %>%
  dplyr::select(locality, premax_protrusion:log_kinesis_coeff_head)

#make disparity results into a list
kinematic_disparity_results <- list(clean_disparities, disparity_sig)
names(kinematic_disparity_results) <- c("disparities","pairwise_p_values")

#examples of how to extract the table for premax protrusion pairwise p-values in base r
pairwise <- kinematic_disparity_results$pairwise_p_values
pairwise[pairwise$trait == "premax_protrusion", ]



#################################################
# RUN UNIVARIATE RATES FOR KINEMATIC TRAITS     #
#################################################


#set up data here
continuous_traits <- all_kinematic_traits %>% 
  dplyr::select(tree_name, premax_protrusion:log_kinesis_coeff_head)
discrete_regimes <- all_kinematic_traits %>% dplyr::select(tree_name, locality)


discrete <- as.vector(t(discrete_regimes[,2]))
names(discrete) <- t(discrete_regimes[,1])
#set up continuous data as a named data_frame
first_title <- names(continuous_traits[,1]) #set of vector to exclude species names column
continuous <- as.data.frame(continuous_traits %>%
                              dplyr::select(everything(), -one_of(first_title)))
rownames(continuous) <- t(continuous_traits[,1])

groups <- levels(as.factor(discrete))
trait_vec <- names(continuous)
#great! our data is all set up and ready to go

rate_results <- data.frame(matrix(ncol = ncol(continuous), nrow = length(groups)))
colnames(rate_results) <- trait_vec
rownames(rate_results) <- groups
rate_results$discrete <- groups

#set up data frame to catch pairwise comparisons significance results. 
num_comparisons <- length(groups)*(length(groups)-1)/2
rate_sig <- data.frame(matrix(ncol = length(groups) + 2,nrow = 0))
colnames(rate_sig) <- c("comparison_group", groups, "trait")

i = 1 #iterates through continuous traits

#sequential_rates
print("starting to run rates")
for(i in 1:ncol(continuous)){
  
  j = 1 #iterates through groups within rate analysis
  
  cont_trait <- continuous[,i]
  names(cont_trait) <- rownames(continuous)
  
  run_rates <- compare.evol.rates(A = cont_trait, phy = cichlid_tree, gp = discrete, iter = 10000, method = "permutation", print.progress = T)
  
  rate_sig_res <- tibble(melt(as.matrix(run_rates$pairwise.pvalue), varnames = c("locality_1","locality_2"))) %>%
    mutate(trait = rep(trait_vec[i]))
  
  for (j in 1:length(groups)){
    rate_results[j,i] <- run_rates$sigma.d.gp[j]
    j <- j + 1
  }
  
  rate_sig <- rbind(rate_sig, rate_sig_res)
}

#returns 2 sets of results
#1) rate_results, with rates
#2) rate_sig, which has pairwise p-values

#clean up and save results
rate_sig_clean <- as.data.frame(rate_sig %>%
                                  rename(pairwise_p = value))

rate_results_clean <- rate_results %>%
  rename(locality = discrete) %>%
  dplyr::select(locality, premax_protrusion:log_kinesis_coeff_head)

kinematic_rate_results <- list(rate_results_clean, rate_sig_clean)
names(kinematic_rate_results) <- c("rates","pairwise_p_values")


#################################################
# RUN MULTIVARIATE DISPARITIES & RATES          #
#################################################

#RUN MULTIVARIATE KINEMATIC COMPONENTS DISPARITIES

require(reshape2)
#run kinematic rates and disparities
#set up data
kinematic_trait_data <- as.matrix(all_kinematic_traits %>% 
                                    dplyr::select(premax_protrusion:maximum_gape))
rownames(kinematic_trait_data) <- all_kinematic_traits$tree_name

location_regimes <- as.vector(all_kinematic_traits$locality)
names(location_regimes) <- all_kinematic_traits$tree_name

#run
kinematic_rates <- compare.evol.rates(A = kinematic_trait_data, 
                                      phy = cichlid_tree, method = "permutation", 
                                      gp = location_regimes, iter = 10000)

kinematic_disparities <- morphol.disparity(f1 = kinematic_trait_data ~ location_regimes, 
                                           groups = ~ location_regimes, partial = FALSE, 
                                           iter = 10000, seed = NULL, print.progress = TRUE)


#Run head shape rates and disparities
#set up data (match to all_kinematic_traits so we can make rownames match tree names)
matrix_start_shape <- as.matrix(start_head_shapes_alignment %>%
                                  arrange(match(row.names(.), all_kinematic_traits$Genus_species)))
rownames(matrix_start_shape) <- all_kinematic_traits$tree_name

#run
start_shape_rates <- compare.evol.rates(A = matrix_start_shape, phy = cichlid_tree, 
                                        gp = location_regimes, iter = 10000, 
                                        seed = NULL, method = c("permutation"), 
                                        print.progress = TRUE)

start_shape_disparities <- morphol.disparity(f1 = matrix_start_shape ~ location_regimes, 
                                             groups = ~ location_regimes, 
                                             partial = FALSE, iter = 10000, 
                                             seed = NULL, print.progress = TRUE)



#run trajectories analyses
#set up data (match to all_kinematic_traits so we can make rownames match tree names)
matrix_trajectories <- as.matrix(trajectories_alignment %>%
                                   arrange(match(row.names(.), all_kinematic_traits$Genus_species)))
rownames(matrix_trajectories) <- all_kinematic_traits$tree_name

#run
trajectories_rates <- compare.evol.rates(A = matrix_trajectories, phy = cichlid_tree, 
                                         gp = location_regimes, iter = 10000, 
                                         seed = NULL, method = c("permutation"), 
                                         print.progress = TRUE)

trajectories_disparities <- morphol.disparity(f1 = matrix_trajectories ~ location_regimes, 
                                              groups = ~ location_regimes, 
                                              partial = FALSE, iter = 10000, 
                                              seed = NULL, print.progress = TRUE)


#################################################
# COLLATE MULTIVARIATE DISPARITIES              #
#################################################

extract_disp <- function(groups, disparity_res){
  res <- data.frame(row.names = groups)
  i = 1
  for (i in 1:length(groups)){
    res$locality[i] <- groups[i]
    res$disparity[i] <- disparity_res$Procrustes.var[i]
  }
  return(tibble(res))
}

loc <- c("malawi", "new_world", "tanganyika","victoria")

#put together disparity results

kinematic_disp_res <- extract_disp(groups = loc, 
                                   disparity_res = kinematic_disparities) %>%
  rename(kinematic_disparity = disparity)

trajectory_disp_res <- extract_disp(groups = loc, disparity_res = trajectories_disparities) %>%
  rename(trajectory_disparity = disparity)

all_disp_res <- extract_disp(groups = loc, disparity_res = start_shape_disparities) %>%
  rename(head_shape_disparity = disparity) %>%
  full_join(kinematic_disp_res, by = "locality") %>%
  full_join(trajectory_disp_res, by = "locality")


#get disparities pairwise p-values into a neat data frame
kinematic_disparities_pvalues <- melt(kinematic_disparities$PV.dist.Pval) %>% 
  mutate(traits = rep("kinematic_components")) %>%
  rename(locality_1 = Var1, locality_2 = Var2, pairwise_disparity_p = value)

shape_disparities_pvalues <- melt(start_shape_disparities$PV.dist.Pval) %>% 
  mutate(traits = rep("head_shape")) %>%
  rename(locality_1 = Var1, locality_2 = Var2, pairwise_disparity_p = value)

trajectories_disparities_pvalues <- melt(trajectories_disparities$PV.dist.Pval) %>% 
  mutate(traits = rep("trajectories")) %>%
  rename(locality_1 = Var1, locality_2 = Var2, pairwise_disparity_p = value)

all_disparity_pvalues <- bind_rows(kinematic_disparities_pvalues, shape_disparities_pvalues, trajectories_disparities_pvalues)

#integrate it into a list
all_multivariate_disparity_results <- list(all_disp_res, all_disparity_pvalues)
names(all_multivariate_disparity_results) <- c("disparities","pairwise_p_values")

#################################################
# COLLATE MULTIVARIATE RATES                    #
#################################################

#combine all these data

extract_rates <- function(groups, rate_res){
  res <- data.frame(row.names = groups)
  i = 1
  for (i in 1:length(groups)){
    res$locality[i] <- groups[i]
    res$rate[i] <- rate_res$sigma.d.gp[i]
  }
  return(tibble(res))
}

loc <- c("malawi", "new_world", "tanganyika","victoria")

kinematic_results <- extract_rates(groups = loc,rate_res = kinematic_rates) %>%
  rename(kinematic_rate = rate)

trajectory_results <- extract_rates(groups = loc,rate_res = trajectories_rates) %>%
  rename(trajectory_rate = rate)

all_rate_results <- extract_rates(groups = loc, rate_res = start_shape_rates) %>%
  rename(head_shape_rate = rate) %>%
  full_join(kinematic_results, by = "locality") %>%
  full_join(trajectory_results, by = "locality")


#get rate p-values into a single data frame

kinematic_rate_pvalues <- melt(as.matrix(kinematic_rates$pairwise.pvalue)) %>% 
  mutate(traits = rep("kinematic_components")) %>%
  rename(locality_1 = Var1, locality_2 = Var2, pairwise_rate_p = value)

shape_rate_pvalues <- melt(as.matrix(start_shape_rates$pairwise.pvalue)) %>% 
  mutate(traits = rep("head_shape")) %>%
  rename(locality_1 = Var1, locality_2 = Var2, pairwise_rate_p = value)

trajectories_rate_pvalues <- melt(as.matrix(trajectories_rates$pairwise.pvalue)) %>% 
  mutate(traits = rep("trajectories")) %>%
  rename(locality_1 = Var1, locality_2 = Var2, pairwise_rate_p = value)

all_rates_pvalues <- bind_rows(kinematic_rate_pvalues, shape_rate_pvalues, trajectories_rate_pvalues)


#collect overall statistics about multivariate disparities
kinematic_overall_statistics <- c("kinematics", kinematic_rates$sigma.d.ratio, kinematic_rates$P.value, kinematic_rates$Z)
shape_overall_statistics <- c("shape", start_shape_rates$sigma.d.ratio, start_shape_rates$P.value, start_shape_rates$Z)
trajectory_overall_statistics <- c("trajectories", trajectories_rates$sigma.d.ratio, trajectories_rates$P.value, trajectories_rates$Z)

rates_overall_statistics <- as.data.frame(rbind(kinematic_overall_statistics, shape_overall_statistics, trajectory_overall_statistics)) %>%
  rename(traits = V1, overall_sigma2_ratio = V2, overall_pvalue = V3, Z_score = V4)

#integrate it into a list
all_multivariate_rate_results <- list(as.data.frame(all_rate_results), all_rates_pvalues, rates_overall_statistics)
names(all_multivariate_rate_results) <- c("rates","pairwise_p_values", "overall_statistics")

#################################################
# DISPARITY THROUGH TIME ANALYSES               #
#################################################

#################################################
# SET UP DATA MATRICES AND SUBTREES             #
#################################################

#set up df to catch MDI statistics for kinematic data
kinematic_mdi <- as.data.frame(matrix(nrow = 5, ncol = 3))
colnames(kinematic_mdi) <- c("group","MDI","MDI_p_value")

#set up df to catch MDI statistics for shape data
shape_mdi <- as.data.frame(matrix(nrow = 5, ncol = 3))
colnames(shape_mdi) <- c("group","MDI","MDI_p_value")

#set up df to catch MDI statistics for trajectory data
trajectory_mdi <- as.data.frame(matrix(nrow = 5, ncol = 3))
colnames(trajectory_mdi) <- c("group","MDI","MDI_p_value")


#################################################
# WHOLE DATASET DISPARITY THROUGH TIME          #
#################################################

#base data is all_kinematic_traits

#set up the kinematic data
kinematic_trait_matrix <- as.matrix(all_kinematic_traits %>% 
                                      dplyr::select(premax_protrusion:maximum_gape))
rownames(kinematic_trait_matrix) <- all_kinematic_traits$tree_name

#set up the head shape pc data & rename rows to tree_name
allgroups_shape_matrix <- as.matrix(start_head_shapes_4pcs %>%
                                      arrange(match(Genus_species, all_kinematic_traits$Genus_species)) %>%
                                      dplyr::select(PC1:PC4))
rownames(allgroups_shape_matrix) <- all_kinematic_traits$tree_name

#set up the trajectory data & rename rows to tree_name
allgroups_trajectory_matrix <- as.matrix(trajectories_4pcs %>%
                                           arrange(match(Genus_species, all_kinematic_traits$Genus_species)) %>%
                                           dplyr::select(PC1:PC4))
rownames(allgroups_trajectory_matrix) <- all_kinematic_traits$tree_name

#tree is base tree

kinematic_allgroups_dtt <- dtt(phy = cichlid_tree, data = kinematic_trait_matrix, index = "avg.sq",
                               mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                               nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

kinematic_mdi[1,] <- c("all_groups", kinematic_allgroups_dtt$MDI, kinematic_allgroups_dtt$MDIpVal)


#head shape data
allgroups_shape_dtt <- dtt(phy = cichlid_tree, data = allgroups_shape_matrix, index = "avg.sq",
                           mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                           nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)
shape_mdi[1,] <- c("allgroups", allgroups_shape_dtt$MDI, allgroups_shape_dtt$MDIpVal)


#trajectory data
allgroups_trajectory_dtt <- dtt(phy = cichlid_tree, data = allgroups_trajectory_matrix, index = "avg.sq",
                                mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                                nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

trajectory_mdi[1,] <- c("allgroups", allgroups_trajectory_dtt$MDI, allgroups_trajectory_dtt$MDIpVal)



#################################################
# TANGANYIKA DISPARITY THROUGH TIME             #
#################################################

#Set up reference dataset
tanganyika_kinematic_data <- all_kinematic_traits %>% filter(locality == "tanganyika")

#set up the kinematic data
tanganyika_kinematic_matrix <- as.matrix(tanganyika_kinematic_data %>% 
                                           dplyr::select(premax_protrusion:maximum_gape))
row.names(tanganyika_kinematic_matrix) <- tanganyika_kinematic_data$tree_name

#set up the head shape pc data & rename rows to tree_name
tanganyika_shape_matrix <- as.matrix(start_head_shapes_4pcs %>%
                                       filter(Genus_species %in% tanganyika_kinematic_data$Genus_species) %>%
                                       arrange(match(Genus_species, tanganyika_kinematic_data$Genus_species)) %>%
                                       dplyr::select(PC1:PC4))
rownames(tanganyika_shape_matrix) <- tanganyika_kinematic_data$tree_name

#set up the trajectory data & rename rows to tree_name
tanganyika_trajectory_matrix <- as.matrix(trajectories_4pcs %>%
                                            filter(Genus_species %in% tanganyika_kinematic_data$Genus_species) %>%
                                            arrange(match(Genus_species, tanganyika_kinematic_data$Genus_species)) %>%
                                            dplyr::select(PC1:PC4))
rownames(tanganyika_trajectory_matrix) <- tanganyika_kinematic_data$tree_name

#set up the tree
tanganyika_tree <- drop.tip(phy = cichlid_tree, setdiff(cichlid_tree$tip.label, tanganyika_kinematic_data$tree_name))

#run disparity through time
#kinematic data
tanganyika_kinematic_dtt <- dtt(phy = tanganyika_tree, data = tanganyika_kinematic_matrix, index = "avg.sq",
                                mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                                nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

kinematic_mdi[2,] <- c("tanganyika", tanganyika_kinematic_dtt$MDI, tanganyika_kinematic_dtt$MDIpVal)


#head shape data
tanganyika_shape_dtt <- dtt(phy = tanganyika_tree, data = tanganyika_shape_matrix, index = "avg.sq",
                            mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                            nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

shape_mdi[2,] <- c("tanganyika", tanganyika_shape_dtt$MDI, tanganyika_shape_dtt$MDIpVal)

#trajectory data
tanganyika_trajectory_dtt <- dtt(phy = tanganyika_tree, data = tanganyika_trajectory_matrix, index = "avg.sq",
                                 mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                                 nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

trajectory_mdi[2,] <- c("tanganyika", tanganyika_trajectory_dtt$MDI, tanganyika_trajectory_dtt$MDIpVal)


#################################################
# MALAWI DISPARITY THROUGH TIME                 #
#################################################

#Set up reference dataset
malawi_kinematic_data <- all_kinematic_traits %>% filter(locality == "malawi")

#set up the kinematic data
malawi_kinematic_matrix <- as.matrix(malawi_kinematic_data %>% 
                                       dplyr::select(premax_protrusion:maximum_gape))
row.names(malawi_kinematic_matrix) <- malawi_kinematic_data$tree_name

#set up the head shape pc data & rename rows to tree_name
malawi_shape_matrix <- as.matrix(start_head_shapes_4pcs %>%
                                   filter(Genus_species %in% malawi_kinematic_data$Genus_species) %>%
                                   arrange(match(Genus_species, malawi_kinematic_data$Genus_species)) %>%
                                   dplyr::select(PC1:PC4))
rownames(malawi_shape_matrix) <- malawi_kinematic_data$tree_name

#set up the trajectory data & rename rows to tree_name
malawi_trajectory_matrix <- as.matrix(trajectories_4pcs %>%
                                        filter(Genus_species %in% malawi_kinematic_data$Genus_species) %>%
                                        arrange(match(Genus_species, malawi_kinematic_data$Genus_species)) %>%
                                        dplyr::select(PC1:PC4))
rownames(malawi_trajectory_matrix) <- malawi_kinematic_data$tree_name

#set up the tree
malawi_tree <- drop.tip(phy = cichlid_tree, setdiff(cichlid_tree$tip.label, malawi_kinematic_data$tree_name))

#run disparity through time
#kinematic data
malawi_kinematic_dtt <- dtt(phy = malawi_tree, data = malawi_kinematic_matrix, index = "avg.sq",
                            mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                            nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

kinematic_mdi[3,] <- c("malawi", malawi_kinematic_dtt$MDI, malawi_kinematic_dtt$MDIpVal)


#head shape data
malawi_shape_dtt <- dtt(phy = malawi_tree, data = malawi_shape_matrix, index = "avg.sq",
                        mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                        nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

shape_mdi[3,] <- c("malawi", malawi_shape_dtt$MDI, malawi_shape_dtt$MDIpVal)

#trajectory data
malawi_trajectory_dtt <- dtt(phy = malawi_tree, data = malawi_trajectory_matrix, index = "avg.sq",
                             mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                             nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

trajectory_mdi[3,] <- c("malawi", malawi_trajectory_dtt$MDI, malawi_trajectory_dtt$MDIpVal)


#################################################
# VICTORIA DISPARITY THROUGH TIME               #
#################################################

#Set up reference dataset
victoria_kinematic_data <- all_kinematic_traits %>% filter(locality == "victoria")

#set up the kinematic data
victoria_kinematic_matrix <- as.matrix(victoria_kinematic_data %>% 
                                         dplyr::select(premax_protrusion:maximum_gape))
row.names(victoria_kinematic_matrix) <- victoria_kinematic_data$tree_name

#set up the head shape pc data & rename rows to tree_name
victoria_shape_matrix <- as.matrix(start_head_shapes_4pcs %>%
                                     filter(Genus_species %in% victoria_kinematic_data$Genus_species) %>%
                                     arrange(match(Genus_species, victoria_kinematic_data$Genus_species)) %>%
                                     dplyr::select(PC1:PC4))
rownames(victoria_shape_matrix) <- victoria_kinematic_data$tree_name

#set up the trajectory data & rename rows to tree_name
victoria_trajectory_matrix <- as.matrix(trajectories_4pcs %>%
                                          filter(Genus_species %in% victoria_kinematic_data$Genus_species) %>%
                                          arrange(match(Genus_species, victoria_kinematic_data$Genus_species)) %>%
                                          dplyr::select(PC1:PC4))
rownames(victoria_trajectory_matrix) <- victoria_kinematic_data$tree_name

#set up the tree
victoria_tree <- drop.tip(phy = cichlid_tree, setdiff(cichlid_tree$tip.label, victoria_kinematic_data$tree_name))

#run disparity through time
#kinematic data
victoria_kinematic_dtt <- dtt(phy = victoria_tree, data = victoria_kinematic_matrix, index = "avg.sq",
                              mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                              nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

kinematic_mdi[4,] <- c("victoria", victoria_kinematic_dtt$MDI, victoria_kinematic_dtt$MDIpVal)


#head shape data
victoria_shape_dtt <- dtt(phy = victoria_tree, data = victoria_shape_matrix, index = "avg.sq",
                          mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                          nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

shape_mdi[4,] <- c("victoria", victoria_shape_dtt$MDI, victoria_shape_dtt$MDIpVal)

#trajectory data
victoria_trajectory_dtt <- dtt(phy = victoria_tree, data = victoria_trajectory_matrix, index = "avg.sq",
                               mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                               nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

trajectory_mdi[4,] <- c("victoria", victoria_trajectory_dtt$MDI, victoria_trajectory_dtt$MDIpVal)


#################################################
# NEW WORLD DISPARITY THROUGH TIME              #
#################################################

#Set up reference dataset
new_world_kinematic_data <- all_kinematic_traits %>% filter(locality == "new_world")

#set up the kinematic data
new_world_kinematic_matrix <- as.matrix(new_world_kinematic_data %>% 
                                          dplyr::select(premax_protrusion:maximum_gape))
row.names(new_world_kinematic_matrix) <- new_world_kinematic_data$tree_name

#set up the head shape pc data & rename rows to tree_name
new_world_shape_matrix <- as.matrix(start_head_shapes_4pcs %>%
                                      filter(Genus_species %in% new_world_kinematic_data$Genus_species) %>%
                                      arrange(match(Genus_species, new_world_kinematic_data$Genus_species)) %>%
                                      dplyr::select(PC1:PC4))
rownames(new_world_shape_matrix) <- new_world_kinematic_data$tree_name

#set up the trajectory data & rename rows to tree_name
new_world_trajectory_matrix <- as.matrix(trajectories_4pcs %>%
                                           filter(Genus_species %in% new_world_kinematic_data$Genus_species) %>%
                                           arrange(match(Genus_species, new_world_kinematic_data$Genus_species)) %>%
                                           dplyr::select(PC1:PC4))
rownames(new_world_trajectory_matrix) <- new_world_kinematic_data$tree_name

#set up the tree
new_world_tree <- drop.tip(phy = cichlid_tree, setdiff(cichlid_tree$tip.label, new_world_kinematic_data$tree_name))

#run disparity through time
#kinematic data
new_world_kinematic_dtt <- dtt(phy = new_world_tree, data = new_world_kinematic_matrix, index = "avg.sq",
                               mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                               nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

kinematic_mdi[5,] <- c("new_world", new_world_kinematic_dtt$MDI, new_world_kinematic_dtt$MDIpVal)


#head shape data
new_world_shape_dtt <- dtt(phy = new_world_tree, data = new_world_shape_matrix, index = "avg.sq",
                           mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                           nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

shape_mdi[5,] <- c("new_world", new_world_shape_dtt$MDI, new_world_shape_dtt$MDIpVal)

#trajectory data
new_world_trajectory_dtt <- dtt(phy = new_world_tree, data = new_world_trajectory_matrix, index = "avg.sq",
                                mdi.range = c(0, 0.75), #hash this line if you want to run the MDI calculations on the whole tree
                                nsim = 10000, CI = 0.95, plot = T, calculateMDIp = T)

trajectory_mdi[5,] <- c("new_world", new_world_trajectory_dtt$MDI, new_world_trajectory_dtt$MDIpVal)


#################################################
# COLLATE MDI SCORES                            #
#################################################

all_MDI_table <- list(kinematic_mdi, shape_mdi, trajectory_mdi)
names(all_MDI_table) <- c("kinematic_MDI", "shape_MDI", "trajectory_MDI")

# <- end disparity through time ->


###########################################################
# HYPERVOLUMES - iterate the following for each data type #
###########################################################

##########################################
# SET UP THE REAL HYPERVOLUMES           #
##########################################

require(hypervolume)

#TANGANYIKA
tanganyika_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality == "tanganyika") %>%
                                               dplyr::select(PC1:PC4)),
                                      name = "tanganyika")

not_tanganyika_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality != "tanganyika") %>%
                                                   dplyr::select(PC1:PC4)),
                                          name = "not_tanganyika")


#NEW WORLD
new_world_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality == "new_world") %>%
                                              dplyr::select(PC1:PC4)),
                                     name = "new_world")

not_new_world_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality != "new_world") %>%
                                                  dplyr::select(PC1:PC4)),
                                         name = "not_new_world")


#MALAWI
malawi_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality == "malawi") %>%
                                           dplyr::select(PC1:PC4)),
                                  name = "malawi")

not_malawi_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality != "malawi") %>%
                                               dplyr::select(PC1:PC4)),
                                      name = "not_malawi")


#VICTORIA
victoria_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality == "victoria") %>%
                                             dplyr::select(PC1:PC4)),
                                    name = "victoria")

not_victoria_hv <- hypervolume_gaussian(data= (start_head_shapes_4pcs %>% filter(locality != "victoria") %>%
                                                 dplyr::select(PC1:PC4)),
                                        name = "not_victoria")



##########################################
# SET UP SIMULATED HYPERVOLUMES          #
##########################################

#Note April 2024.
#You do not need to do permutations of hypervolumes to assess significance with 
  #this script with hypervolume v3.0+
#This script was made before the hypervolume package had its own permutation fuctions
    #You can run permuted hypervolumes natively in the package now.

require(foreach)
require(parallel)
require(doParallel)

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

run_hvs_par <- function(head_shape_pca_data_hvs){
  
  #set up data frame to catch results
  hv_results <- data.frame(matrix(nrow = 10, ncol = 11))
  colnames(hv_results) <- c("hv1","hv2","hv1_centroid","hv2_centroid","hv1_volume","hv2_volume","dist_centr",
                            "overlap_jaccard","overlap_sorensen","frac_unique_1","frac_unique_2")
  
  #set up vectors to pull comparisons
  hv_1_vector_names <- c("tanganyika_hv", "new_world_hv", "malawi_hv","victoria_hv", 
                         "tanganyika_hv", "tanganyika_hv", "tanganyika_hv","malawi_hv", "malawi_hv", "victoria_hv")
  
  hv_2_vector_names <- c("not_tanganyika_hv", "not_new_world_hv","not_malawi_hv", "not_victoria_hv", "new_world_hv", "malawi_hv", "victoria_hv", "victoria_hv", "new_world_hv", "new_world_hv")
  
  r = 1
  
  permuted_samples <- head_shape_pca_data_hvs %>%
    dplyr::select(PC1:PC4) 
  permuted_samples$permuted_locality <- head_shape_pca_data_hvs[,j+6]
  
  
  #TANGANYIKA
  tanganyika_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality == "tanganyika") %>%
                                                 dplyr::select(PC1:PC4)),
                                        name = "tanganyika", verbose = FALSE)
  
  not_tanganyika_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality != "tanganyika") %>%
                                                     dplyr::select(PC1:PC4)),
                                            name = "not_tanganyika", verbose = FALSE)
  #MALAWI
  malawi_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality == "malawi") %>%
                                             dplyr::select(PC1:PC4)),
                                    name = "malawi", verbose = FALSE)
  
  not_malawi_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality != "malawi") %>%
                                                 dplyr::select(PC1:PC4)),
                                        name = "not_malawi", verbose = FALSE)
  #VICTORIA
  victoria_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality == "victoria") %>%
                                               dplyr::select(PC1:PC4)),
                                      name = "victoria", verbose = FALSE)
  
  not_victoria_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality != "victoria") %>%
                                                   dplyr::select(PC1:PC4)),
                                          name = "not_victoria", verbose = FALSE)
  #NEW WORLD
  new_world_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality == "new_world") %>%
                                                dplyr::select(PC1:PC4)),
                                       name = "new_world", verbose = FALSE)
  
  not_new_world_hv <- hypervolume_gaussian(data= (permuted_samples %>% filter(permuted_locality != "new_world") %>%
                                                    dplyr::select(PC1:PC4)),
                                           name = "not_new_world", verbose = FALSE)
  
  #add new hvs into vectors
  hv_1_vector <- c(tanganyika_hv, new_world_hv, malawi_hv, victoria_hv, tanganyika_hv, tanganyika_hv, 
                   tanganyika_hv, malawi_hv, malawi_hv, victoria_hv)
  
  hv_2_vector <- c(not_tanganyika_hv, not_new_world_hv,not_malawi_hv,not_victoria_hv, new_world_hv, 
                   malawi_hv, victoria_hv, victoria_hv, new_world_hv, new_world_hv)
  
  #ok we've generated our new hypervolumes. now we're going to run through our other loop to create all the comparisons
  i = 1
  
  
  for(i in 1:10){
    
    print(paste("comparing",hv_1_vector_names[i],"and",hv_2_vector_names[i]))
    print(paste("this is comparison",i,"of 10"))
    
    hv1 <- hv_1_vector[[i]]
    hv2 <- hv_2_vector[[i]]
    
    dist_min <- hypervolume_distance(hv1, hv2, type = "minimum", num.points.max = 10000, check.memory = F)
    dist_centr <- hypervolume_distance(hv1, hv2, type = "centroid", num.points.max = 10000, check.memory = F)
    
    hvset <- hypervolume_set(hv1, hv2, num.points.max = 100000, check.memory = F, verbose = FALSE)
    
    hv_overlap <- hypervolume_overlap_statistics(hvset)
    
    hv_results$hv1[r] <- hv_1_vector_names[i]
    hv_results$hv2[r] <- hv_2_vector_names[i]
    hv_results$hv1_centroid[r] <- list(get_centroid(hv1))
    hv_results$hv2_centroid[r] <- list(get_centroid(hv2))
    hv_results$hv1_volume[r] <- get_volume(hv1)
    hv_results$hv2_volume[r] <- get_volume(hv2)
    hv_results$dist_centr[r] <- dist_centr
    hv_results$overlap_jaccard[r] <- hv_overlap[[1]] # Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)
    hv_results$overlap_sorensen[r] <- hv_overlap[[2]] # Sorensen similarity (twice the volume of intersection of 1 and 2 divided by vol- ume of 1 plus volume of 2)
    hv_results$frac_unique_1[r] <- hv_overlap[[3]] #Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
    hv_results$frac_unique_2[r] <- hv_overlap[[4]] #Unique fraction 2 (volume of unique component of 2 divided by volume of 2))
    
    i = i + 1
    r = r + 1
    
  }
  
  #index = index + 1
  
  #            }
  return(hv_results)
  
}

head_shape_pca_data_hvs <- start_head_shapes_4pcs %>% dplyr::select(Genus_species,locality,PC1:PC4)

permutation = 1
for(permutation in 1:10000){
  head_shape_pca_data_hvs[,permutation+6] <- sample(head_shape_pca_data_hvs$locality)
}

j = 1

parallel_hvs <- foreach(j = 1:10000, combine = data.frame) %dopar% {
  run_hvs_par(head_shape_pca_data_hvs = head_shape_pca_data_hvs)
}

merge_hvs <- parallel_hvs %>%
  bind_rows(., .id = "column_label") %>%
  rename(permutation_count = column_label)

parallel::stopCluster(cl = my.cluster)

##########################################
# COMPARE HYPERVOLUMES TO PERMUTATIONS   #
##########################################

library(readr)

hv_1_vector <- c(tanganyika_hv,
                 new_world_hv,
                 malawi_hv,
                 victoria_hv,
                 tanganyika_hv,
                 tanganyika_hv,
                 tanganyika_hv,
                 malawi_hv,
                 malawi_hv,
                 victoria_hv)

hv_2_vector <- c(not_tanganyika_hv,
                 not_new_world_hv,
                 not_malawi_hv,
                 not_victoria_hv,
                 new_world_hv,
                 malawi_hv,
                 victoria_hv,
                 victoria_hv,
                 new_world_hv,
                 new_world_hv)

hv_1_vector_names <- c("tanganyika_hv",
                       "new_world_hv",
                       "malawi_hv",
                       "victoria_hv",
                       "tanganyika_hv",
                       "tanganyika_hv",
                       "tanganyika_hv",
                       "malawi_hv",
                       "malawi_hv",
                       "victoria_hv")

hv_2_vector_names <- c("not_tanganyika_hv",
                       "not_new_world_hv",
                       "not_malawi_hv",
                       "not_victoria_hv",
                       "new_world_hv",
                       "malawi_hv",
                       "victoria_hv",
                       "victoria_hv",
                       "new_world_hv",
                       "new_world_hv")

hv_results <- tibble(id = seq(1:10))
hv_results$hv1 <- rep(NA, 10)
hv_results$hv2 <- rep(NA, 10)
hv_results$hv1_centroid <- rep(NA, 10)
hv_results$hv2_centroid <- rep(NA, 10)
hv_results$hv1_volume <- rep(NA, 10)
hv_results$hv2_volume <- rep(NA, 10)
hv_results$dist_centr <- rep(NA, 10)
hv_results$dist_min <- rep(NA, 10)
hv_results$overlap_jaccard <- rep(NA, 10) # Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)
hv_results$overlap_sorensen <- rep(NA, 10) # Sorensen similarity (twice the volume of intersection of 1 and 2 divided by vol- ume of 1 plus volume of 2)
hv_results$frac_unique_1 <- rep(NA, 10) #Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
hv_results$frac_unique_2 <- rep(NA, 10) #Unique fraction 2 (volume of unique component of 2 divided by volume of 2))


i = 1
r = 1


for(i in 1:10){
  hv1 <- hv_1_vector[[i]]
  hv2 <- hv_2_vector[[i]]
  
  dist_min <- hypervolume_distance(hv1, hv2, type = "minimum", num.points.max = 10000, check.memory = F)
  dist_centr <- hypervolume_distance(hv1, hv2, type = "centroid", num.points.max = 10000, check.memory = F)
  
  hvset <- hypervolume_set(hv1, hv2, num.points.max = 100000, check.memory = F)
  
  hv_overlap <- hypervolume_overlap_statistics(hvset)
  
  hv_results$hv1[r] <- hv_1_vector_names[i]
  hv_results$hv2[r] <- hv_2_vector_names[i]
  hv_results$hv1_centroid[r] <- list(get_centroid(hv1))
  hv_results$hv2_centroid[r] <- list(get_centroid(hv2))
  hv_results$hv1_volume[r] <- get_volume(hv1)
  hv_results$hv2_volume[r] <- get_volume(hv2)
  hv_results$dist_centr[r] <- dist_centr
  hv_results$dist_min[r] <- dist_min
  hv_results$overlap_jaccard[r] <- hv_overlap[[1]] # Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)
  hv_results$overlap_sorensen[r] <- hv_overlap[[2]] # Sorensen similarity (twice the volume of intersection of 1 and 2 divided by vol- ume of 1 plus volume of 2)
  hv_results$frac_unique_1[r] <- hv_overlap[[3]] #Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
  hv_results$frac_unique_2[r] <- hv_overlap[[4]] #Unique fraction 2 (volume of unique component of 2 divided by volume of 2))

  i = i + 1
  r = r + 1
  
}

real_shape_hvs <- hv_results  %>%
  mutate(permutation_count = rep("real")) %>%
  dplyr::select(everything(), -id, -dist_min)  %>%
  mutate(permuted = rep(FALSE))

hv_results_focal <- real_shape_hvs #update this to real dataset of variable of interest
hv_results_permuted <- merge_hvs #update this to simulations of variable of interest 

all_hv_results <- bind_rows(hv_results_focal, hv_results_permuted)



hv_1_vector_names <- c("tanganyika_hv", "new_world_hv", "malawi_hv","victoria_hv", "tanganyika_hv", 
                       "tanganyika_hv", "tanganyika_hv","malawi_hv", "malawi_hv", "victoria_hv")

hv_2_vector_names <- c("not_tanganyika_hv", "not_new_world_hv","not_malawi_hv", "not_victoria_hv",
                       "new_world_hv", "malawi_hv", "victoria_hv", "victoria_hv", "new_world_hv", "new_world_hv")

list_of_test_statistics <- c("overlap_sorensen","frac_unique_1","frac_unique_2")

p_table <- data.frame(matrix(data = NA, nrow = length(hv_1_vector_names), ncol = (length(list_of_test_statistics)+2)))
colnames(p_table) <- c("hv1", "hv2", list_of_test_statistics)

comparison = 1

for (comparison in 1:length(hv_1_vector_names)) {
  
  hv_1 = hv_1_vector_names[comparison]
  hv_2 = hv_2_vector_names[comparison]
  
  hv_data <- all_hv_results %>% filter(hv1 == hv_1 & hv2 == hv_2)
  unpermuted_values <- hv_data %>% filter(permuted == F)
  permuted_data <- hv_data %>% filter(permuted == T)
  
  p_table[comparison,"hv1"] <- hv_1
  p_table[comparison,"hv2"] <- hv_2
  
  stat = 1
  
  for (stat in 1:length(list_of_test_statistics)){
    
    test_statistic <- list_of_test_statistics[stat]
    unpermuted_value <- unlist(unpermuted_values[1,test_statistic])
    
    if(test_statistic == "overlap_jaccard" | test_statistic == "overlap_sorensen"){
      
      #prop_extreme <- 
      num_less_extreme <- table((permuted_data[,test_statistic]) < unpermuted_value)[[1]]
      
      prop_extreme <- num_less_extreme/nrow(hv_data)
      rounded_prop_extreme <- floor(prop_extreme*1000)/1000
      
    } else {
      
      num_less_extreme <- table((permuted_data[,test_statistic]) > unpermuted_value)[[1]]
      
      prop_extreme <- (num_less_extreme/nrow(hv_data))
      rounded_prop_extreme <- floor(prop_extreme*1000)/1000
      
    }
    
    p_table[comparison,test_statistic] <- rounded_prop_extreme
    
    stat <- stat + 1
    
  }
  
  comparison <- comparison + 1
  
}



