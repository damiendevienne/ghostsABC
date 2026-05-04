## HERE IS THE TRUE ABC PROCESS. 

source("functions4abc.R")

args <- commandArgs(trailingOnly = TRUE)

d <- as.numeric(args[1])
t <- as.numeric(args[2])
l <- as.numeric(args[3])
rho <- as.numeric(args[4])

print(paste(d,t,l,rho, sep=" - "))
# start from an observed gene tree and an observed set of gene trees. 
# WE DO A SELF-TEST, so to get these trees for now we do simulations with known parameters
N_extant <- 25L #L is to force it to be an integer (not a double)
target_rho <- rho # proportion of extant species that are sampled (i.e. that are tips of the observed species tree)
#dtl <- sample_theta()
dtl <- c(d=d,t=t,l=l)
N_gene_trees <- 40L # Number of observed gene trees

# SIMULATE SPECIES TREE
birth <- 1.0
death <- 0
sp_tree <- simulate_species_tree(N_extant, birth, death) ## -> PASS BY RPHYLO MAYBE? 
sp_tree_ape <- read.tree(text=tree_to_newick(sp_tree))       # ape format
sp_tree_ape$root.edge <- 0 #no branch before root to ensure mrca(aug.tree) = mrca(tree)
sp_tree_ape$tip.label <- paste0("e", sp_tree_ape$tip.label)
#plot(sp_tree_ape)
# We augment the tree to add ghosts (here for now only extant but unsampled species)
augmented_sp_tree_ape <- augment_tree(sp_tree_ape, 4,0,target_rho, 1,1)
augmented_sp_tree_rustree <- parse_newick(write.tree(augmented_sp_tree_ape))
true_rho <- N_extant / Ntip(augmented_sp_tree_ape)
observed_gene_trees <- sim_gene_trees_fun(augmented_sp_tree_rustree, N_gene_trees, dtl["d"], dtl["t"], dtl["l"])
obs_SCORES <- do.call(rbind, lapply(observed_gene_trees, compute_informative_scores))






smc_abc_eval <- function(THETAS, sp_tree_ape, N_gene_trees, obs_SCORES) {
  p <- 11 #11 scores
  print("Augmenting trees...")
  multi_augmented_sp_trees_ape <- augment_tree_batch(
    sp_tree_ape, 4, 0,
    THETAS[, "rho"], 1, 1
  )
  print("Augmenting trees... DONE")
  ks_matrix <- matrix(NA_real_, nrow(THETAS), p)
  for (i in seq_len(nrow(THETAS))) {
    theta <- THETAS[i, ]
    if (i %% 50 == 0) {
      gc()
      cat(paste(i,paste(theta,collapse="-"), sep=" -> "), "\n")
    }

    sp_tree_aug <- parse_newick(write.tree(
      multi_augmented_sp_trees_ape$trees[[i]]
    ))
    simulated_trees <- sim_gene_trees_fun(
      sp_tree_aug,
      N_gene_trees,
      theta[1], theta[2], theta[3]
    )
    scores <- t(vapply(simulated_trees,compute_informative_scores,numeric(p)))
    ks_matrix[i, ] <- compute_ks_vector(obs_SCORES, scores)
    rm(simulated_trees, scores)
  }
  list(
    ks_matrix = ks_matrix,
    d_abc = rowMeans(ks_matrix)
  )
}


get_new_thetas <- function(THET, dabc, epsilon=0.05, N_keep=40, N_simulations=400) {
  THETAS_sel <- THET[order(dabc)[1:N_keep],] #approximate posterior sample
  THETAS_new <- NULL
  for (i in seq_len(N_simulations)) {
    base <- THETAS_sel[sample(nrow(THETAS_sel), 1), ] # sample one retained particle
    noise <- runif(length(base), -epsilon, epsilon) # add random perturbation in [-epsilon, epsilon]
    theta <- base + noise
    theta[theta < 0] <- -theta[theta < 0] #reflective boundary at 0
    theta[theta > 1] <- 2 - theta[theta > 1] #reflective boundary at 1
    #theta <- pmin(pmax(theta, 0), 1) # prevent negative values and values > 1 
    THETAS_new <- rbind(THETAS_new, theta)
  }
  THETAS_new
}

N_simulations <- 400
# ABC STEP 1
THETAS_init <- sample_theta(N_simulations) # Sample DTL values
THETAS_init$rho <- runif(N_simulations, 0, 1) # Sample RHO values#
res1 <- smc_abc_eval(THETAS_init, sp_tree_ape, N_gene_trees, obs_SCORES)
# ABC STEP 2
THETAS1 <- get_new_thetas(THETAS_init, res1$d_abc, epsilon=0.15)
{
  res2 <- smc_abc_eval(THETAS1, sp_tree_ape, N_gene_trees, obs_SCORES)
# ABC STEP 3
THETAS2 <- get_new_thetas(THETAS1, res2$d_abc, epsilon=0.15)
res3 <- smc_abc_eval(THETAS2, sp_tree_ape, N_gene_trees, obs_SCORES)
# ABC STEP 4
THETAS3 <- get_new_thetas(THETAS2, res3$d_abc, epsilon=0.1)
res4 <- smc_abc_eval(THETAS3, sp_tree_ape, N_gene_trees, obs_SCORES)
# ABC STEP 5
THETAS4 <- get_new_thetas(THETAS3, res4$d_abc, epsilon=0.1)
res5 <- smc_abc_eval(THETAS4, sp_tree_ape, N_gene_trees, obs_SCORES)
# ABC STEP 6
THETAS5 <- get_new_thetas(THETAS4, res5$d_abc, epsilon=0.05)
res6 <- smc_abc_eval(THETAS5, sp_tree_ape, N_gene_trees, obs_SCORES)
# ABC STEP 7
THETAS6 <- get_new_thetas(THETAS5, res6$d_abc, epsilon=0.05)
res7 <- smc_abc_eval(THETAS6, sp_tree_ape, N_gene_trees, obs_SCORES)
THETAS7 <- get_new_thetas(THETAS6, res7$d_abc, epsilon=0.05)
}




######################
## PLOT             ##
######################

theta_list <- list(
  THETAS_init,
  THETAS1,
  THETAS2,
  THETAS3,
  THETAS4,
  THETAS5,
  THETAS6, 
  THETAS7
)

d_list <- list(
  res1$d_abc,
  res2$d_abc,
  res3$d_abc,
  res4$d_abc,
  res5$d_abc,
  res6$d_abc,
  res7$d_abc
)

pdf(sprintf(
  "results7/res_d%.3f_t%.3f_l%.3f_rho%.3f.pdf",
  d,t,l,rho
), width=8, height=8)


par(mfrow = c(2,2))

params <- c("d","t","l","rho")
true_vals <- c(dtl[1], dtl[2], dtl[3], true_rho)

for (j in 1:4) {

  boxplot(
    lapply(theta_list, function(x) x[, j]),
    main = params[j],
    xlab = "Iteration",
    ylab = "Value"
  )

  abline(h = true_vals[j], col = "red", lty = 2, lwd = 2)
}

dev.off()

save(theta_list,d_list, dtl,true_rho,
  file = sprintf(
    "results7/res_d%.3f_t%.3f_l%.3f_rho%.3f.Rdata",
    d,t,l,rho
  )
)

