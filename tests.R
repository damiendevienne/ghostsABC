# let's see the effect of different parameters on the distribution of summary statistics. 

source("functions4abc.R")




N_extant <- 25L #L is to force it to be an integer (not a double)
birth <- 1.0
sp_tree <- simulate_species_tree(N_extant, birth, death, 33L)
sp_tree_ape <- read.tree(text=tree_to_newick(sp_tree))       # ape format
plot(sp_tree_ape)

# Paramètres
d_max <- 1
t_max <- 1
l_max <- 1
death_max <- 1
# Tirage aléatoire de theta
sample_theta <- function() {
  c(
    d = runif(1, 0, d_max),
    t = runif(1, 0, t_max),
    l = runif(1, 0, l_max)
  )
}


meanSCORES  <- NULL
THETAS <- NULL
for (i in 1:100) {
  print(i)
  theta <- sample_theta()
  THETAS <- rbind(THETAS, theta)
  gene_trees_obs <- simulate_dtl_batch( #_obs is for "observed"
    species_tree = sp_tree,
    n = 20L,          # Number of gene trees
    lambda_d = theta["d"],
    lambda_t = theta["t"],
    lambda_l = theta["l"],
    transfer_alpha = 0,
    require_extant = TRUE,
    seed = 123L
  )
  gene_trees_obs_ape_clean <- lapply(gene_trees_obs, keep_extant_tips_gn_full)
  # remove gene trees with less than 4 tips (some summary statistics cannot be computed on such small trees)
  gene_trees_obs_ape_clean <- gene_trees_obs_ape_clean[lapply(gene_trees_obs_ape_clean, Ntip) >= 4]
  SCORES <- do.call(rbind, lapply(gene_trees_obs_ape_clean, compute_all_scores2))
  meanSCORES <- rbind(meanSCORES, apply(SCORES, 2, mean))
}

colnames(THETAS) <- c("d", "t", "l")

meanSCORES_norm <- as.matrix(
  apply(meanSCORES, 2, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0, length(x)))
    (x - min(x)) / rng
  })
)

rownames(meanSCORES_norm) <- rownames(meanSCORES)


cor_mat <- sapply(1:ncol(THETAS), function(i) {
  apply(meanSCORES_norm, 2, function(x) {
    cor(THETAS[, i], x, method = "spearman", use = "pairwise.complete.obs")
  })
})

rownames(cor_mat) <- colnames(meanSCORES_norm)
colnames(cor_mat) <- colnames(THETAS)




#gene_trees_obs_ape <- lapply(gene_trees_obs, function(gt) read.tree(text = tree_to_newick(gt))) # ape format

# CLEAN THE GENE TREES TO ONLY KEEP VISIBLE (extant) TIPS (i.e. NO GHOSTS)



# -----------------------------
# 2. PERFORM SIMULATIONS AND COMPARE SIMULATD AND 
#   "observed" datasets
# -----------------------------

