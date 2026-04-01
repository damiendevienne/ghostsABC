### GHERE WE WILL DO ABC approach to recover DTL rates (first). 
### WHEN WORKING, WE CAN THEN TRY TO RECOVER THE PRESENCE OF GHOSTS (i.e. EXTINCT/UNSAMPLED LINEAGES) IN THE DATASET. 


#SOURCE ALL USEFUL FUNCTIONS
source("functions4abc.R")

# -----------------------------
# 1. simulate the "observed" dataset
#  -> one species tree and 100 gene trees
#  -> no ghosts, i.e. no extinct/unsampled lineages
#  -> Similar to what Baudet et al. refer to as "self-test" 
# -----------------------------

N_extant <- 10L #L is to force it to be an integer (not a double)
birth <- 1.0
death <- 0

# Simulate a species tree
sp_tree <- simulate_species_tree(N_extant, birth, death, 42L)
sp_tree_ape <- read.tree(text=tree_to_newick(sp_tree))       # ape format

# Simulate 50 gene trees within the species tree
# Batch of gene trees (more efficient)
gene_trees_obs <- simulate_dtl_batch( #_obs is for "observed"
  species_tree = sp_tree,
  n = 10L,          # Number of gene trees
  lambda_d = 0.2,
  lambda_t = 0.3,
  lambda_l = 0.2,
  transfer_alpha = 0,
  require_extant = TRUE,
  seed = 123L
)
#gene_trees_obs_ape <- lapply(gene_trees_obs, function(gt) read.tree(text = tree_to_newick(gt))) # ape format

# CLEAN THE GENE TREES TO ONLY KEEP VISIBLE (extant) TIPS (i.e. NO GHOSTS)
gene_trees_obs_ape_clean <- lapply(gene_trees_obs, keep_extant_tips_gn_full)
class(gene_trees_obs_ape_clean) <- "multiPhylo" # to make it a multiPhylo object (list of phylo objects)



# -----------------------------
# 2. PERFORM SIMULATIONS AND COMPARE SIMULATD AND 
#   "observed" datasets
# -----------------------------

# Paramètres
d_max <- 1
t_max <- 1
l_max <- 1
N_sim <- 5000       # nombre de simulations
N_keep <- 200       # nombre de points acceptés
K_diff <- 20        # nombre de différenciations aléatoires

# Tirage aléatoire de theta
sample_theta <- function() {
  c(
    d = runif(1, 0, d_max),
    t = runif(1, 0, t_max),
    l = runif(1, 0, l_max)
  )
}


