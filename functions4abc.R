## D de Vienne 30/01/2026

require(ape)
require(treestats)

# -----------------------------
# 1. Load Rust library and R wrappers
# -----------------------------

# Load the library
dyn.load("/home/ddevienne/Documents/CRETE/local-work/DTL/rustree/target/release/librustree.so")

# Source the wrapper functions
source("/home/ddevienne/Documents/CRETE/local-work/DTL/rustree/R/rustree.R")


###########
### ABC ###
###########


# -----------------------------
# 2. Fonctions
# -----------------------------
# return TRUE or FALSE depending on whether the gene tree made it to present time
gene_tree_extant <- function(tree, pres_time) {
  tol <- 1e-4  # small tolerance
  return(any(tree$depth >= pres_time - tol))
}

# Extract species from gene trees tip labels
get_species <- function(labels) {
  sub("_.*", "", labels)
}

# Function to only keep extant tips in specis and gene trees
keep_extant_tips_sp <- function(tree, treeape) {
    tips2keep <- as.character(tree$events$node_id[tree$events$event_type=="Leaf"])
    tree.out <- keep.tip(treeape, tips2keep)
    tree.out
}

keep_extant_tips_gn <- function(tree, treeape) {
    leavesinsp <- tree$species_tree$name[tree$species_tree$bd_event=="Leaf"]
    leavesingn <- tree$name[tree$event=="Leaf"]
    tips2keep <- leavesingn[is.element(get_species(leavesingn), leavesinsp)]
    tree.out <- keep.tip(treeape, tips2keep)
    tree.out
}

keep_extant_tips_gn_full <- function(tree) {
    leavesinsp <- tree$species_tree$name[tree$species_tree$bd_event=="Leaf"]
    leavesingn <- tree$name[tree$event=="Leaf"]
    tips2keep <- leavesingn[is.element(get_species(leavesingn), leavesinsp)]
    treeape <- read.tree(text = tree_to_newick(tree))
    tree.out <- keep.tip(treeape, tips2keep)
    tree.out
}


###############################################
## FUNCTION TIP-BLIND and SPTREE-INDEPENDENT ##
###############################################

# Number of tips in the  gene tree
n_tips_in_tree <- function(gt) {
  length(gt$tip.label)
}

# Number of species in a gene tree (different from number tips) -- here gene tree is in phylo format.
n_species_in_tree <- function(gt) {
  sp <- get_species(gt$tip.label)
  length(unique(sp))
}

compute_all_scores <- function(tr) {
  c(
    colless                = colless(tr),
    colless_corr           = colless_corr(tr),
    colless_quad           = colless_quad(tr),
    blum                   = blum(tr),
    rogers                 = rogers(tr),
    ew_colless             = ew_colless(tr),
    pitchforks             = pitchforks(tr), #remove!
    cherries               = cherries(tr),
    double_cherries        = double_cherries(tr),
    four_prong             = four_prong(tr),
    stairs                 = stairs(tr),
    stairs2                = stairs2(tr),
    rquartet               = rquartet(tr),
    i_stat                 = i_stat(tr),
    j_one                  = j_one(tr),
    avg_ladder             = avg_ladder(tr),
    max_ladder             = max_ladder(tr),
    root_imbalance         = root_imbalance(tr),
    max_depth              = max_depth(tr),
    sackin                 = sackin(tr),
    average_leaf_depth     = average_leaf_depth(tr),
    tot_path_length        = tot_path_length(tr),
    tot_internal_path      = tot_internal_path(tr),
    avg_vert_depth         = avg_vert_depth(tr),
    b1                     = b1(tr),
    b2                     = b2(tr),
    sym_nodes              = sym_nodes(tr),
    max_width              = max_width(tr),
    max_del_width          = max_del_width(tr),
    mw_over_md             = mw_over_md(tr),
    tot_coph               = tot_coph(tr),
    area_per_pair          = area_per_pair(tr),
    mean_pair_dist         = mean_pair_dist(tr),
    var_pair_dist          = var_pair_dist(tr),
    psv                    = psv(tr),
    wiener                 = wiener(tr),
    diameter               = diameter(tr),
    max_betweenness        = max_betweenness(tr),
    crown_age              = crown_age(tr),
    tree_height            = tree_height(tr),
    gammaStat              = gammaStat(tr),
    pigot_rho              = pigot_rho(tr),
    phylogenetic_diversity = phylogenetic_diversity(tr),
    mean_branch_length     = mean_branch_length(tr),
    var_branch_length      = var_branch_length(tr),
    var_branch_length_ext  = var_branch_length_ext(tr),
    var_branch_length_int  = var_branch_length_int(tr),
    treeness               = treeness(tr),
    Ntip                   = Ntip(tr),
    max_closeness          = max_closeness(tr, weight = FALSE),
    max_closeness_w        = max_closeness(tr, weight = TRUE)
  )
}


compute_all_scores2 <- function(tr) {
  c(
    mean_branch_length     = mean_branch_length(tr),
    var_branch_length      = var_branch_length(tr),
    var_branch_length_ext  = var_branch_length_ext(tr),
    var_branch_length_int  = var_branch_length_int(tr),
    stairs                 = stairs(tr),
    j_one                  = j_one(tr),
    var_pair_dist          = var_pair_dist(tr),
    max_closeness          = max_closeness(tr, weight = FALSE),
    max_closeness_w        = max_closeness(tr, weight = TRUE),
    colless_corr           = colless_corr(tr),
    ew_colless             = ew_colless(tr),
    crown_age              = crown_age(tr),
    average_leaf_depth     = average_leaf_depth(tr),
    mean_pair_dist         = mean_pair_dist(tr),
    gammaStat              = gammaStat(tr),
    phylogenetic_diversity = phylogenetic_diversity(tr),
    pigot_rho              = pigot_rho(tr),
    avg_ladder             = avg_ladder(tr),
    max_ladder             = max_ladder(tr),
    tree_height            = tree_height(tr)
  )
}


compute_all_scores2_timed <- function(tr) {
  
  timed <- function(expr) {
    t <- system.time(val <- eval(expr))
    list(value = val, time = unname(t["elapsed"]))
  }
  
  res <- list(
    mean_branch_length     = timed(quote(mean_branch_length(tr))),
    var_branch_length      = timed(quote(var_branch_length(tr))),
    var_branch_length_ext  = timed(quote(var_branch_length_ext(tr))),
    var_branch_length_int  = timed(quote(var_branch_length_int(tr))),
    stairs                 = timed(quote(stairs(tr))),
    j_one                  = timed(quote(j_one(tr))),
    var_pair_dist          = timed(quote(var_pair_dist(tr))),
    max_closeness          = timed(quote(max_closeness(tr, weight = FALSE))),
    max_closeness_w        = timed(quote(max_closeness(tr, weight = TRUE))),
    colless_corr           = timed(quote(colless_corr(tr))),
    ew_colless             = timed(quote(ew_colless(tr))),
    crown_age              = timed(quote(crown_age(tr))),
    average_leaf_depth     = timed(quote(average_leaf_depth(tr))),
    mean_pair_dist         = timed(quote(mean_pair_dist(tr))),
    gammaStat              = timed(quote(gammaStat(tr))),
    phylogenetic_diversity = timed(quote(phylogenetic_diversity(tr))),
    pigot_rho              = timed(quote(pigot_rho(tr))),
    avg_ladder             = timed(quote(avg_ladder(tr))),
    max_ladder             = timed(quote(max_ladder(tr))),
    tree_height            = timed(quote(tree_height(tr)))
  )
  
  values <- sapply(res, `[[`, "value")
  times  <- sapply(res, `[[`, "time")
  
  list(values = values, times = times)
}




library(corrplot)
## plot correations 
plot_corr <- function(score) {
  cor_mat <- cor(score, use = "pairwise.complete.obs")
  corrplot(cor_mat, 
           method = "color",
           type = "upper",
           addCoef.col = "black",
           number.cex = 0.7,
           tl.cex = 0.8,
           tl.col = "black")
}
# + All measures included from the treestats package
# calc_all_stats(gt)

###############################
## FUNCTION SPTREE-DEPENDENT ##
###############################
# these functions compare each gene tree to the species tree. The go a bit closer to reconciliation but without reconciliation. 















# Création d'une différenciation aléatoire pour una paire d'arbres
# Création d'une différenciation aléatoire pour una paire d'arbres
RF_on_random_differentiation <- function(t1, t2, ndiff, method="PH85",FUN=min, norm=FALSE) {
  sp1 <- get_species(t1$tip.label)
  sp2 <- get_species(t2$tip.label)
  common_species <- intersect(sp1,sp2)
  if (length(common_species) <= 3) {
    return(NA)  # Pas assez d'espèces communes pour calculer une distance RF
  }
  d <- vapply(seq_len(ndiff), function(i) {
    p1 <- sample.int(length(sp1))
    p2 <- sample.int(length(sp2))

    i1 <- p1[match(common_species, sp1[p1])]
    i2 <- p2[match(common_species, sp2[p2])]

    t1.diff <- keep.tip(t1, t1$tip.label[i1])
    t2.diff <- keep.tip(t2, t2$tip.label[i2])

    t1diffsp <- t1.diff
    t1diffsp$tip.label <- get_species(t1diffsp$tip.label)
    t2diffsp <- t2.diff
    t2diffsp$tip.label <- get_species(t2diffsp$tip.label)
    dtopo <- dist.topo(unroot(t1diffsp), unroot(t2diffsp), method=method)
    if (norm) dtopo / (2*(length(common_species) - 3)) else dtopo
  }, numeric(1))
  FUN(d)
}


# Compare distance to neirest neighbors in a tree
dist_to_nn <- function(tree) {
  mat <- cophenetic(tree)
  spnames <- get_species(rownames(mat))
  idx_by_species <- split(seq_along(spnames),
                          factor(spnames, levels = unique(spnames)))
  res <- vapply(seq_along(spnames), function(i) {
    same_sp_idx <- idx_by_species[[ spnames[i] ]]
    min(mat[i, -same_sp_idx])
  }, numeric(1))
  minperspecies <- unlist(lapply(idx_by_species, function(x,y) min(y[x]), y= res))
  return(minperspecies[order(names(minperspecies))]) #ordering for easy comparison later
}

# # Distance pour un jeu de n arbres
# distance_dataset <- function(G_sim, G_obs) {
#   mean(mapply(wRF_diff_min, G_sim, G_obs))
# }

# # -----------------------------
# # 3. Jeu de données observé (pour test)
# # -----------------------------
# S <- "SpeciesTreePlaceholder"
# G_obs <- simulate_gene_trees(S, c(d=0.05, t=0.02, l=0.03), n_genes=100)

# # -----------------------------
# # 4. Boucle ABC par rejet
# # -----------------------------
# results <- matrix(NA, nrow = N_sim, ncol = 4)
# colnames(results) <- c("d", "t", "l", "dist")

# for (i in seq_len(N_sim)) {

#   theta <- sample_theta()

#   G_sim <- simulate_gene_trees(S, theta)

#   dist <- distance_dataset(G_sim, G_obs)

#   results[i, ] <- c(theta, dist)
# }

# # -----------------------------
# # 5. Sélection des meilleurs
# # -----------------------------
# results <- results[order(results[, "dist"]), ]
# posterior <- results[1:N_keep, 1:3]

# # -----------------------------
# # 6. Diagnostic
# # -----------------------------
# pairs(posterior, main="Posterior of D,T,L")
# summary(posterior)
