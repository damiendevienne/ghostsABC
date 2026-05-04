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


library(ape)
library(phangorn)

# Function get_overlap_scores computes, for each node of the gene tree, a measure of how much overlap there is in the identity of species on the left and right of the node. 
# The score is the intersection in species over the union of species on the left and right.
#It is a (rough) measure of how much the gene tree node corresponds to a speciation event (low score) or to a duplication event (high score).
fake_bootstrap <-function(sptr, gntr, N=100) {
  get_N_difftrees <- differentiation_fast(gntr, 100)
  sp_to_keep_index <- match(get_N_difftrees[[1]]$tip.label, sptr$tip.label)
  sptr2 <- keep.tip(sptr, sp_to_keep_index) 
  mat <- do.call(rbind, lapply(get_N_difftrees, function(x) prop.clades(sptr, x)))
  mat[is.na(mat)] <- 0
  colMeans(mat)
}



get_overlap_scores <- function(tree, return_mean=FALSE) {
  # extract species (prefix before "_")
  sp <- get_species(tree$tip.label)
  sp_levels <- unique(sp)
  k <- length(sp_levels)
  ntips <- Ntip(tree)
  nnodes <- Nnode(tree)
  # map species to indices
  sp_id <- match(sp, sp_levels)
  n <- ntips + nnodes
  pres <- matrix(FALSE, n, k)
  # initialize tips
  pres[1:ntips, ] <- diag(k)[sp_id, ]
  overlap <- rep(NA_real_, n)
  for (node in n:(ntips+1)) {
    ch <- Children(tree, node)
    v1 <- pres[ch[1], ]
    v2 <- pres[ch[2], ]
    # propagate
    pres[node, ] <- v1 | v2
    # Jaccard overlap
    inter <- sum(v1 & v2)
    uni   <- sum(v1 | v2)
    overlap[node] <- if (uni == 0) 0 else inter / uni
  }
  if (return_mean) mean(overlap[(ntips+1):n], na.rm = TRUE) else overlap[(ntips+1):n]
}

# for each node in the species tree 
# find the node in the gene tree that best matchs it (in terms of overlap of species)
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}


best_match_split <- function(sptr, gntr) {
  gntr$tip.label <- get_species(gntr$tip.label) # we only keep the species name (prefix before "_") in the tip labels of the gene tree
  sptip_gn <- gntr$tip.label
  gene_desc <- Descendants(gntr)
  gene_desc_sp <- lapply(gene_desc, function(x) {
    unique(sptip_gn[x])
  })
  sptip_sp <- sptr$tip.label
  species_desc <- Descendants(sptr)
  species_desc_sp <- lapply(species_desc, function(x) {
    unique(sptip_sp[x])
  })
  #remove tips
  species_desc_sp <- species_desc_sp[(Ntip(sptr)+1):length(species_desc_sp)]
  ns <- length(species_desc_sp)
  ng <- length(gene_desc_sp)
  M <- matrix(0, ns, ng)  
  for (i in seq_len(ns)) {
    Si <- species_desc_sp[[i]]
    for (j in seq_len(ng)) {
      M[i, j] <- jaccard(Si, gene_desc_sp[[j]])
    }
  }
  best_score <- apply(M, 1, max)
  mean(best_score)
}

sample_theta <- function(N = 1) {

  d_max <- 0.6
  t_max <- 0.6
  l_max <- 0.6

  THETAS <- matrix(NA_real_, N, 3)
  colnames(THETAS) <- c("d","t","l")

  i <- 1

  while (i <= N) {

    d <- runif(1, 0, d_max)
    t <- runif(1, 0, t_max)
    l <- runif(1, 0, l_max)

    if ((d + t) < (l + 0.1)) {

      THETAS[i, ] <- c(d, t, l)
      i <- i + 1

    }
  }

  as.data.frame(THETAS)
}




# Tirage aléatoire de theta
sample_theta_old <- function(N = 1) {
  d_max <- 0.8
  t_max <- 0.8
  l_max <- 0.8
  
  data.frame(
    d = runif(N, 0, d_max),
    t = runif(N, 0, t_max),
    l = runif(N, 0, l_max)
  )

}

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
toape <- function(tree) {
    treeape <- read.tree(text = tree_to_newick(tree))
    treeape 
}


#############UPDATED
keep_extant_tips_gn_full_da_trees <- function(tree) {
#    treeape <- read.tree(text = tree_to_newick(tree))
    tree.out <- keep.tip(tree, grep("e",get_species(tree$tip.label))) #nope... We still miss some extinct tips that are named e_becau
    if (is.null(tree.out)) return(NULL) # nothing left in the tree after pruning 
    if (Ntip(tree.out)<2) return(NULL) # if there are less than 2 tips, we cannot compute summary statistics, so we return NULL.
    x <- plotPhyloCoor(tree.out)[1:Ntip(tree.out),1]
    tol <- 1e-3  # adjust if needed
    keep <- which(x >= (max(x) - tol))
    tree.out.final <- keep.tip(tree.out, keep)
    tree.out.final
}

augment_tree <- function(sptr, lambda, mu, rho, niter, nthin) {
  # Write the tree to a temporary Newick file
  temp_file <- tempfile(fileext = ".tre")
  write.tree(sptr, file = temp_file)
  
  # Call the Julia script to augment the tree
  result <- system2("julia", args = c("AugmenTree.jl", temp_file, lambda, mu, rho, niter, nthin), stdout = TRUE)
  
  # Read the augmented tree from the output
  augmented_tree <- read.tree(file=paste0(result, ".trees"))
  
  # Clean up the temporary file
  unlink(temp_file)
  
  return(augmented_tree)
}

augment_tree_batch <- function(sptr, lambda, mu, rho_values, niter = 1, nthin = 1) {
  temp_tree <- tempfile(fileext = ".tre")
  write.tree(sptr, file = temp_tree)
  rho_file <- tempfile()
  writeLines(as.character(rho_values), rho_file)
  result <- system2(
    "julia",
    args = c("AugmenTree_batch.jl",
             temp_tree,
             lambda, mu,
             rho_file,
             niter, nthin),
    stdout = TRUE
  )
  # Read the augmented tree from the output
  tree_files <- paste0(result, ".trees")
  augmented_trees <- sapply(tree_files, read.tree, simplify = FALSE)
  class(augmented_trees) <- "multiPhylo"
  # rho_values_out <- as.numeric(gsub(".*rho_(.*?)_rep.*", "\\1", result))
  rho_values_out <- rho_values
  unlink(c(temp_tree, rho_file, tree_files))
  list(trees = augmented_trees,
       rhos = rho_values_out)
}


augment_tree_batch_old <- function(sptr, lambda, mu, n_rho_replicates, niter=1, nthin=1, target_ntips_max = 1000) {
  # Write the tree to a temporary Newick file
  temp_file <- tempfile(fileext = ".tre")
  write.tree(sptr, file = temp_file)
  
  # Call the Julia script to augment the tree
  result <- system2("julia", args = c("AugmenTree_batch.jl", temp_file, lambda, mu, n_rho_replicates, niter, nthin, target_ntips_max), stdout = TRUE) # nolint
  
  # Read the augmented tree from the output
  augmented_trees <- sapply(paste0(result, ".trees"), read.tree, simplify = FALSE)
  class(augmented_trees) <- "multiPhylo" # to make it a multiPhylo object (list of phylo objects)

  # extract rho values from the file names and assign them as names to the trees
  rho_values <- as.numeric(gsub(".*rho_(.*?)_rep.*","\\1",result))
  # Clean up the temporary file
  unlink(temp_file)
  
  return(list(trees = augmented_trees, rhos = rho_values))
}

toape <- function(tr) {
  read.tree(text = tree_to_newick(tr))
}

sim_gene_trees_fun <- function(sptr, n_gene_trees, d_val, t_val, l_val) {
  gene_trees <- simulate_dtl_batch_ape ( #_obs is for "observed"
    species_tree = sptr,
    n = n_gene_trees,          # Number of gene trees
    lambda_d = d_val,
    lambda_t = t_val,
    lambda_l = l_val,
    transfer_alpha = 0,
    require_extant = TRUE,
    seed <- sample.int(.Machine$integer.max, 1)
  )
  gene_trees_ape_clean <- lapply(gene_trees, keep_extant_tips_gn_full_da_trees) #remove unobservable tips
  # remove NULL trees (no extant tips after pruning)
  gene_trees_ape_clean <- gene_trees_ape_clean[!sapply(gene_trees_ape_clean, is.null)]
  # remove gene trees with less than 4 tips (some summary statistics cannot be computed on such small trees)
  gene_trees_ape_clean <- gene_trees_ape_clean[lapply(gene_trees_ape_clean, Ntip) >= 4] 
  class(gene_trees_ape_clean) <- "multiPhylo"
  return(gene_trees_ape_clean)
}


sim_gene_trees_fun_old <- function(sptr, n_gene_trees, d_val, t_val, l_val, seed=123L) {
  gene_trees <- simulate_dtl_batch ( #_obs is for "observed"
    species_tree = sptr,
    n = n_gene_trees,          # Number of gene trees
    lambda_d = d_val,
    lambda_t = t_val,
    lambda_l = l_val,
    transfer_alpha = 0,
    require_extant = TRUE,
    seed = seed
  )
  gene_trees_ape_clean <- lapply(gene_trees, keep_extant_tips_gn_full_da_trees) #remove unobservable tips
  # remove NULL trees (no extant tips after pruning)
  gene_trees_ape_clean <- gene_trees_ape_clean[!sapply(gene_trees_ape_clean, is.null)]
  # remove gene trees with less than 4 tips (some summary statistics cannot be computed on such small trees)
  gene_trees_ape_clean <- gene_trees_ape_clean[lapply(gene_trees_ape_clean, Ntip) >= 4] 
  class(gene_trees_ape_clean) <- "multiPhylo"
  return(gene_trees_ape_clean)
}

# Compute ks statistic between observed and simulated summary statistics for each summary statistic (column of the SCORES matrix)
compute_ks_vector <- function(obs_SCORES, sim_SCORES) {
  p <- ncol(obs_SCORES)
  ks <- numeric(p)
  for (j in 1:p) {
    x <- obs_SCORES[, j]
    y <- sim_SCORES[, j] 
    # KS statistic
    ks[j] <- suppressWarnings(ks.test(x, y)$statistic)
  }
  names(ks) <- colnames(obs_SCORES)
  return(ks)
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
    max_width              = treestats::max_width(tr), #need to specify treestats:: because of conflict with max_width from ggplot2
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
    phylogenetic_diversity = phylogenetic_diversity(tr, extinct_tol= 1e-8),
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
compute_all_scores_nonredondant <- function(tr) { #only scores that are not too much correlated
  c(
    colless                = colless(tr, normalization="yule"),
#    colless_corr           = colless_corr(tr, normalization = "yule"), # gives Inf values... 
    colless_quad           = colless_quad(tr,  normalization="yule"),
    ew_colless             = ew_colless(tr),
    double_cherries        = double_cherries(tr),
    four_prong             = four_prong(tr),
    stairs                 = stairs(tr),
    stairs2                = stairs2(tr),
    rquartet               = rquartet(tr, normalization="yule"),
    i_stat                 = i_stat(tr),
    j_one                  = j_one(tr),
    avg_ladder             = avg_ladder(tr),
    max_ladder             = max_ladder(tr),
    root_imbalance         = root_imbalance(tr),
    average_leaf_depth     = average_leaf_depth(tr),
    b2                     = b2(tr),
    mean_pair_dist         = mean_pair_dist(tr),
#    var_pair_dist          = var_pair_dist(tr), #does not correlate with any parameter value
    diameter               = diameter(tr, weight = FALSE),
    diameter_w             = diameter(tr, weight = TRUE),
    max_betweenness        = max_betweenness(tr, normalization="tips"),
    tree_height            = tree_height(tr),
    # gammaStat              = gammaStat(tr), # does not correlate with any parameter value
    pigot_rho              = pigot_rho(tr),
    phylogenetic_diversity = phylogenetic_diversity(tr, extinct_tol= 1e-8),
    mean_branch_length     = mean_branch_length(tr),
    var_branch_length      = var_branch_length(tr),
    var_branch_length_ext  = var_branch_length_ext(tr),
    var_branch_length_int  = var_branch_length_int(tr),
    treeness               = treeness(tr),
    Ntip                   = Ntip(tr),
    max_closeness          = max_closeness(tr, weight = FALSE, normalization = "tips"),
    max_closeness_w        = max_closeness(tr, weight = TRUE, normalization = "tips")
  )
}

compute_all_scores_nonredondant_nonorm <- function(tr) { #only scores that are not too much correlated
  c(
    colless                = colless(tr),
    colless_corr           = colless_corr(tr),
    colless_quad           = colless_quad(tr),
    ew_colless             = ew_colless(tr), ###YYEESS
    double_cherries        = double_cherries(tr),
    four_prong             = four_prong(tr),
    stairs                 = stairs(tr), ###YYEESS
    stairs2                = stairs2(tr), ###YYEESS
    rquartet               = rquartet(tr),
    i_stat                 = i_stat(tr), ###YYEESS
    j_one                  = j_one(tr), ###YYEESS
    avg_ladder             = avg_ladder(tr), ###YYEESS
    max_ladder             = max_ladder(tr),
    root_imbalance         = root_imbalance(tr),
    average_leaf_depth     = average_leaf_depth(tr),
    b2                     = b2(tr),
    mean_pair_dist         = mean_pair_dist(tr),
#    var_pair_dist          = var_pair_dist(tr), #does not correlate with any parameter value
    diameter               = diameter(tr, weight = FALSE),
    diameter_w             = diameter(tr, weight = TRUE), ###YYEESS
    max_betweenness        = max_betweenness(tr),
    tree_height            = tree_height(tr), ###YYEESS
    gammaStat              = gammaStat(tr), ### YYEESS
    pigot_rho              = pigot_rho(tr),
    phylogenetic_diversity = phylogenetic_diversity(tr, extinct_tol= 1e-8),
    mean_branch_length     = mean_branch_length(tr), ### YYEESS
    var_branch_length      = var_branch_length(tr), 
    var_branch_length_ext  = var_branch_length_ext(tr), ### YYEESS
    var_branch_length_int  = var_branch_length_int(tr),
    treeness               = treeness(tr), ###YYEESS
    Ntip                   = Ntip(tr),
    max_closeness          = max_closeness(tr, weight = FALSE),
    max_closeness_w        = max_closeness(tr, weight = TRUE)
  )
}

compute_informative_scores <- function(tr) {
  c(
    ew_colless             = ew_colless(tr),
    stairs                 = stairs(tr),
    stairs2                = stairs2(tr),
    i_stat                 = i_stat(tr),
    j_one                  = j_one(tr),
    avg_ladder             = avg_ladder(tr),
    tree_height            = tree_height(tr),
    gammaStat              = gammaStat(tr),
    mean_branch_length     = mean_branch_length(tr),
    treeness               = treeness(tr), 
    ntip                   = Ntip(tr)
  )
}

compute_all_scores_safe <- function(tr) {
  fns <- list(
    colless, colless_corr, colless_quad, blum, rogers, ew_colless,
    pitchforks, cherries, double_cherries, four_prong, stairs, stairs2,
    rquartet, i_stat, j_one, avg_ladder, max_ladder, root_imbalance,
    max_depth, sackin, average_leaf_depth, tot_path_length, tot_internal_path,
    avg_vert_depth, b1, b2, sym_nodes, max_width, max_del_width, mw_over_md,
    tot_coph, area_per_pair, mean_pair_dist, var_pair_dist, psv, wiener,
    diameter, max_betweenness, crown_age, tree_height, gammaStat, pigot_rho,
    phylogenetic_diversity, mean_branch_length, var_branch_length,
    var_branch_length_ext, var_branch_length_int, treeness, Ntip,
    function(x) max_closeness(x, weight = FALSE),
    function(x) max_closeness(x, weight = TRUE)
  )
  
  names(fns) <- c(
    "colless","colless_corr","colless_quad","blum","rogers","ew_colless",
    "pitchforks","cherries","double_cherries","four_prong","stairs","stairs2",
    "rquartet","i_stat","j_one","avg_ladder","max_ladder","root_imbalance",
    "max_depth","sackin","average_leaf_depth","tot_path_length","tot_internal_path",
    "avg_vert_depth","b1","b2","sym_nodes","max_width","max_del_width","mw_over_md",
    "tot_coph","area_per_pair","mean_pair_dist","var_pair_dist","psv","wiener",
    "diameter","max_betweenness","crown_age","tree_height","gammaStat","pigot_rho",
    "phylogenetic_diversity","mean_branch_length","var_branch_length",
    "var_branch_length_ext","var_branch_length_int","treeness","Ntip",
    "max_closeness","max_closeness_w"
  )
  
  res <- sapply(names(fns), function(nm) {
    tryCatch(fns[[nm]](tr), error = function(e) { message("FAILED: ", nm); NA })
  })
  return(res)
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

MAST<-function(t1,t2) {
    ##random name for file
    randname<-paste(sample(c(letters, LETTERS,0:9), 10, replace=TRUE),collapse="")
    infile<-paste(randname, ".nex", sep="")
    outfile<-paste(randname, ".out", sep="")
    trashfile<-paste(randname, ".trash", sep="")
    exec.agree<-function(t1,t2, infi, outfi,trashfi) {
        ##we must create a temp file with a unique name
        write.nexus(t1,t2,file=infile)
        cat(paste("BEGIN PAUP;\nAgree /all=sets showtree=no treefile=",outfi,";\nEND;\n",sep=""),file=infi, append=TRUE)
        system(paste("./paup -n ",infi," -u > ",trashfi,sep=""))
        system(paste("rm ",infi,sep=""))
       system(paste("rm ",trashfi,sep=""))
    }
    step1<-exec.agree(t1,t2,infile, outfile, trashfile)
    a<-read.nexus(outfile)
    system(paste("rm ",outfile, sep=""))
    if (class(a)=="multiPhylo") MASTsize<-Ntip(a[[1]])
    else MASTsize<-Ntip(a)
    return(list(size=MASTsize,trees=a))
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
#for a tree gntr get ndiff differentiations 
differentiation0 <- function(gntr, ndiff) {
  species_gn <- get_species(gntr$tip.label)
  groups <- split(gntr$tip.label, species_gn)
  spunik <- unique(species_gn)
  # common <- intersect(sptr$tip.label, names(groups))
  # sptr2 <- keep.tip(sptr, common)
  if (!length(spunik)>2) return(rep(NA, ndiff)) # not enough common species to compute distance)
  res <- replicate(ndiff, {
    # sample one per species
    chosen <- unlist(lapply(groups, sample, 1))
    # drop unused tips
    gntr_sampled <- keep.tip(gntr, chosen)
    gntr_sampled$tip.label <- get_species(gntr_sampled$tip.label)
    gntr_sampled
  }, simplify=FALSE)
  class(res) <- "multiPhylo"
  res
}

differentiation <- function(gntr, ndiff) {
  species_gn <- get_species(gntr$tip.label)
  spunik <- unique(species_gn)
  
  if (length(spunik) <= 2) return(rep(NA, ndiff))
  
  # precompute groups as indices (faster than labels)
  groups <- split(seq_along(gntr$tip.label), species_gn)
  ng <- length(groups)
  
  res <- vector("list", ndiff)
  
  for (i in seq_len(ndiff)) {
    # sample one index per species (no lapply/unlist)
    chosen_idx <- integer(ng)
    j <- 1
    for (g in groups) {
      chosen_idx[j] <- g[sample.int(length(g), 1)]
      j <- j + 1
    }
    
    # subset tree
    tr <- keep.tip(gntr, gntr$tip.label[chosen_idx])
    
    # relabel once
    tr$tip.label <- species_gn[chosen_idx]
    
    res[[i]] <- tr
  }
  
  class(res) <- "multiPhylo"
  res
}

fast_keep_tip <- function(phy, keep) {
  
  ape:::drop.tip.phylo(phy, drop, trim.internal = TRUE, collapse.singles = TRUE)
} 

# differentiation_fast <- function(gntr, ndiff) {
#   species_gn <- get_species(gntr$tip.label)
#   if (length(unique(species_gn)) <= 2) return(rep(NA, ndiff))
#   groups <- split(seq_along(gntr$tip.label), species_gn)
#   res <- vector("list", ndiff)  
#   for (i in seq_len(ndiff)) {
#     chosen_idx <- vapply(groups, function(g) g[sample.int(length(g), 1)], integer(1))
#     keep <- gntr$tip.label[chosen_idx]
#     drop <- setdiff(gntr$tip.label, keep)
#     tr <- drop.tip.phylo(gntr, drop, trim.internal = TRUE, collapse.singles = TRUE)
#     tr$tip.label <- get_species(tr$tip.label)    
#     res[[i]] <- tr
#   }
#   class(res) <- "multiPhylo"
#   res
# }
differentiation_fast <- function(gntr, ndiff) {
  species_gn <- get_species(gntr$tip.label)
  if (length(unique(species_gn)) <= 2) return(rep(NA, ndiff))  
  n <- length(gntr$tip.label)
  groups <- split(seq_len(n), species_gn)
  res <- vector("list", ndiff)  
  for (i in seq_len(ndiff)) {
    chosen_idx <- vapply(groups, function(g) g[sample.int(length(g), 1)], 1L)
    keep_mask <- logical(n)
    keep_mask[chosen_idx] <- TRUE    
    tr <- ape:::drop.tip.phylo(
      gntr,
      gntr$tip.label[!keep_mask],
      trim.internal = TRUE,
      collapse.singles = TRUE
    )
    # **FIX: match order after pruning**
    new_labels <- species_gn[match(tr$tip.label, gntr$tip.label)]
    tr$tip.label <- new_labels    
    res[[i]] <- tr
  }
  class(res) <- "multiPhylo"
  res
}




# these functions compare each gene tree to the species tree. The go a bit closer to reconciliation but without reconciliation. 
##add the mean depth of nodes as well. 
RF_multi_differentiation <- function(gntr, sptr, ndiff, method="PH85",FUN=min, norm=FALSE) {
  species_gn <- get_species(gntr$tip.label)
  groups <- split(gntr$tip.label, species_gn)
  common <- intersect(sptr$tip.label, names(groups))
  sptr2 <- keep.tip(sptr, common)
  if (!length(common)>2) return(rep(NA, ndiff)) # not enough common species to compute distance)
  res <- replicate(ndiff, {
    # sample one per species
    chosen <- unlist(lapply(groups, sample, 1))
    # drop unused tips
    gntr_sampled <- keep.tip(gntr, chosen)
    gntr_sampled$tip.label <- get_species(gntr_sampled$tip.label)    
    dtopo <- dist.topo(unroot(sptr2), unroot(gntr_sampled), method=method)
    if (norm) dtopo <- dtopo / (2*(length(common) - 3)) ## SOLVE NORMALIZATION ISSUE IF METHOD="SCORE" 
    #in addition to dtopo, we compare node depth between the two tres: 
    h_sp <- node.depth.edgelength(sptr2)
    maxh_sp <- max(h_sp)
    h_gn <- node.depth.edgelength(gntr_sampled)
    maxh_gn <- max(h_gn)
    h_sp <- h_sp[(Ntip(sptr2)+1):(Ntip(sptr2)+Nnode(sptr2))]/maxh_sp
    h_gn <- h_gn[(Ntip(gntr_sampled)+1):(Ntip(gntr_sampled)+Nnode(gntr_sampled))]/maxh_gn
    corr_depth <- ks.test(h_sp, h_gn)$statistic
    return(c(dtopo, corr_depth, mean(h_gn)))
  })
  res
}


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
plot_augmented_tree <- function(augmented_species_tree_ape) {

  # tips beginning with "e"
  e_tips <- grep("^e", augmented_species_tree_ape$tip.label)

  # root node
  root <- Ntip(augmented_species_tree_ape) + 1

  # edges to highlight
  edges_to_keep <- c()

  for (tip in e_tips) {

    path <- nodepath(augmented_species_tree_ape, from = root, to = tip)

    for (i in 1:(length(path) - 1)) {

      parent <- path[i]
      child  <- path[i + 1]

      edge_id <- which(
        augmented_species_tree_ape$edge[,1] == parent &
        augmented_species_tree_ape$edge[,2] == child
      )

      edges_to_keep <- c(edges_to_keep, edge_id)
    }
  }

  edges_to_keep <- unique(edges_to_keep)

  # edge styles
  edge_col <- rep("#929292", nrow(augmented_species_tree_ape$edge))
  edge_lwd <- rep(1.5, nrow(augmented_species_tree_ape$edge))

  edge_col[edges_to_keep] <- "black"
  edge_lwd[edges_to_keep] <- 3

  # labels: only e-tips
  tip_labels <- rep("", Ntip(augmented_species_tree_ape))
  tip_labels[e_tips] <- augmented_species_tree_ape$tip.label[e_tips]

  # temporary modified tree
  tr <- augmented_species_tree_ape
  tr$tip.label <- tip_labels

  plot(
    tr,
    edge.color = edge_col,
    edge.width = edge_lwd, 
    cex=2
  )
}


sim_gene_trees_fun_full <- function(sptr, n_gene_trees, d_val, t_val, l_val) {
  gene_trees <- simulate_dtl_pr_species_batch_ape ( #_obs is for "observed"
    species_tree = sptr,
    n = n_gene_trees,          # Number of gene trees
    lambda_d = d_val,
    lambda_t = t_val,
    lambda_l = l_val,
    transfer_alpha = 0,
    require_extant = TRUE,
    seed <- sample.int(.Machine$integer.max, 1)
  )
  gene_trees_ape_clean <- lapply(gene_trees, keep_extant_tips_gn_full_da_trees) #remove unobservable tips
  # remove NULL trees (no extant tips after pruning)
  gene_trees_ape_clean <- gene_trees_ape_clean[!sapply(gene_trees_ape_clean, is.null)]
  # remove gene trees with less than 4 tips (some summary statistics cannot be computed on such small trees)
  gene_trees_ape_clean <- gene_trees_ape_clean[lapply(gene_trees_ape_clean, Ntip) >= 4] 
  class(gene_trees_ape_clean) <- "multiPhylo"
  return(list(a = gene_trees, b = gene_trees_ape_clean))
}
