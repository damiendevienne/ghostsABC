# ghostsABC

This repository contains scripts and experiments for an ABC (Approximate Bayesian Computation) approach to ghost estimation along branches of a species tree.

## Overview

Implements methods to infer ghost lineages (unsampled or extinct populations) in phylogenetic species trees using approximate Bayesian computation techniques.

## General pipeline

### Protocol without ghosts, following ABC papers by Blerina

- Start from a given species tree and a collection of gene trees
- Simulate gene trees along the branches of the species tree with various parameters
- Compare simulated and observed gene trees
- Modify parameters slightly
- Repeat the two last steps until simulated and observed gene trees agree. This happens when parameters are the good ones 

For all this it is straightforward to make simulations so that we know the true parameteres and try to recover them.


### Protocol with ghosts, extension of Blerina's ABC approach

- Start from a given species tree and a collection of gene trees

## Contents

- Scripts for ABC analysis
- Experimental pipelines
- ...

## Usage

See individual scripts for documentation and usage examples.



# MEASURES TO USE FOR ABC

Protocol used for choosing features: 
I did simulations with varying d,t,l (uniformly sampled in [0,1]), 1 species tree, with 20 tips, birth=1, death=0, 10 tips, 20 gene trees, 100 times
I computed all measures present in the treestats package 
I removed those being too long to compute
I removed redundant one by computing correlations between meaures: everything above 0.95 was collapsed (keeping 1)
I selected measures for which at least one of d, t or l correlated (by eye) with the score.
I computed spearman correlation between d,t,l and scores
I removed scores where best spearman was...