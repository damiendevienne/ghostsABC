# ghostsABC

This repository contains scripts and experiments for an ABC (Approximate Bayesian Computation) approach to ghost estimation along branches of a species tree.

## Overview

Implements methods to infer ghost lineages (unsampled or extinct populations) in phylogenetic species trees using approximate Bayesian computation techniques.

## General pipeline


### Protocol without ghosts, following ABC papers by Blerina

- Start from a given species tree and a collection of gene trees
- Simulate gene trees along the branches of the species tree with various parameters
- Compare simulated and observed gene trees
- Modify parameters slightly
- Repeat the two last steps until simulated and observed gene trees agree. This happens when parameters are the good ones 

For all this it is straightforward to make simulations so that we know the true parameteres and try to recover them.


### Protocol with ghosts, extension of Blerina's ABC approach

- Start from a given species tree and a collection of gene trees



## Contents

- Scripts for ABC analysis
- Experimental pipelines
- ...

## Usage

See individual scripts for documentation and usage examples.

