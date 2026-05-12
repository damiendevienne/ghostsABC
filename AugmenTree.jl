# Julia function to augment a tree with ghosts using Tapestree, and save the resulting tree(s) in Newick format.
# Usage: julia script.jl <input_tree_file> <lambda_value> <mu_value> <rho_value> <n_iter_value> <n_thin_value>


using Tapestree

if length(ARGS) < 6
    error("Usage: julia script.jl <input_tree_file> <lambda_value> <mu_value> <rho_value> <n_iter_value> <n_thin_value>")
end

input_file = ARGS[1]
lambda_value = parse(Float64, ARGS[2])
mu_value = parse(Float64, ARGS[3])
rho_value = parse(Float64, ARGS[4])
n_iter_value = parse(Int, ARGS[5])
n_thin_value = parse(Int, ARGS[6])

tree = read_newick(input_file)

r, tv = insane_cbd(
    tree,
    λi = lambda_value,
    μi = mu_value,
    pupdp = (0.0, 0.0, 0.2),
    tρ = Dict("" => rho_value), 
    niter=n_iter_value,
    nthin=n_thin_value
)

relabeled_tree = map(x -> sT_label(x, tree), tv)

base = replace(input_file, r"\.[^.]+$" => "")
rho_str = replace(string(rho_value), "." => "_")
output_file = "$(base)_rho$(rho_str)"

write_newick(relabeled_tree, output_file)

println(output_file)