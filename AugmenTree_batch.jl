# Julia function to augment a tree with ghosts using Tapestree
# Usage:
# julia script.jl <input_tree_file> <lambda_value> <mu_value> <n_rho> <n_iter_value> <n_thin_value>

using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn)) #remove Info prined by Tapestree

using Tapestree
using Random
using Distributions 

if length(ARGS) < 6
    error("Usage: julia script.jl <input_tree_file> <lambda_value> <mu_value> <rho_file> <n_iter_value> <n_thin_value>")
end

input_file   = ARGS[1]
lambda_value = parse(Float64, ARGS[2])
mu_value     = parse(Float64, ARGS[3])
rho_file     = ARGS[4]
n_iter_value = parse(Int, ARGS[5])
n_thin_value = parse(Int, ARGS[6])

# target_ntips_max = parse(Int, ARGS[7])

tree = read_newick(input_file)
base = replace(input_file, r"\.[^.]+$" => "")


# nt = ntips(tree)
# rho_min = nt / target_ntips_max
# rho_max = 1.0   # corresponds to N = nt

# u = rand(n_rho)

# rho_values = 1.0 ./ (
#     (1.0 / rho_min) .-
#     u .* ((1.0 / rho_min) - (1.0 / rho_max))
# )
#We sample rho values so that the number of tips at the end is uniformly distributed between nt and target_ntips_max (or close to this, but hard to get in practice)

# N_values = rand(Uniform(nt, target_ntips_max), n_rho)
# rho_values = nt ./ N_values

rho_values = parse.(Float64, readlines(rho_file))

# min, max = 0, 1
# rho_values = min .+ (max - min) .* rand(n_rho)

output_files = String[]

for (i, rho_value) in enumerate(rho_values)

    r, tv = insane_cbd(
        tree,
        λi = lambda_value,
        μi = mu_value,
        pupdp = (0.0, 0.0, 0.2),
        tρ = Dict("" => rho_value),
        niter = n_iter_value,
        nthin = n_thin_value
    )

    relabeled_trees = map(x -> sT_label(x, tree), tv)

    rho_str = round(rho_value, digits=6) |> string
    output_file = "$(base)_rho_$(rho_str)_rep$(i)" # .trees will be added by Tapestree automatically

    write_newick(relabeled_trees, output_file)

    push!(output_files, output_file)
end

for f in output_files
    println(f)
end

