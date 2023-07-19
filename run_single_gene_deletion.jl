# Install COBREXA if not installed yet, and load it
import Pkg
Pkg.add("COBREXA")
using COBREXA

Pkg.add("GLPK")
using GLPK

# Let's download and open the big model again

path_to_model = "data/E_coli_iJO1366.json"
model = load_model(StandardModel, path_to_model);

Pkg.add("Distributed")
using Distributed
addprocs(4)  # only add processes if you are sure that you have sufficient resources available!

# load our stuff on the small cluster
# only necessary if you added the extra processes
@everywhere using COBREXA, GLPK 

# `screen` function allows us to run many analyses on a model with parallel,
# with many optimizations related for distributed processing (e.g., data are
# only moved once).

# knockout_fluxes = screen(
#     model, # the model which we process
#     args = tuple.(genes(model)), # all argument lists for the analyses
#     analysis = (model, gene) -> # the analysis functionthat we want to run on the cluster for each item from the argument list
#         flux_balance_analysis_dict(model, GLPK.Optimizer, modifications = [knockout(gene)]),
#     workers = workers(), # this gives it the list of worker nodes to use
# )

# Exercise: Try processing more genes with and without the workers parameter to
# see the speed-up.

# Now that we see that it works, let's post-process the results a little, and
# also add more genes:
knockout_fluxes = screen(
    model,
    args = tuple.(genes(model)),
    analysis = (model, gene) -> begin
        res = flux_balance_analysis_dict(model, GLPK.Optimizer, modifications = [knockout(gene)])
        if !isnothing(res)
            gene => res["BIOMASS_Ec_iJO1366_core_53p95M"]
        else
            gene => 0.0
        end
    end,
    workers = workers(),
)

# After everything works, you can erase the limit to the first 50 genes and see
# a complete result.

# Let's create a CSV with a report, as always
Pkg.add(["DataFrames","CSV"])
using DataFrames, CSV

df = DataFrame(gene = first.(knockout_fluxes))
df.name = gene_name.(Ref(model), df.gene)
df.fluxes = last.(knockout_fluxes)
df

# Typically we want to mark the genes that changed something. Let's mark the
# genes that are required for growth as essential, and the ones that reduce the
# growth somehow (but not fatally) as interesting.
best_result = maximum(last.(knockout_fluxes))
essential_threshold = 0.01 * best_result
df.essential = df.fluxes .<= essential_threshold
df.interesting = (df.fluxes .< best_result * 0.999) .&& .!df.essential
df

CSV.write("out/ko_report.csv", df)
