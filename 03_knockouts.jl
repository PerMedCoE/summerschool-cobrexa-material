
using COBREXA

# Let's download and open the big model again
import Downloads
Downloads.download("http://bigg.ucsd.edu/static/models/iJO1366.json", "ecoli.json")
m = load_model(StandardModel, "ecoli.json");

# Main reference is this:
# https://lcsb-biocore.github.io/COBREXA.jl/stable/examples/07_gene_deletion/
#
# Each reaction has a gene association, which dictates the gene products that
# need to be available so that the reaction can "run".
genes(m)

grr = reaction_gene_association(m, "PFK")

# The result is in DNF for (computational) simplicity; the rules can be
# converted e.g. to Strings for reading:
COBREXA._unparse_grr(String, grr)

# We might knock out genes by running through the reactions and evaluating DNF.
# The knockout is available as a modification for simplicity:

gene_name(m, "b0720")

using GLPK
sol = flux_balance_analysis_dict(m, GLPK.Optimizer, modifications = [knockout("b0720")])
sol["BIOMASS_Ec_iJO1366_core_53p95M"]
# ...the model is still feasible but growth is basically zero.

# We can screen through all genes. One could simply write:
[
    flux_balance_analysis_dict(m, GLPK.Optimizer, modifications = [knockout(g)]) for
    g in genes(m)[1:10]
]
# ...but that would run for quite a bit of time, and does not always even
# return a solution! (for some knockouts, there's even no feasible solution)

# First, let's use COBREXA parallelization capabilities to make this bearably
# fast. We use Distributed package to run this over multiple processes:

using Distributed
addprocs(8)  # you may add more depending on your machine or cluster size

# load our stuff on the small cluster
@everywhere using COBREXA, GLPK

# screen function allows us to run many analyses on a model with parallel, with
# many optimizations related for distributed processing (e.g., data are only
# moved once).

knockout_fluxes = screen(
    m, # the model which we process
    args = tuple.(genes(m)[1:10]), # all argument lists for the analyses
    analysis = (m, gene) -> # the analysis function ("lambda") that we want to run on the cluster for each item from the argument list
        flux_balance_analysis_dict(m, GLPK.Optimizer, modifications = [knockout(gene)]),
    workers = workers(), # this gives it the list of worker nodes to use
)

# let's preprocess the results a little, and add more genes
knockout_fluxes = screen(
    m,
    args = tuple.(genes(m)[1:50]),
    analysis = (m, gene) -> begin
        res = flux_balance_analysis_dict(m, GLPK.Optimizer, modifications = [knockout(gene)])
        if !isnothing(res)
            gene => res["BIOMASS_Ec_iJO1366_core_53p95M"]
        else
            gene => 0.0
        end
    end,
    workers = workers(),
)

# after debugging, you can erase the limit to the first 50 genes

# create a CSV with a report
using DataFrames, CSV

df = DataFrame(gene = first.(knockout_fluxes))
df.fluxes = last.(knockout_fluxes)

best_result = maximum(last.(knockout_fluxes))
essential_threshold = 0.01 * best_result
df.essential = df.fluxes .<= essential_threshold
df.interesting = (df.fluxes .< best_result * 0.999) .&& .!df.essential

CSV.write("ko_report.csv", df)
