
# Install COBREXA if not installed yet, and load it
import Pkg
Pkg.add("COBREXA")
using COBREXA

# Let's download and open the big model again
import Downloads
Downloads.download("http://bigg.ucsd.edu/static/models/iJO1366.json", "ecoli.json")
model = load_model(StandardModel, "ecoli.json");

# Main reference is this:
# https://lcsb-biocore.github.io/COBREXA.jl/stable/examples/07_gene_deletion/
#
# Each reaction has a gene association, which dictates the gene products that
# need to be available so that the reaction can "run".
genes(model)

gene_rules_dnf = reaction_gene_association(model, "PFK")

# The result is in DNF for (computational) simplicity; the rules can be
# converted e.g. to Strings which are more suitable for reading:
COBREXA._unparse_grr(String, gene_rules_dnf)

# We might knock out genes by running through the reactions and evaluating DNF.
# The knockout is available as a modification for simplicity:
gene_name(model, "b0720")

# We need a solver
Pkg.add("GLPK")

# Run the knockout (implemented as an analysis modification for convenience)
using GLPK
sol = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer,
    modifications = [knockout("b0720")],
)
sol["BIOMASS_Ec_iJO1366_core_53p95M"]
# ...the model is still feasible (so it can "sustain itself"), but growth is
# basically zero.

# We can screen through all genes. One could simply write something like:
[
    flux_balance_analysis_dict(model, GLPK.Optimizer, modifications = [knockout(g)]) for
    g in genes(model)[1:10]
]
# ...but that might be slow for larger amounts of genes, and we would like to
# add some special handling for knockouts where there is no feasible solution
# (and the function returns `nothing`).

# First, let's use COBREXA parallelization capabilities to make this faster. We
# use Distributed package to run this over multiple processes:

using Distributed
#addprocs(8)  # you may add more depending on your machine or cluster size

# load our stuff on the small cluster
#@everywhere using COBREXA, GLPK

# `screen` function allows us to run many analyses on a model with parallel,
# with many optimizations related for distributed processing (e.g., data are
# only moved once).

knockout_fluxes = screen(
    model, # the model which we process
    args = tuple.(genes(model)[1:10]), # all argument lists for the analyses
    analysis = (model, gene) -> # the analysis function ("lambda") that we want to run on the cluster for each item from the argument list
        flux_balance_analysis_dict(model, GLPK.Optimizer, modifications = [knockout(gene)]),
    workers = workers(), # this gives it the list of worker nodes to use
)

# Exercise: Try processing more genes with and without the workers parameter to
# see the speed-up.

# Now that we see that it works, let's post-process the results a little, and
# also add more genes:
knockout_fluxes = screen(
    model,
    args = tuple.(genes(model)[1:50]),
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

CSV.write("ko_report.csv", df)

# ## Doing the knockouts manually, the hard way
#
# Now, let's have a look at how the knockouts are computed.
#
# Each reaction has a Gene-Reaction Rule (GRR) that marks the genes required
# for it to actually work in the organism. These are generally any Boolean
# expressions, but in COBREXA we tend to store them in disjunctive normal form
# (DNF, see https://en.wikipedia.org/wiki/Disjunctive_normal_form) which
# closely corresponds to the biological meaning of gene units that form
# interchangeable complexes. You can access them using the `grr` field in
# Reaction structures:
model.reactions["RNDR2"].grr

# Here, the reaction can be supported by either of the 2 possibilities (enzyme
# complexes) where the first possibility is built from gene products of genes
# `b2234`, `b2235`, and `b2582`; and as the second possibility it can also use
# `b3781` instead of the `b2582`.

# We may list all GRRs simply by iterating through the model reactions:
[rid => r.grr for (rid,r) in model.reactions]

# It is often interesting to ask which reactions may depend on which gene, we
# can make a convenience function for that:
reactions_of_gene(model, gene) =
  [rid for (rid,r) in model.reactions if !isnothing(r.grr) && any(complex -> any(contains(gene), complex), r.grr)]

reactions_of_gene(model, "b1064")

# Using the vector notation is quite convenient for creating lists that allow
# us to get an overview of the situation:
gene_name.(Ref(model), genes(model)) .=> reactions_of_gene.(Ref(model), genes(model))

# Now, given a set of genes that we want to knock out, we can manually find if
# a given reaction will still work or not. Let's try on RNDR1:

grr = model.reactions["RNDR1"].grr

ko_genes = ["b2234"]

# We can transform the `grr` to a form where it says which genes are present
# and which genes are not:
grr_available = map(c -> map(!in(ko_genes), c), grr)

# To determine if the reaction _can_ work, at least one ("any") of the
# complexes must have all gene products available:
any(all, grr_available)

# Since `b2234` is essential for RNDR1 (it needs to be present in all complexes
# that may run the reaction), the reaction is effectively disabled by knocking
# out `b2234`.

# What if we knock out `b2582`?
ko_genes = ["b2582"]
grr_available = map(c -> map(!in(ko_genes), c), grr)

any(all, grr_available)

# ...the reaction may still work with just `b2582` knocked out.

# Anyway, if we knock out multiple genes, the reaction will cease to work again:
ko_genes = ["b2582", "b3781"]
grr_available = map(c -> map(!in(ko_genes), c), grr)
any(all, grr_available)

# We can formalize the knockout evaluation in a function
function is_reaction_knocked_out(model, reaction, ko_genes)
    grr = model.reactions[reaction].grr
    if isnothing(grr)
        return false # reactions without a gene-reaction rule happen spontaneously and cannot be knocked out
    end
    grr_available = map(c -> map(!in(ko_genes), c), grr)
    !any(all, grr_available)
end

# Now, we can manually modify the model to disable the reactions that would be
# knocked out by a certain gene combination:
ko_genes = ["b2582", "b3781"]
for (rid, r) = model.reactions
    if is_reaction_knocked_out(model, rid, ko_genes)
        r.lb = r.ub = 0.0
    end
end

# Does the model still grow?
sol = flux_balance_analysis_dict(model, GLPK.Optimizer)
sol["BIOMASS_Ec_iJO1366_core_53p95M"]

# ...which seems like the combination of the 2 genes was not essential at all.
