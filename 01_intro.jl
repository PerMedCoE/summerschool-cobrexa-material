
using COBREXA

# Let's first get a simple model to have a look at
import Downloads
Downloads.download(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
)

# Load the model
ecoli = load_model("e_coli_core.json");

# have a look at what the model contains. First, metabolites
metabolites(ecoli)

# reactions
reactions(ecoli)

# we can have a look at what each reaction does
reaction_stoichiometry(ecoli, "PFK")

# or look at which metabolites belong to which compartments (Julia looping!)
[m => metabolite_compartment(ecoli, m) for m in metabolites(ecoli)]

# or try to get human-readable names from the reactions
sort([r => reaction_name(ecoli, r) for r in reactions(ecoli)])

# this is useful to e.g. quickly look up dehumanized gene IDs
gs = genes(ecoli)
sort(gene_name.(Ref(ecoli), gs) .=> gs)

# the model's internals are described by a bipartite graph between reactions
# and metabolites, commonly stoichiometry
stoichiometry(ecoli)
# (the matrix is likely zoomed out, but otherwise it's a normal matrix)

# Let's try to build a tiny custom model
# ref: https://lcsb-biocore.github.io/COBREXA.jl/stable/examples/04b_standardmodel_construction/

# The model type for "manually constructed" models we call a StandardModel
# because it is kinda standard way to handle stuff in cobra community.
# Technically we can convert the above model to StandardModel and play with it:
ecoli = convert(StandardModel, ecoli)
ecoli.reactions["PFK"]
ecoli.reactions["PFK"].lb = 0.5

# But we want to make a completely custom model
m = StandardModel("MyModel")

# let's make some metabolites
a = Metabolite("a", name = "molecule A", formula = "H2O", compartment = "outside")
b = Metabolite("b") #details omitted for demonstration
c = Metabolite("c")

# Push the prepared metabolites into the model (the ! stands for "execute!", it
# is a syntactic convention for warning that the function changes the model in
# place)
add_metabolites!(m, [a, b, c])

# let's make some reactions
bs = [Reaction("b$i", lb = 0.0, ub = 10.0) for i = 1:3]
bs[1].metabolites = Dict("a" => 1)
bs[2].metabolites = Dict("b" => -1)
bs[3].metabolites = Dict("c" => -1)
vs = [Reaction("v$i", lb = 0.0, ub = 10.0) for i = 1:3]
vs[1].metabolites = Dict("a" => -1, "b" => 1)
vs[2].metabolites = Dict("a" => -1, "c" => 1)
vs[3].metabolites = Dict("a" => 1, "c" => -1)

# add the reactions to the model
add_reactions!(m, [bs; vs])

# let's have a look at what it looks like as a matrix
stoichiometry(m)
