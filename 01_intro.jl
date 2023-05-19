
# Installation. In Julia REPL this can be done easier using the packaging mode
# (typing `]add COBREXA`)
import Pkg
Pkg.add("COBREXA")

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
# is a syntactic convention for warning that the function changes some of the
# parameters (the model) in place)
add_metabolites!(m, [a, b, c])

# let's make some reactions
b1 = Reaction("b1", lb = 0.0, ub = 10.0);
b1.metabolites = Dict("a" => 1);

b2 = Reaction("b2", lb = 0.0, ub = 10.0);
b2.metabolites = Dict("b" => -1);

b3 = Reaction("b3", lb = 0.0, ub = 10.0);
b3.metabolites = Dict("c" => -1);

# add the reactions to the model
add_reactions!(m, [b1, b2, b3])

# shortcut for the above, makes an array of 3 reactions at once
vs = [Reaction("v$i", lb = 0.0, ub = 10.0) for i = 1:3]

# ...and we can fill them in 
vs[1].metabolites = Dict("a" => -1, "b" => 1)
vs[2].metabolites = Dict("a" => -1, "c" => 1)
vs[3].metabolites = Dict("a" => 1, "c" => -1)

# add the reactions to the model
add_reactions!(m, vs)

# let's have a look at what it looks like as a matrix
stoichiometry(m)

# Now let's try to numerically solve the model. First, we need a linear solver;
# GLPK will do a good job for this purpose. Other alternatives include Tulip
# (native interior-point solver), OSQP and Clarabel (for quadratic problems),
# SCIP (free and fast), Clp, and very good commercial solvers include Gurobi
# and CPLEX (you will need to install the licenses for these manually).

Pkg.add("GLPK")
using GLPK

# "Solving" the model is finding optimum with respect to some reaction, thus we
# first need to choose the objective that the solver should actually optimize.
# Here, let's maximize the flux through the reaction "v1".
change_objective!(m, "v1")

# You can observe the result in the model reactions' objective coefficients
# (and fine-tune that manually if needed):
m.reactions["v1"].objective_coefficient

# Flux balance analysis finds a steady state flux through the model which is
# within the reaction flux lower and upper bounds:
solution = flux_balance_analysis(m, GLPK.Optimizer)

# The above returned a solved optimization model description. From that we can
# extract the actual value of the objective function:
solved_objective_value(solution)

# ...and the complete description of the flux through the reactions
flux_vector(m, solution)

# ...or a little better as a dictionary
flux_dict(m, solution)

# typically, one does all of this in one step with a shortcut function:
flux = flux_balance_analysis_dict(m, GLPK.Optimizer)

# As a check, we can verify that the amount of metabolites in the model indeed
# stays balanced; if the solver worked, this vector should be zero (or within
# the numerical tolerance of zero):
stoichiometry(m) * flux_vector(m, solution)
