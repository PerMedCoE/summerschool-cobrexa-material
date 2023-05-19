
# Install the packages if they are not installed yet
import Pkg
Pkg.add(["COBREXA", "GLPK"])

using COBREXA

# Let's download and open a big model
import Downloads
Downloads.download("http://bigg.ucsd.edu/static/models/iJO1366.json", "ecoli.json")
m = load_model(StandardModel, "ecoli.json");

# ...we've opened it as StandardModel right away to allow easy manual
# modifications.

# see the objective function
objective(m) # gives a sparse vector that maximizes a single reaction
reactions(m)[19] # the number may differ

# See the exchanges.
# Unfortunately, instead of proper annotation the "standard" for identifying
# exchange reactions now depends on prefixing the reaction ID with EX_, so at
# least we can filter them out by checking the prefixes.
exchanges = filter(startswith("EX_"), reactions(m))

exchanges .=> reaction_name.(Ref(m), exchanges)

# Let's run the FBA. First we need everyone's favorite linear solver
# Ref: https://lcsb-biocore.github.io/COBREXA.jl/stable/examples/05a_fba/
using GLPK

# For simplicity, let's ask for a dictionary right away.
sol = flux_balance_analysis_dict(m, GLPK.Optimizer)

# we can have a look at our objective reaction
sol["BIOMASS_Ec_iJO1366_core_53p95M"]  # (it is possible to tab through the dictionary keys)

# we can have a look at the main things that happen
flux_summary(sol)
# (this guesses the exchange/biomass status based on reaction IDs)

# ...or write the result to a file for future use. First, let's create a data frame:
Pkg.add(["DataFrames", "CSV"])
using DataFrames, CSV
df = DataFrame(reaction = collect(keys(sol)), flux = collect(values(sol)))

# ...and write the CSV:
CSV.write("mysolution.csv", df)

# Let's choke the model a bit, reducing the availbale oxygen and sugar
m.reactions["EX_o2_e"].lb = -20
m.reactions["EX_glc__D_e"].lb = -2
sol = flux_balance_analysis_dict(m, GLPK.Optimizer)
sol["BIOMASS_Ec_iJO1366_core_53p95M"] # less growth

# allow eating acetate instead
m.reactions["EX_ac_e"].lb = -100
sol = flux_balance_analysis_dict(m, GLPK.Optimizer)
sol["BIOMASS_Ec_iJO1366_core_53p95M"] # a lot of growth again

# At this point, the original model data has been overwritten and there's no
# telling which bounds are still from the original model or which have been
# modified. For many reasons it is better to do this stuff without breaking the
# model internals manually, and COBREXA has a system of "analysis
# modifications" for that purpose:
m = load_model(StandardModel, "ecoli.json");
sol = flux_balance_analysis_dict(
    m,
    GLPK.Optimizer,
    modifications = [
        change_constraint("EX_o2_e", lb = -20.0),
        change_constraint("EX_glc__D_e", lb = -2.0),
    ],
);
sol["BIOMASS_Ec_iJO1366_core_53p95M"]
#...this scales much better if you need to try more stuff. Other modifications
# include e.g. `change_objective`, `silence` for shutting down the output from
# overly verbose solvers (such as OSQP), and `change_optimizer_attribute` for
# tuning the optimizer behavior.

# There is a nice app at https://escher.github.io/ that allows us to visualize
# and browse the solutions to metabolic models. You can load the visualization
# of this model as Map: "Core metabolism (e_coli_core)" and Tool: "Viewer".
# 
# Let's produce a JSON file with our solution that we can upload:

Pkg.add("JSON")
using JSON

write(
    "mysolution.json", 
    JSON.json(
        sol, # the solution
        2, # make the JSON human-readable by using 2-space indentation, instead of optimizing for size
    )
)

# You can now upload the file `mysolution.json` to the Escher viewer via
# `Data → Load reaction data`
# to see the fluxes visualized.

# More configurable Escher plotting directly from Julia to files (PDF, PNG) is
# available via https://github.com/stelmo/Escher.jl
