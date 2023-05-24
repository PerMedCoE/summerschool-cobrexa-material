
# Install the packages if they are not installed yet
import Pkg
Pkg.add(["COBREXA", "GLPK"])

using COBREXA

# Let's download and open a big model
import Downloads
Downloads.download("http://bigg.ucsd.edu/static/models/iJO1366.json", "ecoli.json")
model = load_model(StandardModel, "ecoli.json");

# ...we've opened it as StandardModel right away to allow easy manual
# modifications.

# see the objective function
objective(model) # gives a sparse vector that maximizes a single reaction
reactions(model)[19] # the number may differ

# See the exchanges.
# Unfortunately, instead of proper annotation the "standard" for identifying
# exchange reactions now depends on prefixing the reaction ID with EX_, so at
# least we can filter them out by checking the prefixes.
exchanges = filter(startswith("EX_"), reactions(model))

exchanges .=> reaction_name.(Ref(model), exchanges)

# Let's run the FBA. First we need everyone's favorite linear solver
# Ref: https://lcsb-biocore.github.io/COBREXA.jl/stable/examples/05a_fba/
using GLPK

# For simplicity, let's ask for a dictionary right away.
solution_dict = flux_balance_analysis_dict(model, GLPK.Optimizer)

# we can have a look at our objective reaction
solution_dict["BIOMASS_Ec_iJO1366_core_53p95M"]  # (it is possible to tab through the dictionary keys)

# we can have a look at the main things that happen
flux_summary(solution_dict)
# (this guesses the exchange/biomass status based on reaction IDs)

# ...or write the result to a file for future use. First, let's create a data frame:
Pkg.add(["DataFrames", "CSV"])
using DataFrames, CSV
df = DataFrame(reaction = collect(keys(solution_dict)), flux = collect(values(solution_dict)))

# ...and write the CSV:
CSV.write("mysolution.csv", df)

# Let's choke the model a bit, reducing the availbale oxygen and sugar
model.reactions["EX_o2_e"].lb = -20
model.reactions["EX_glc__D_e"].lb = -2
solution_dict = flux_balance_analysis_dict(model, GLPK.Optimizer)
solution_dict["BIOMASS_Ec_iJO1366_core_53p95M"] # less growth

# allow eating acetate instead
model.reactions["EX_ac_e"].lb = -100
solution_dict = flux_balance_analysis_dict(model, GLPK.Optimizer)
solution_dict["BIOMASS_Ec_iJO1366_core_53p95M"] # a lot of growth again

# At this point, the original model data has been overwritten and there's no
# telling which bounds are still from the original model or which have been
# modified. For many reasons it is better to do this stuff without breaking the
# model internals manually, and COBREXA has a system of "analysis
# modifications" for that purpose:
model = load_model(StandardModel, "ecoli.json");
solution_dict = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer,
    modifications = [
        change_constraint("EX_o2_e", lb = -20.0),
        change_constraint("EX_glc__D_e", lb = -2.0),
    ],
);
solution_dict["BIOMASS_Ec_iJO1366_core_53p95M"]
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
        solution_dict, # the solution
        2, # enforce nice 2-space indentation instead of squashing the JSON to optimize for size
    )
)

# You can now upload the file `mysolution.json` to the Escher viewer via
# `Data → Load reaction data`
# to see the fluxes visualized.

# More configurable Escher plotting directly from Julia to files (PDF, PNG) is
# available via https://github.com/stelmo/Escher.jl
