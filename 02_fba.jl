
using COBREXA

# Let's download and open a big model
import Downloads
Downloads.download("http://bigg.ucsd.edu/static/models/iJO1366.json", "ecoli.json")
m = load_model(StandardModel, "ecoli.json");

# ...we've opened it as StandardModel right away to allow easy manual
# modifications.

# let's have a look at the genes
m.genes["b0722"]

# see the objective function
objective(m) # gives a sparse vector that maximizes a single reaction
reactions(m)[19] # the number may differ

# see the exchanges.
# Unfortunately the standard depends on prefixing the reaction ID with EX_, so
# at least we can filter the reactions to them.
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

# or write the result to a file for future use
using DataFrames, CSV
df = DataFrame(reaction = collect(keys(sol)), flux = collect(values(sol)))
CSV.write("mysolution.csv", df)

# Let's choke the model a bit, reducing the availbale oxygen and sugar
m.reactions["EX_o2_e"].lb = -20
m.reactions["EX_glc__D_e"].lb = -2
sol = flux_balance_analysis_dict(m, GLPK.Optimizer)
sol["BIOMASS_Ec_iJO1366_core_53p95M"] # less growth

# allow eating acetate instead
m.reactions["EX_ac_e"].lb = -100;
sol = flux_balance_analysis_dict(m, GLPK.Optimizer);
sol["BIOMASS_Ec_iJO1366_core_53p95M"] # a lot of growth again

# at this point, it is clear that the original model data has been overwritten
# so it might be much better to do this stuff reproducibly. COBREXA has
# modifications for that purpose.
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
#...this scales much better if you need to do more stuff.

# Escher plotting is available via https://github.com/stelmo/Escher.jl
