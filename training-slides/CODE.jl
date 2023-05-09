# import the necessary packages
using COBREXA, GLPK

download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

m = load_model("e_coli_core.json")

flux_balance_analysis_dict(m, GLPK.Optimizer)

m2 = convert(StandardModel, m)

m2.reactions["EX_o2_e"].lb=-10.0
m2.reactions["EX_o2_e"].ub=0.0
flux_balance_analysis_dict(m2, GLPK.Optimizer)["BIOMASS_Ecoli_core_w_GAM"]

println.(reactions(m));


flux_balance_analysis_dict(
  m, GLPK.Optimizer,
  modifications=[
    change_constraint("EX_o2_e", lb=0.0, ub=0.0),
    change_constraint("EX_h2o_e", lb=-5.0, ub=5.0),
  ],
)["BIOMASS_Ecoli_core_w_GAM"]

using JuMP

flux_balance_analysis_dict(
  m, GLPK.Optimizer,
  modifications=[
    change_constraint("EX_o2_e", lb=0.0, ub=0.0),
    change_objective("EX_co2_e", sense=JuMP.MAX_SENSE),
  ],
)["EX_co2_e"]

reactions(m) .=> flux_balance_analysis_vec(m, GLPK.Optimizer)

flux_variability_analysis(m,
  GLPK.Optimizer,
  bounds = objective_bounds(0.9),
  modifications=[
    change_constraint("EX_o2_e", lb=-10.0, ub=0.0),
    change_constraint("EX_glc__D_e", lb=-1.0, ub=0.0),
    knockout("b0727"),
  ],
)

flux_variability_analysis(m,
  GLPK.Optimizer,
  bounds = o -> (0.66*o, 0.66*o),
  modifications=[],
)

choke_oxygen_exchange(choke) = (mod, opt) -> begin
  r = first(filter(startswith("EX_o2"), reactions(mod)))
	change_constraint(r, lb=-choke, ub=0.0)(mod, opt)
end

flux_balance_analysis_dict(m,
  GLPK.Optimizer,
  modifications=[
    choke_oxygen_exchange(5),
        knockout("b0727"),
  ],
)["BIOMASS_Ecoli_core_w_GAM"]

screen(
  m,
  args = tuple.(collect(-19:0)),
  analysis = (m, arg) ->
  	arg => flux_balance_analysis_dict(m,
      GLPK.Optimizer,
      modifications=[
        change_constraint("EX_o2_e", lb=arg, ub=0.0),
      ],
    )["BIOMASS_Ecoli_core_w_GAM"],
)

screen(
  m,
  args = tuple.(reactions(m)),
  analysis = (m, reac) -> begin
  	r = flux_balance_analysis_dict(m,
      GLPK.Optimizer,
      modifications=[
        change_constraint(reac, lb=0.0, ub=0.0),
      ],
    )
    reac => if isnothing(r)
      nothing
    else 
      r["BIOMASS_Ecoli_core_w_GAM"]
    end
  end,
)

res = screen(
  m,
  args = tuple.(genes(m)),
  analysis = (m, gene) -> begin
  	r = flux_balance_analysis_dict(m,
      GLPK.Optimizer,
      modifications=[
        knockout(gene),
      ],
    )
    gene_name(m, gene) => if isnothing(r)
      nothing
    else 
      r["BIOMASS_Ecoli_core_w_GAM"]
    end
  end,
)

[ gene for (gene,objv) in res if isnothing(objv) ]



using Distributed
addprocs(5)
@everywhere using COBREXA, GLPK

@everywhere choke_oxygen_exchange(choke) = (mod, opt) -> begin
  r = first(filter(startswith("EX_o2"), reactions(mod)))
	change_constraint(r, lb=-choke, ub=0.0)(mod, opt)
end

res = screen(
  m,
  args = [ (i,j) for i=genes(m), j=genes(m)],
  analysis = (m, gene1, gene2) -> begin
  	r = flux_balance_analysis_dict(m,
      GLPK.Optimizer,
      modifications=[
        knockout([gene1, gene2]),
        choke_oxygen_exchange(5),
      ],
    )
    if isnothing(r)
      nothing
    else 
      r["BIOMASS_Ecoli_core_w_GAM"]
    end
  end,
  workers=workers(),
)


# choking all of the reactions to 90%, 80%, ...
res = screen(
  m,
  args = [ (i,j) for i=reactions(m), j=10:-1:0],
  analysis = (m, rid, amount) -> begin
  	r = flux_balance_analysis_dict(m,
      GLPK.Optimizer,
      modifications=[
        change_constraint(rid, lb=-amount, ub=amount),
      ],
    )
    if isnothing(r)
      nothing
    else 
      r["BIOMASS_Ecoli_core_w_GAM"]
    end
  end,
)

reactions(m) .=> res[:,end]

using UnicodePlots
heatmap([ isnothing(x) ? 0 : x for x=res ])


struct MyTrivialModel <: ModelWrapper
  inner :: MetabolicModel
end

COBREXA.unwrap_model(x::MyTrivialModel) = x.inner
function COBREXA.bounds(x::MyTrivialModel)
  lb,ub = bounds(x.inner)
  (0.5 .* lb, 0.5 .* ub)
end

m = MyTrivialModel(load_model("e_coli_core.json"))
flux_balance_analysis(m, GLPK.Optimizer)

struct ScaledModel <: ModelWrapper
	inner::MetabolicModel
	ratio::Float64
end

COBREXA.unwrap_model(x::ScaledModel) = x.inner
COBREXA.bounds(x::ScaledModel) =
	let (lbs, ubs) = bounds(x.inner)
		(x.ratio .* lbs, x.ratio .* ubs)
	end

COBREXA.balance(x::ScaledModel) = x.ratio .* balance(x.inner)

m = ScaledModel(load_model("e_coli_core.json"), 0.3)

flux_balance_analysis_dict(
  	m, GLPK.Optimizer
)["BIOMASS_Ecoli_core_w_GAM"]

save_model(convert(JSONModel, m), "test.json")

struct MagicModel <: ModelWrapper
	inner::MetabolicModel
	amount::Float64
end

using SparseArrays

COBREXA.unwrap_model(x::MagicModel) = x.inner
COBREXA.reactions(x::MagicModel) =
	[reactions(x.inner); "MagicReaction"]
COBREXA.bounds(x::MagicModel) =
	let (lbs, ubs) = bounds(x.inner)
		([lbs; -x.amount], [ubs; -x.amount])
	end
function COBREXA.stoichiometry(x::MagicModel)
	orig = stoichiometry(x.inner)
  nm = n_metabolites(x.inner)
  atp_idx = first(indexin(["atp_c"], metabolites(m.inner)))
  adp_idx = first(indexin(["adp_c"], metabolites(m.inner)))
  pi_idx = first(indexin(["pi_c"], metabolites(m.inner)))
  column = spzeros(nm)
  column[atp_idx] = 1
  column[adp_idx] = -1
  column[pi_idx] = -1
  [orig column]
end
COBREXA.objective(x::MagicModel) =
	sparse([objective(x.inner); 0])[:,begin]

# This was the mistake I did during the course -- completely forgot about there needs
# to be a mapping of "fluxes" (actual human-readable fluxes of things) to "reactions"
# (low-level variables that may not correspond to actual "reactions" as in certain
# models there's e.g. splitting to isozymes and similar things)
using LinearAlgebra
COBREXA.reaction_flux(x::MagicModel) = sparse(I,n_reactions(x), n_reactions(x))
COBREXA.fluxes(x::MagicModel) = reactions(x)
# For simplicity, here we just default to fluxes == reaction variables.

m = MagicModel(load_model("e_coli_core.json"), -10)

flux_balance_analysis_dict(
  	m, GLPK.Optimizer
)["BIOMASS_Ecoli_core_w_GAM"]


# equal flux through 2 reactions (incomplete)
struct EqualFlux <: ModelWrapper
  inner :: MetabolicModel
  rid1 :: String
  rid2 :: String
end

m = load_model("e_coli_core.json")
m = MagicModel(m, 10)
m = EqualFlux(m, "EX_o2_e", "EX_co2_e")
...
flux_balance_analysis_dict(
  	m, GLPK.Optimizer
)["BIOMASS_Ecoli_core_w_GAM"]

# "how to model a community in MetabolicModel"
struct CommunityOfTwo <:ModelWrapper
  inner :: MetabolicModel
end

COBREXA.unwrap_model(x::CommunityOfTwo) = x.inner

COBREXA.reactions(x::CommunityOfTwo) =
	[ "1_" .* reactions(x.inner); "2_" .* reactions(x.inner)]

COBREXA.metabolites(x::CommunityOfTwo) =
  [ "1_" .* metabolites(x.inner); "2_" .* metabolites(x.inner)]

COBREXA.bounds(x::CommunityOfTwo) =
	let (lbs, ubs) = bounds(x.inner)
		([lbs; lbs], [ubs; ubs])
	end

function COBREXA.stoichiometry(x::CommunityOfTwo)
	s = stoichiometry(x.inner)
  m,n = size(s)
  [s spzeros(m,n);
   spzeros(m,n) s]
end

COBREXA.objective(x::CommunityOfTwo) =
	sparse([objective(x.inner); objective(x.inner)])[:,begin]

m = CommunityOfTwo(load_model("e_coli_core.json"))






