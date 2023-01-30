# this script requires https://github.com/SciML/SBMLToolkit.jl/pull/108

using OrdinaryDiffEq, ModelingToolkit, EasyModelAnalysis, SBML, SBMLToolkit, UnPack, Plots
function sub_cont_ev(ev, rules)
    ModelingToolkit.SymbolicContinuousCallback(substitute(ev.eqs, rules),
                                               substitute(ev.affect, rules))
end
fn = "Giordano2020.xml"

myread(fn) = readSBML(fn, doc -> begin
                          set_level_and_version(3, 2)(doc)
                          convert_promotelocals_expandfuns(doc)
                      end)

m = myread(fn)
rn = ReactionSystem(m)
sys = convert(ODESystem, rn)
eqs = equations(sys)
defs_ = ModelingToolkit.defaults(sys)
defs = deepcopy(defs_)
evs = ModelingToolkit.get_continuous_events(sys)

# these are the constant=false params 
params_to_sub = unique(ModelingToolkit.lhss(vcat(map(x -> x.affect, evs)...)))

@unpack alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta = sys
ps = [alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta]

subd = Dict(params_to_sub .=> ps)
newevs = map(x -> sub_cont_ev(x, subd), evs)

sys2 = ODESystem(eqs, ModelingToolkit.get_iv(sys), states(sys), parameters(sys);
                 continuous_events = newevs, defaults = defs, name = nameof(sys))
prob = ODEProblem(sys2, [], (0.0, 1000.0))
sol = solve(prob, Tsit5())
plot(sol)
