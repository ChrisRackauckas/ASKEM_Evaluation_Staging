# Evaluation Scenario 2

```@example scenario2
cd(@__DIR__)
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
sys2 = structural_simplify(sys2)
```

```@example carcione
prob = ODEProblem(sys2, [], (0.0, 1000.0))
sol = solve(prob, Tsit5())

@unpack Infected = sys2
plot(sol, idxs = Infected)
```

### Sensitivity Analysis

```@example carcione
pbounds = [
    alpha => [0.003, 0.006],
    epsilon => [1 / 6, 1 / 2],
    gamma => [0.1, 0.2],
    mu => [0.01, 0.02],
    beta => [0.7, 0.9],
]
create_sensitivity_plot(prob, 100.0, Deceased, pbounds; samples = 2000)
```
