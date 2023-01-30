# this script requires https://github.com/SciML/SBMLToolkit.jl/pull/108
# This is a model created on COPASI 4.27 (Build 217) which reproduces the Figures 2b, 2d, 3b, 3d, 4b, 4d in the article - https://www.nature.com/articles/s41591-020-0883-7

# To reproduce Fig 2b and 2d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1, alpha_modifer = 1 and run for t = 45 d (for Fig 2b) and t = 350 days (for Fig 2d).
# Set alpha_modifier = 0 for the remaining 4 cases
# To reproduce Fig 3b, set Event_trigger_Fig3b = 1, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1 and run for t = 350 days.
# To reproduce Fig 3d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 1, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1 and run for t = 350 days.
# To reproduce Fig 4b, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 1, Event_trigger_Fig4d = 0, epsilon_modifer = 0 and run for t = 350 days.
# To reproduce Fig 4d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 1, epsilon_modifer = 0 and run for t = 350 days.

using OrdinaryDiffEq, ModelingToolkit, EasyModelAnalysis, SBML, SBMLToolkit, UnPack, Test
import Plots

@info "usings"

function sub_cont_ev(ev, rules)
    ModelingToolkit.SymbolicContinuousCallback(substitute(ev.eqs, rules),
                                               substitute(ev.affect, rules))
end
fn = "Scenario2/Giordano2020.xml"

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

@parameters t
D = Differential(t)
eqs2 = deepcopy(eqs)
append!(eqs2, D.(ps) .~ 0)

sys2 = ODESystem(eqs2, ModelingToolkit.get_iv(sys), states(sys), parameters(sys);
                 continuous_events = evs, defaults = defs, name = nameof(sys))
ssys = structural_simplify(sys2)
prob = ODEProblem(ssys, [], (0.0, 100.0))

sol = solve(prob, Tsit5())
@test sol.retcode == ReturnCode.Success
Plots.plot(sol)

sys = structural_simplify(odesys)
@unpack Infected, Exposed, Deceased, Recovered, Total_population, Susceptible = sys

pbounds = [
    alpha => [0.003, 0.006],
    epsilon => [1 / 6, 1 / 2],
    gamma => [0.1, 0.2],
    mu => [0.01, 0.02],
    beta => [0.7, 0.9],
]
create_sensitivity_plot(prob, 100.0, Deceased, pbounds; samples = 2000)

# Unit Test #1:
u0s = [I => 200/60e6, D => 20/60e6, A => 1/60e6, R => 2/60e6, T => 0, H => 0, E => 0; S => 1]
pars = [alpha => 0.570, beta => 0.011, delta = 0.011, gamma => 0.456, epsilon => 0.171, theta => 0.371,
        zeta => 0.125, eta => 0.125, mu => 0.017, nu => 0.027, tau => 0.01,
        lambda => 0.034, rho => 0.034,  kappa => 0.017,  xi => 0.017,  sigma => 0.017] # The resulting basic reproduction number is R0 = 2.38.

prob_test1 = remake(prob, u0 = u0s, p = pars)
# maybe we should allow get_max_t to take a term and have it generate an observable 
xmax, xmaxval = get_max_t(prob_test1, sum([I, D, A, R, T]))
