# this script requires https://github.com/SciML/SBMLToolkit.jl/pull/108
# This is a model created on COPASI 4.27 (Build 217) which reproduces the Figures 2b, 2d, 3b, 3d, 4b, 4d in the article - https://www.nature.com/articles/s41591-020-0883-7

# To reproduce Fig 2b and 2d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1, alpha_modifer = 1 and run for t = 45 d (for Fig 2b) and t = 350 days (for Fig 2d).
# Set alpha_modifier = 0 for the remaining 4 cases
# To reproduce Fig 3b, set Event_trigger_Fig3b = 1, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1 and run for t = 350 days.
# To reproduce Fig 3d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 1, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1 and run for t = 350 days.
# To reproduce Fig 4b, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 1, Event_trigger_Fig4d = 0, epsilon_modifer = 0 and run for t = 350 days.
# To reproduce Fig 4d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 1, epsilon_modifer = 0 and run for t = 350 days.

# The total population is partitioned into eight stages of dis- ease: 
# S, susceptible (uninfected);
# I, infected (asymptomatic or pauci-symptomatic infected, undetected); -> ND AS
# D, diagnosed (asymp- tomatic infected, detected); D AS
# A, ailing (symptomatic infected, undetected); 
# R, recognized (symptomatic infected, detected); 
# T, threatened (infected with life-threatening symptoms, detected); 
# H, healed (recovered); E, extinct (dead). 

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

@unpack Infected, Healed, Extinct, Diagnosed, Ailing, Recognized, Susceptible, Threatened = sys
@unpack alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta, theta, tau = sys
ps = [alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta]

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

# Unit Test #1:
ITALY_POPULATION = 60e6
u0s = [
    Infected => 200 / ITALY_POPULATION,
    Diagnosed => 20 / ITALY_POPULATION,
    Ailing => 1 / ITALY_POPULATION,
    Recognized => 2 / ITALY_POPULATION,
    Threatened => 0,
    Healed => 0,
    Extinct => 0,
]
push!(u0s, Susceptible => 1 - sum(last.(u0s)))

# The resulting basic reproduction number is R0 = 2.38.
pars = [alpha => 0.570, beta => 0.011, delta => 0.011, gamma => 0.456, epsilon => 0.171,
    theta => 0.371,
    zeta => 0.125, eta => 0.125, mu => 0.017, nu => 0.027, tau => 0.01,
    lambda => 0.034, rho => 0.034, kappa => 0.017, xi => 0.017, sigma => 0.017]

prob_test1 = remake(prob, u0 = u0s, p = pars)
solt1 = solve(prob_test1, Tsit5(); saveat = 0:100)
# solt1 = solve(prob_test1, Tsit5(); saveat = 0:45)
og_states = states(sys)[1:8]
idart = [Infected, Diagnosed, Ailing, Recognized, Threatened]
plot(solt1; idxs = Infected)
plot(solt1; idxs = Diagnosed)
plot(solt1; idxs = idart)
@test solt1[Infected + Healed] == solt1[Infected] + solt1[Healed]
plot(solt1.t, solt1[sum(idart)] * ITALY_POPULATION; label = "IDART absolute")
plot(solt1.t, solt1[sum(idart)]; label = "IDART percent")

xmax, xmaxval = get_max_t(prob_test1, sum(idart) * ITALY_POPULATION)

@test isapprox(xmax, 47; atol = 5)
@test isapprox(xmaxval, 0.6)

# Unit test #2:
# To reproduce Fig 2b and 2d, set Event_trigger_Fig3b = 0, Event_trigger_Fig3d = 0, Event_trigger_Fig4b = 0, Event_trigger_Fig4d = 0, epsilon_modifer = 1, alpha_modifer = 1 and run for t = 45 d
# -> These are already the settings in the Model
sysv = eval(quote
                var"##iv#608" = (@variables(t))[1]
                var"##sts#609" = (collect)(@variables(Infected(t), Healed(t), Extinct(t),
                                                      Diagnosed(t), Ailing(t),
                                                      Recognized(t), Susceptible(t),
                                                      Threatened(t), Vaccinated(t),
                                                      alpha(t), epsilon(t), gamma(t),
                                                      beta(t), delta(t), mu(t), nu(t),
                                                      lambda(t), rho(t), kappa(t), xi(t),
                                                      sigma(t), zeta(t), eta(t)))
                var"##ps#610" = (collect)(@parameters(ModelValue_21, epsilon_modifier, tau,
                                                      theta, ModelValue_19, ModelValue_20,
                                                      Event_trigger_Fig4b,
                                                      Event_trigger_Fig3d, ModelValue_16,
                                                      Event_trigger_Fig3b, alpha_modifier,
                                                      Event_trigger_Fig4d, ModelValue_17,
                                                      ModelValue_18, Italy, tau1, phi))
                var"##eqs#611" = [(~)((Differential(t))(Infected),
                                      (+)((*)((+)((/)((*)(Infected, alpha), Italy),
                                                  (/)((*)(Diagnosed, beta), Italy),
                                                  (/)((*)(Ailing, gamma), Italy),
                                                  (/)((*)(Recognized, delta), Italy)),
                                              Susceptible), (*)(-1 // 1, Infected, epsilon),
                                          (*)(-1 // 1, Infected, lambda),
                                          (*)(-1 // 1, Infected, zeta)))
                                  (~)((Differential(t))(Healed),
                                      (+)((*)(Ailing, kappa), (*)(Diagnosed, rho),
                                          (*)(Infected, lambda), (*)(Recognized, xi),
                                          (*)(Threatened, sigma)))
                                  (~)((Differential(t))(Extinct),
                                      (*)(tau, Threatened) + tau1 * Recognized)
                                  (~)((Differential(t))(Diagnosed),
                                      (+)((*)(Infected, epsilon),
                                          (*)(-1 // 1, Diagnosed, eta),
                                          (*)(-1 // 1, Diagnosed, rho)))
                                  (~)((Differential(t))(Ailing),
                                      (+)((*)(Infected, zeta), (*)(-1 // 1, theta, Ailing),
                                          (*)(-1 // 1, Ailing, kappa),
                                          (*)(-1 // 1, Ailing, mu)))
                                  (~)((Differential(t))(Recognized),
                                      (+)((*)(theta, Ailing), (*)(Diagnosed, eta),
                                          (*)(-1 // 1, Recognized, nu),
                                          (*)(-1 // 1, Recognized, xi)))
                                  (~)((Differential(t))(Susceptible),
                                      (*)((+)((/)((*)(-1, Infected, alpha), Italy),
                                              (/)((*)(-1, Diagnosed, beta), Italy),
                                              (/)((*)(-1, Ailing, gamma), Italy),
                                              (/)((*)(-1, Recognized, delta), Italy)) -
                                          phi * Susceptible,
                                          Susceptible))
                                  (~)((Differential(t))(Threatened),
                                      (+)((*)(Ailing, mu), (*)(Recognized, nu),
                                          (*)(-1 // 1, tau, Threatened),
                                          (*)(-1 // 1, Threatened, sigma)))
                                  Differential(t)(Vaccinated) ~ phi * Susceptible
                                  (~)((Differential(t))(alpha), -0.0);
                                  (~)((Differential(t))(epsilon), -0.0);
                                  (~)((Differential(t))(gamma), -0.0);
                                  (~)((Differential(t))(beta), -0.0);
                                  (~)((Differential(t))(delta), -0.0);
                                  (~)((Differential(t))(mu), -0.0);
                                  (~)((Differential(t))(nu), -0.0);
                                  (~)((Differential(t))(lambda), -0.0);
                                  (~)((Differential(t))(rho), -0.0);
                                  (~)((Differential(t))(kappa), -0.0);
                                  (~)((Differential(t))(xi), -0.0);
                                  (~)((Differential(t))(sigma), -0.0);
                                  (~)((Differential(t))(zeta), -0.0);
                                  (~)((Differential(t))(eta), -0.0)]
                var"##defs#612" = (Dict)((Pair)(delta, 0.011), (Pair)(xi, 0.017),
                                         (Pair)(Diagnosed, 3.33333333e-7),
                                         (Pair)(Event_trigger_Fig3b, 0.0),
                                         (Pair)(Extinct, 0.0), (Pair)(kappa, 0.017),
                                         (Pair)(zeta, 0.125), (Pair)(eta, 0.125),
                                         (Pair)(nu, 0.027), (Pair)(Healed, 0.0),
                                         (Pair)(Infected, 3.33333333e-6),
                                         (Pair)(ModelValue_16, 0.0),
                                         (Pair)(alpha_modifier, 1.0), (Pair)(Italy, 1.0),
                                         (Pair)(Event_trigger_Fig3d, 0.0),
                                         (Pair)(ModelValue_20, 1.0), (Pair)(sigma, 0.017),
                                         (Pair)(Threatened, 0.0), (Pair)(lambda, 0.034),
                                         (Pair)(alpha, 0.57),
                                         (Pair)(Event_trigger_Fig4b, 0.0),
                                         (Pair)(ModelValue_17, 0.0),
                                         (Pair)(Event_trigger_Fig4d, 0.0),
                                         (Pair)(Susceptible, 0.9999963),
                                         (Pair)(beta, 0.011),
                                         (Pair)(Recognized, 3.33333333e-8),
                                         (Pair)(rho, 0.034), (Pair)(mu, 0.017),
                                         (Pair)(epsilon, 0.171),
                                         (Pair)(Ailing, 1.66666666e-8),
                                         (Pair)(gamma, 0.456), (Pair)(ModelValue_19, 0.0),
                                         (Pair)(ModelValue_21, 1.0), (Pair)(theta, 0.371),
                                         (Pair)(epsilon_modifier, 1.0), (Pair)(tau, 0.01),
                                         (Pair)(ModelValue_18, 0.0),
                                         Vaccinated => 0,
                                         tau1 => 0.0200,
                                         phi => 0.0025)
                var"##iv#613" = (@variables(t))[1]
                (ODESystem)(var"##eqs#611", var"##iv#613", var"##sts#609", var"##ps#610";
                            defaults = var"##defs#612", name = Symbol("##SBML#530"),
                            checks = false)
            end)
# todo set the event flags
# todo validate the new params 
prob = ODEProblem(sysv, [], (0, 100))
sol = solve(prob, Tsit5())
plot(sol)
plot(sol, idxs=[og_states; Vaccinated])

# @unpack Infected, Extinct = sys

# pbounds = [
#     alpha => [0.003, 0.006],
#     epsilon => [1 / 6, 1 / 2],
#     gamma => [0.1, 0.2],
#     mu => [0.01, 0.02],
#     beta => [0.7, 0.9],
# ]
# create_sensitivity_plot(prob, 100.0, Extinct, pbounds; samples = 20)
