# Evaluation Scenario 1

## Stratified SIR Example

```@example scenario1
using EasyModelAnalysis
using EasyModelAnalysis.ModelingToolkit: toparam
using EasyModelAnalysis.ModelingToolkit.Symbolics: FnType, variables
@parameters β=0.1 γ=0.2
@variables t S(t) I(t) R(t)
D = Differential(t)
n_stratify = 3
Ns = map(v -> toparam(v), variables(:N, 1:n_stratify))
Ss = map(v -> v(t), variables(:S, 1:n_stratify, T = FnType))
Is = map(v -> v(t), variables(:I, 1:n_stratify, T = FnType))
Rs = map(v -> v(t), variables(:R, 1:n_stratify, T = FnType))
using Random
Random.seed!(123)
defs = Dict(x => rand() for x in Iterators.flatten((Ns, Ss, Is, Rs)))
eqs = [map((s, N) -> D(s) ~ -β / N * I * s, Ss, Ns);
       map((i, N) -> D(i) ~ β / N * i * S - γ * i, Is, Ns);
       map((r, i) -> D(r) ~ γ * i, Rs, Is);
       S ~ sum(Ss)
       I ~ sum(Is)
       R ~ sum(Rs)]
@named model = ODESystem(eqs; defaults = defs)
sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 50.0))
sol = solve(prob)
plot(sol, leg = :topright)
```
