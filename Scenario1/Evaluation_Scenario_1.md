# Evaluation Scenario 1

## Stratified SIR

```@example scenario1
using EasyModelAnalysis
using EasyModelAnalysis.ModelingToolkit: toparam
using EasyModelAnalysis.ModelingToolkit.Symbolics: FnType, variables
@parameters γ=1 / 14 R0=5
β = R0 * γ
@variables t S(t) I(t) R(t)
D = Differential(t)
n_stratify = 3
Ns = map(toparam, variables(:N, 1:n_stratify))
Ss = map(v -> v(t), variables(:S, 1:n_stratify, T = FnType))
Is = map(v -> v(t), variables(:I, 1:n_stratify, T = FnType))
Rs = map(v -> v(t), variables(:R, 1:n_stratify, T = FnType))
C = map(toparam, variables(:C, 1:n_stratify, 1:n_stratify))
k = 1000
contect_matrix = fill(0.33, (3, 3))
for i in 1:n_stratify
    nn = 2k
    defs[Ns[i]] = nn
    defs[Ss[i]] = nn - 1
    defs[Is[i]] = 1
    defs[Rs[i]] = 0
end
for i in eachindex(C)
    defs[C[i]] = contect_matrix[i]
end
eqs = [D.(Ss) .~ -β ./ Ns .* Ss .* (C * Is)
       D.(Is) .~ β ./ Ns .* Ss .* (C * Is) .- γ .* Is
       @. D(Rs) ~ γ * Is
       S ~ sum(Ss)
       I ~ sum(Is)
       R ~ sum(Rs)]
@named model = ODESystem(eqs; defaults = defs)
sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 300.0))
sol = solve(prob)
plot(sol, leg = :topright)
```

We can specify a non-uniform contact matrix

```@example scenario1
using Random
Random.seed!(123)
contect_matrix = 0.33 * rand(3, 3)
prob = ODEProblem(sys, [], (0, 300.0), vec(C .=> contect_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```
