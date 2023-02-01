# Evaluation Scenario 3

```@example evalscenario3
using EasyModelAnalysis, LinearAlgebra, CSV
using Catlab, AlgebraicPetri
using Catlab.CategoricalAlgebra
```

## Question 1

## Setup Model

>  1. Begin with a basic SIR model without vital dynamics. Calibrate the model parameters using data on cases during the ‘training period’. Then evaluate the model during the out-of-sample ‘test period’.

To get started with the code, we implemented the basic SIR without vital dynamics directly in ModelingToolkit.jl. This is a version
that was written by an epidemiologist at Microsoft Pandemic, Simon Frost, who has become a fan of the TA3 automated simulation tools
and wrote an entire repository of tutorials for this software. It is found at https://github.com/epirecipes/sir-julia.

In there is an SIR without vital dynamics which we took in full.

```@example evalscenario3
sir = read_json_acset(LabelledPetriNet, "sir.json")
sys = ODESystem(sir)
sys = complete(sys)
@unpack S, I, R, inf, rec = sys
@parameters N = 1
param_sub = [
    inf => inf / N
]
sys = substitute(sys, param_sub)
defs = ModelingToolkit.defaults(sys)
defs[S] = 990
defs[I] = 10
defs[R] = 0.0
defs[N] = sum(x->defs[x], (S, I, R))
defs[inf] = 0.5
defs[rec] = 0.25
tspan = (0.0, 40.0)
prob = ODEProblem(sys, [], tspan);
sol = solve(prob);
```

```@example evalscenario3
plot(sol)
```

### Perform Model Calibration

#### Model Calibration Unit Test

As a unit test of the model calibration tools, we generated data at the default parameters, then ran the global optimization,
to see how well the parameters were recovered.

```@example evalscenario3
dataset = solve(prob, saveat = 0.1)
t_train = dataset.t[1:201]
t_test = dataset.t[202:end]
data_train = [S => dataset[S][1:201], I => dataset[I][1:201], R => dataset[R][1:201]]

data_test = [S => dataset[S][202:end], I => dataset[I][202:end], R => dataset[R][202:end]]
```

```@example evalscenario3
fitparams = global_datafit(prob, [inf => [0.2, 2.0], rec => [0.05, 0.5]],
                           t_train, data_train)
```

This then gives the forecasts in the test data:

```@example evalscenario3
_prob = remake(prob, p = fitparams)
sol = solve(_prob, saveat = t_test);
plot(sol, idxs = S)
plot!(t_test, data_test[1][2])
```

```@example evalscenario3
plot(sol, idxs = I)
plot!(t_test, data_test[2][2])
```

```@example evalscenario3
plot(sol, idxs = R)
plot!(t_test, data_test[3][2])
```

This looks very good and matches the original data, confirming that the inverse problem functionality is functional.

Now we train on data from June 1 2021 to September 30 2021.


#### Application to Real Data from TA1

```@example evalscenario3
using CSV, DataFrames, Downloads

# Infectious/Recovered day by day:
url = "https://raw.githubusercontent.com/DARPA-ASKEM/program-milestones/4ee24c6a268998148fc635b1e35088f4ebbbb7ce/6-month-milestone/evaluation/scenario_3/ta_4/usa-IRDVHN_age.csv"
file = CSV.File(Downloads.download(url))
df_raw = DataFrame(file)

start_train = 171
stop_train = 171+121
start_test = 171+122
stop_test = 171+122+92

df_train = df_raw[start_train:stop_train, :]
df_test = df_raw[start_test:stop_test, :]

t_train = collect(0:(size(df_train, 1)-1))
t_test = collect(0:(size(df_test, 1)-1))

N_total = 334998398 # assumed to be constant from (https://github.com/DARPA-ASKEM/program-milestones/blob/main/6-month-milestone/evaluation/scenario_3/ta_1/usa-2021-population-age-stratified.csv)
#S = N_total - R - I
data_train = [S => N_total .- df_train.I .-  df_train.R, I => df_train.I, R => df_train.R]
data_test = [S => N_total .- df_test.I .- df_test.R, I => df_test.I, R => df_test.R]

u0s = [S => N_total - df_train.I[1] - df_train.R[1], I => df_train.I[1], R => df_train.R[1]]
_prob = remake(prob, u0 = u0s, tspan = (t_train[1], t_train[end]), p = [N => N_total])

fitparams = global_datafit(_prob, [inf => [0, 1.0], rec => [0.0, 1.0]], t_train, data_train)
```

```@example evalscenario3
# Plot training fit
_prob_train = remake(_prob, p = fitparams)
sol = solve(_prob_train, saveat = t_train);

cs = Plots.distinguishable_colors(10)[end-5:end]
plot(sol, idxs = S, color = cs[1])
plot!(t_train, data_train[1][2], lab = "S_train", color = cs[2])

plot!(sol, idxs = I, color = cs[3])
plot!(t_train, data_train[2][2], lab = "I_train", color = cs[4])

plot!(sol, idxs = R, color = cs[5])
p = plot!(t_train, data_train[3][2], lab = "R_train", color = cs[6], dpi=300)
```
```@example evalscenario3
savefig(p, "train_fit_S3_Q1.png")
```

# Plot test fit
```@example evalscenario3
u0s = [S => N_total - df_test.I[1] - df_test.R[1], I => df_test.I[1], R => df_test.R[1]]
_prob_test = remake(_prob, p = fitparams, u0=u0s, tspan = (t_test[1], t_test[end]))
sol = solve(_prob_test, saveat = t_test);
plot(sol, idxs = S, color = cs[1])
plot!(t_test, data_test[1][2], lab = "S_test", color = cs[2])

plot!(sol, idxs = I, color = cs[3])
plot!(t_test, data_test[2][2], lab = "I_test", color = cs[4])

plot!(sol, idxs = R, color = cs[5])
p = plot!(t_test, data_test[3][2], lab = "R_test", color = cs[6], dpi=300)
```
```@example evalscenario3
savefig(p, "test_fit_S3_Q1.png")
```


* Expect time series data on I + R
* Start with an assumption on the recovery
* Possible additoinal: alternative measure for recovery rate
* Modeling assumption: use total infections from 2 weeks ago as R0, determine I0 and S0 from that
* Need time series for total population of US over time

## Question 2: Add Hospitalizations and Deaths

This expands the original SIR model to explore a model space comprising SIRD, SIRH, and SIRHD.
```@example evalscenario3
sird = read_json_acset(LabelledPetriNet,"sird.json")
sirh = read_json_acset(LabelledPetriNet,"sirh.json")
sirhd = read_json_acset(LabelledPetriNet,"sirhd.json")
sirhd_sys = ODESystem(sirhd)
sirhd_sys = complete(sirhd_sys)
@unpack S, I, R, H, D, inf, rec, ideath, death, hosp, hrec = sirhd_sys
@parameters N = 1
param_sub = [
    inf => inf / N
]
sirhd_sys = substitute(sirhd_sys, param_sub)
defs = ModelingToolkit.defaults(sirhd_sys)
defs[S] = N_total - 10
defs[I] = 10
defs[H] = 0
defs[D] = 0
defs[R] = 0.0
defs[N] = N_total
defs[inf] = 0.5
defs[rec] = 0.25
defs[ideath] = 0.25
defs[death] = 0.25
defs[hosp] = 0.25
defs[hrec] = 0.25
tspan = (0.0, 40.0)
sirhd_prob = ODEProblem(sirhd_sys, [], tspan)
sirhd_sol = solve(sirhd_prob)
plot(sirhd_sol)
```

Question 2 involves doing the same analysis as question one but on the SIR model with hopsitalizations and deaths included.
To establish unit tests, we first showcase building the model and solving inverse problems using the ModelingToolkit version
of the model.

The inverse problem solving is done via the same functionality as before.

```@example evalscenario3
param_bounds = [inf, rec, ideath, death, hosp, hrec] .=> ([0.01, 10.0],)
_prob = remake(sirhd_prob, u0 = u0s, tspan = (t_train[1], t_train[end]), p = [N => N_total])
fitparams2 = global_datafit(_prob, param_bounds, t_train, data_train, maxiters = 200_000)
```

Notice that this fit is not as good. That is to be expected because it's fitting the SIRHD model on the
SIR model's output data. Thus we should expect that it also does not forecast entirely correctly.

```@example evalscenario3
sirhd_prob2 = remake(_prob, p = fitparams2)
sol = solve(sirhd_prob2, saveat = t_test);
plot(sol, idxs = S)
plot!(t_test, data_test[1][2])
```

```@example evalscenario3
plot(sol, idxs = I)
plot!(t_test, data_test[2][2])
```

```@example evalscenario3
plot(sol, idxs = R)
plot!(t_test, data_test[3][2])
```

This checks out.

#### Data Ask

```@example evalscenario3
data_train = [
S => N_total .- df_train.I .-  df_train.R .- df_train.D .- df_train.H,
I => df_train.I, R => df_train.R, H => df_train.H, D => df_train.D
]
data_test = [
S => N_total .- df_test.I .-  df_test.R .- df_test.D .- df_test.H,
I => df_test.I, R => df_test.R, H => df_test.H, D => df_test.D
]

u0s = [
S => N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1],
I => df_train.I[1], R => df_train.R[1], H => df_train.H[1], D => df_train.D[1]
]
_prob2 = remake(sirhd_prob2, u0 = u0s, tspan = (t_train[1], t_train[end]), p = [N => N_total])

param_bounds = [
    inf => [0.0, 70]
    rec => [0.0, 5.0]
    death => [0.0, 5.0]
    ideath => [0.0, 5.0]
    hosp => [0.0, 10.0]
    hrec => [0.0, 10.0]
]
fitparams2 = global_datafit(_prob2, param_bounds, t_train, data_train, maxiters = 200_000)
```
```@example evalscenario3
# Plot training fit
_prob2_train = remake(_prob2, p = fitparams2)
sol = solve(_prob2_train, saveat = t_train);

plot(sol, idxs = S, color = cs[1])
plot!(t_train, data_train[1][2], lab = "S_train", color = cs[2])

plot!(sol, idxs = I, color = cs[3])
plot!(t_train, data_train[2][2], lab = "I_train", color = cs[4])

plot!(sol, idxs = R, color = cs[5])
p = plot!(t_train, data_train[3][2], lab = "R_train", color = cs[6], dpi=300)
```
```@example evalscenario3
savefig(p, "train_fit_S3_Q2.png")
```
```@example evalscenario3
u0s = [
S => N_total - df_test.I[1] - df_test.R[1] - df_test.H[1] - df_test.D[1],
I => df_test.I[1], R => df_test.R[1], H => df_test.H[1], D => df_test.D[1]
]
_prob2 = remake(_prob2, p = fitparams, u0=u0s, tspan = (t_test[1], t_test[end]))
sol = solve(_prob2, saveat = t_test);

plot(sol, idxs = S, color = cs[1])
plot!(t_test, data_test[1][2], lab = "S_test", color = cs[2])

plot!(sol, idxs = I, color = cs[3])
plot!(t_test, data_test[2][2], lab = "I_test", color = cs[4])

plot!(sol, idxs = R, color = cs[5])
p = plot!(t_test, data_test[3][2], lab = "R_test", color = cs[6], dpi=300)
```
```@example evalscenario3
savefig(p, "test_fit_S3_Q2.png")
```

* Daily time series on number of patients admitted to the hospital all US
* time series for mortality
* 10 gig file on whether hospitalized or not => percentage for the difference in parameters
    * Plot the percentage over time by month, see if a constant assumption is okay or not,
    * If not, need to use the time series
* Any factor for underreporting estimate? Wastewater time series

### Evaluate Model Forecasts

In order to evaluate the model forecasts, we developed a functional which does the forecasting part with multiple models
and puts a score on the forecast result. This score is calculated using the L2 norm. It was added to the EasyModelAnalysis.jl
library in https://github.com/SciML/EasyModelAnalysis.jl/pull/129 as part of the evaluation on day 1.

```@example evalscenario3
norm(solve(_prob, saveat = t_test)[S] - data_test[1][2]) +
norm(solve(_prob, saveat = t_test)[I] - data_test[2][2]) +
norm(solve(_prob, saveat = t_test)[R] - data_test[3][2])
```

```@example evalscenario3
norm(solve(_prob2, saveat = t_test)[S] - data_test[1][2]) +
norm(solve(_prob2, saveat = t_test)[I] - data_test[2][2]) +
norm(solve(_prob2, saveat = t_test)[R] - data_test[3][2]) +
norm(solve(_prob2, saveat = t_test)[H] - data_test[4][2]) +
norm(solve(_prob2, saveat = t_test)[D] - data_test[5][2])
```

## Question 3: Add Vaccinations

This expands the previous SIRHD model to add vaccination.
```@example evalscenario3
sirhd_vax = read_json_acset(LabelledPetriNet, "sirhd_vax.json")
sirhd_vax_sys = structural_simplify(ODESystem(sirhd_vax))
```

Question 3 is the same analysis as questions 1 and 2 done on a model with vaccination added. In order to build unit tests for
the analysis and functionality, we started by building the model with vaccine by hand, awaiting a swap to the version from
TA2.

The unit test analysis code is as follows:

```@example evalscenario3

```

#### Data Asks

* Time series of vaccinations
* Hospitalization rate difference due to vaccination?
* Recovery rate difference due to vaccination?
* Mortality rate difference due to vaccination? Hospitalized and not hospitalized

## Question 4: Age-Stratified Model

Question 4 is the same analysis as questions 1, 2, and 3 on a model with age-stratification added. In order to build unit tests for
the analysis and functionality, we started by building the model with vaccine by hand, awaiting a swap to the version from
TA2.
```julia
using ModelingToolkit.Symbolics: variable
using ModelingToolkit: toparam
sirhd_vax_age16 = read_json_acset(LabelledPetriNet,"sirhd_vax_age16.json")
sirhd_vax_age16  = ODESystem(sirhd_vax_age16)
sirhd_vax_age16  = complete(sirhd_vax_age16 )
_names = string.(ModelingToolkit.getname.(states(sirhd_vax_age16)))

sts_names = Symbol.(getindex.(_names, 6), :_, getindex.(_names, 11))
sts = map(n->variable(n, T = SymbolicUtils.FnType)(t), sts_names)

__names = split.(string.(parameters(sirhd_vax_age16)), "\"")
ps_names = Symbol.(getindex.(split.(getindex.(__names, 3), "\\"), 1), :_,
                   getindex.(split.(getindex.(__names, 5), "\\"), 1))
ps = map(n->toparam(variable(n)), ps_names)
subs = [
    parameters(sirhd_vax_age16) .=> ps
    states(sirhd_vax_age16) .=> sts
]
sirhd_vax_age16 = substitute(sirhd_vax_age16, subs)
@unpack S_U, I_U, R_U, H_U, D_U, S_V, I_V, R_V, H_V, D_V = sirhd_vax_sys
@unpack id_vax, inf_infuu, inf_infuv, hosp_id, ideath_id, rec_id, hrec_id, death_id, inf_infvu, inf_infvv = sirhd_vax_sys

@parameters N = 1
param_sub = [
    inf => inf / N
]
sirhd_sys = substitute(sirhd_vax_age16, param_sub)
defs = ModelingToolkit.defaults(sirhd_vax_age16)

# @unpack S, I, R, H, D, inf, rec, ideath, death, hosp, hrec = sirhd_sys # How to do this for stratified model?

defs[S] = N_total - 10
defs[I] = 10
defs[H] = 0
defs[D] = 0
defs[R] = 0.0
defs[N] = N_total
defs[inf] = 0.5
defs[rec] = 0.25
defs[ideath] = 0.25
defs[death] = 0.25
defs[hosp] = 0.25
defs[hrec] = 0.25
tspan = (0.0, 40.0)
sirhd_prob = ODEProblem(sirhd_sys, [], tspan)
sirhd_sol = solve(sirhd_prob)
plot(sirhd_sol)
```

#### Data
```@example evalscenario3
age_groups = ["_0-9", "_10-19", "_20-29", "_30-39", "_40-49", "_50-59",
              "_60-69", "_70-79", "_80-89", "_90-99", "_100+"]
data_train = Dict()
data_test = Dict()
for ag in age_groups
    data_train[ag[2:end]] = [(S+Sv) => N_total .- df_train.I .-  df_train.R .- df_train.H .-  df_train.D,
                (I+Iv) => df_train.I,
                (R+Rv) => df_train.R,
                H => df_train[:, "H_unvac$(ag)"],
                Hv => df_train[:, "H_vac$(ag)"],
                D => df_train[:, "D_unvac$(ag)"],
                Dv => df_train[:, "D_vac$(ag)"]
                ] # Todo: potentially rename symbols depending on what the model looks like
    data_test[ag[2:end]] = [(S+Sv) => N_total .- df_test.I .-  df_test.R .- df_test.H .-  df_test.D,
                (I+Iv) => df_test.I,
                (R+Rv) => df_test.R,
                H => df_test[:, "H_unvac$(ag)"],
                Hv => df_test[:, "H_vac$(ag)"],
                D => df_test[:, "D_unvac$(ag)"],
                Dv => df_test[:, "D_vac$(ag)"]
                ] # Todo: potentially rename symbols depending on what the model looks like
end  # There is also data for I, but not for the oldest 3 age groups. So far we ignore this.
```

* Previous data that is age stratified is cases, and hospitalizations
* 10 stratifications, by 10 years each
* Underreporting over time?
* Data for assumption on recovery rate with respect to age
* Aggregated contact matrix for beta over age, from Scenario 1

## Question 5: Add Reinfection

Question 5 is the same analysis as questions 1, 2, 3, and 4 on a model with reinfection added. In order to build unit tests for
the analysis and functionality, we started by building the model with vaccine by hand, awaiting a swap to the version from
TA2.

#### Data Asks

* Change in hospitalization for people who are reinfected
* State of new york, people who reinfected?
* Median time to reinfection
* It may require R -> S ===> R -> S2
* Maybe model recovered as vaccinated S?

## Question 6: New Data

Question 6 is currently awaiting data from TA3

## Question 7: Analysis

For each model, summarize your conclusions about the following:

 1. Do parameters fit from data seem reasonable and fall within typical ranges you might see in the broader literature? Provide references to support your conclusions.
 2. Describe how well the fitted model compares against historical data, both for the ‘training’ and ‘test’ periods.

### Answer

Question 7 is currently awaiting data from TA3
