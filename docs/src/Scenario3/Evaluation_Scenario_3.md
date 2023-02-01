# Evaluation Scenario 3

```@example evalscenario3
using EasyModelAnalysis, LinearAlgebra
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
tspan = (0.0, 40.0)
u0 = [990, 10, 0]
p = [0.05*10/1000, 0.25]
prob = ODEProblem(sys, u0, tspan, p);
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
fitparams = global_datafit(prob, [β => [0.03, 0.15], c => [9.0, 13.0], γ => [0.05, 0.5]],
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

#### Application to Real Data from TA1

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
tspan = (0.0, 40.0)
u0 = [990, 10, 0, 0, 0]
p = [0.01*10/1000, 0.25, 0.1, 0.1, 0.1, 0.1]
sirhd_prob = ODEProblem(sirhd_sys, u0, tspan, p)
sirhd_sol = solve(sirhd_prob)
plot(sirhd_sol)
```

Question 2 involves doing the same analysis as question one but on the SIR model with hopsitalizations and deaths included.
To establish unit tests, we first showcase building the model and solving inverse problems using the ModelingToolkit version
of the model.

The inverse problem solving is done via the same functionality as before.

```@example evalscenario3
fitparams2 = global_datafit(sirhd_prob, [β => [0.03, 0.15], c => [9.0, 13.0], γ => [0.05, 0.5]],
                            t_train, data_train)
```

Notice that this fit is not as good. That is to be expected because it's fitting the SIRHD model on the
SIR model's output data. Thus we should expect that it also does not forecast entirely correctly.

```@example evalscenario3
sirhd_prob2 = remake(sirhd_prob, p = fitparams2)
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
norm(solve(_prob2, saveat = t_test)[R] - data_test[3][2])
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

This expands the previous SIRHD-Vax model to stratify by age.
```@example evalscenario3
sirhd_vax_age16 = read_json_acset(LabelledPetriNet,"sirhd_vax_age16.json")
```

Question 4 is the same analysis as questions 1, 2, and 3 on a model with age-stratification added. In order to build unit tests for
the analysis and functionality, we started by building the model with vaccine by hand, awaiting a swap to the version from
TA2.

#### Data

* Previous data that is age stratified is cases, and hospitalizations
* 10 stratifications, by 10 years each
* Underreporting over time?
* Data for assumption on recovery rate with respect to age
* Aggregated contact matrix for beta over age, from Scenario 1

## Question 5: Add Reinfection

This expands the previous SIRHD model to add reinfection, and restratifies with vaccination and age.
```@example evalscenario3
sirhd_renew = read_json_acset(LabelledPetriNet,"sirhd_renew.json")
sirhd_renew_vax = read_json_acset(LabelledPetriNet,"sirhd_renew_vax.json")
sirhd_renew_vax_age16 = read_json_acset(LabelledPetriNet,"sirhd_renew_vax_age16.json")
```

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
