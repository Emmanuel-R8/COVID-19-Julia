# Replicating the model presented by [NeherLab](https://neherlab.org/covid19/).

The content of this repo are explained [there](https://emmanuel-r8.github.io/2020/03/25/2020-03-25-forecasting-covid-19.html).

# Optimising the model parameters

## Overview
Parameters fall into 2 categories:

- Disease specific: those are independent from where the disease arises (at least in first approximation).
  Whether someone is in France or Germany, the level of care is essentially identical and
  fatality rates (proportion of death amongst ICU patients) is assumed identical.

- Country-specific: For example how many infected individuals at the beginning of the model,
  what do the mitigation profile looks like...

The optimisation of the model is assessed by the only certain measure: how well the model assesses the
number of deaths (using RMSE). Other published figures are not useful. For example, the number of confirmed
cases is only an unknown portion of the actual total number of cases. The optimisation process is
split into 2 consecutive steps:

  - Optimise disease parameters by modifying them to achieve the lowest possible loss across all
    countries simultaneously (i.e. summing all the losses of all the countries). During this, the
    country-specific parameters remain constant.

  - Optimise each individual country while keeping the disease parameters constant.

## Some results

Those results took less than 10 minutes on my old laptop (Thinkpad x260).

In summary, the results of the optimisation gives an r₀ lower than Neherlab's moderate value of 2.7. This is however
offset by a hug number of non-identified asymptomatic and symptomatic individuals (by a factor of about 20).  The number of days assumed in each state is more or less identical to the Neherlab values apart from the number of days in ICU being about 9 days instead of 14. In addition, mortality doubles for someone in a critical condition without medical attention.


               r₀     tₗ     tᵢ     tₕ     tᵤ     γₑ     γᵢ     γⱼ     γₖ     δₖ     δₗ     δᵤ     Asymp   Sympt   mv0    mv1    mv2    mv3    mv4    mv5    mv6    mv7    mv8    mv9

France        [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 11.673, 12.232, 1.000, 0.589, 0.119, 0.171, 0.375, 1.470, 1.489, 1.247, 1.334, 0.708],
Germany       [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 34.646, 24.006, 1.000, 0.671, 0.187, 0.252, 1.178, 0.774, 1.415, 0.279, 0.196, 1.079],
Italy         [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 36.232, 17.295, 1.000, 0.661, 0.471, 0.116, 0.504, 0.623, 0.363, 0.794, 0.256, 0.452],
Spain         [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 38.428, 11.369, 1.000, 0.803, 0.299, 0.114, 1.285, 1.091, 0.104, 0.619, 1.485, 1.061],
Switzerland   [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 35.814, 14.229, 1.000, 1.014, 0.207, 0.452, 1.130, 0.346, 1.342, 1.220, 0.619, 0.160],
UK            [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 20.750, 23.243, 1.000, 0.873, 0.196, 0.323, 0.863, 1.089, 0.506, 0.577, 0.336, 0.932],
USA           [2.235, 4.410, 4.925, 2.992, 8.821, 0.255, 1.000, 0.824, 1.980, 1.119, 1.925, 1.000, 34.370, 15.416, 1.000, 0.669, 0.257, 0.451, 1.097, 0.753, 0.395, 1.093, 0.901, 1.051]]


Legend:
  - r₀: average number of people infected by a symptomatic individual
  - tₗ     tᵢ     tₕ     tᵤ: average number of days as asymptomatic, symptomatic, hospitalised (severe) and in ICU (critical).
  - γₑ     γᵢ     γⱼ     γₖ: multiplier of the infectiousness. γᵢ is the infectiousness of a symptomatic individual and therefore set at 1. The other multiplier apply to asymptomatic, severe (when
out of hospital) and critical (when out of hospital)
  - δₖ     δₗ     δᵤ: multiplier of the fatality rate as compare to someone in ICU (δᵤ set at 1) for critical individuals in normal hospital beds or out of hospital.
  - Asymp   Sympt: The data provides a number of confirmed cases at the start of the model. Those are multipliers to get the number of asymptomatic and symptomatic individuals (only a portion of
them are actually confirmed).
  - mv...: Multiplier of r₀ reflecting the effectiveness of mitigation measures at weekly intervals.


# Source code

It is organised in 4 files, plus the `runMe.jl` that provides results.

- `COVID-19-model.jl`  sets a number of constants and the differential equations.

- `COVID-19-utils.jl` is a single function that is probably provided in the standard librabries...

- `COVID-19-run-model.jl`  sets the initial default value of the parameters to be optimised,
  calculates solutions and corresponding losses.

- `COVID-19-data.jl` updates data from the repos, populates dictionaries

The optimisation is done with the `BlackBoxOptim.jl`. It is just magic...
