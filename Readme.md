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

In summary, the results of the optimisation gives an r₀ slightly lower than Neherlab's moderate value of 2.7. This is however offset by a hug number of non-identified asymptomatic and symptomatic individuals (by a factor of about 40).  The number of days assumed in each state is more or less identical to the Neherlab values.



| Country     | r₀    | tₗ    | tᵢ    | tₕ    | tᵤ     | γₑ    | γᵢ    | γⱼ    | γₖ    | δₖ    | δₗ    | δᵤ    | Asymp  | Sympt  | mv0   | mv1   | mv2   | mv3   | mv4   | mv5   | mv6   | mv7   | mv8   | mv9   |
|-------------|-------|-------|-------|-------|--------|-------|-------|-------|-------|-------|-------|-------|--------|--------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| France      | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.073 | 2.048 | 0.487 | 1.175 | 1.216 | 0.833 | 40.808 | 44.067 | 0.992 | 1.020 | 0.197 | 0.549 | 0.963 | 0.387 | 0.539 | 0.899 | 0.560 | 1.092 |
| Germany     | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.256 | 2.048 | 0.487 | 1.175 | 1.216 | 0.822 | 46.740 | 43.408 | 1.042 | 0.896 | 0.172 | 0.552 | 0.904 | 0.375 | 0.509 | 0.974 | 0.523 | 1.100 |
| Italy       | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.167 | 2.048 | 0.487 | 1.175 | 1.216 | 0.848 | 42.731 | 42.805 | 0.996 | 0.965 | 0.173 | 0.597 | 1.063 | 0.374 | 0.464 | 0.850 | 0.616 | 1.234 |
| Spain       | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.192 | 2.048 | 0.487 | 1.175 | 1.216 | 0.868 | 47.949 | 40.095 | 0.917 | 0.957 | 0.184 | 0.519 | 0.938 | 0.383 | 0.474 | 0.957 | 0.582 | 1.078 |
| Switzerland | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.237 | 2.048 | 0.487 | 1.175 | 1.216 | 0.857 | 46.312 | 42.157 | 1.080 | 0.890 | 0.193 | 0.554 | 1.083 | 0.377 | 0.555 | 0.857 | 0.625 | 1.048 |
| UK          | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.277 | 2.048 | 0.487 | 1.175 | 1.216 | 0.819 | 45.160 | 41.339 | 0.961 | 0.921 | 0.167 | 0.515 | 0.994 | 0.385 | 0.495 | 0.871 | 0.551 | 1.149 |
| USA         | 2.493 | 5.795 | 5.597 | 3.179 | 13.933 | 0.231 | 1.201 | 2.048 | 0.487 | 1.175 | 1.216 | 0.807 | 43.230 | 39.758 | 1.074 | 0.965 | 0.178 | 0.517 | 0.953 | 0.374 | 0.460 | 0.985 | 0.625 | 1.224 |


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
