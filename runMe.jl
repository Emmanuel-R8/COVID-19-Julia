# The [Neherlab COVID-19](https://neherlab.org/covid19/) forecast model

# Runs with 4 processes (if 4 cores CPU) - Only for optimised libraries
#using Distributed
#TOTAL_PROCS = 4
#Distributed.addprocs(TOTAL_PROCS - Distributed.nprocs())

include("COVID-19-parameters.jl")
include("COVID-19-utils.jl")
include("COVID-19-model.jl")
include("COVID-19-run-model.jl")
include("COVID-19-data.jl")
PRINT_DEBUG = false

using BlackBoxOptim

# Global variable that holds the optimised disease parameters
DiseaseParameters = DISEASE_INIT

# If running for the first time, or no updates for a long time
# updateData()

#------------------------------------------------------------------------------------------------

#-- FOR DEBUGGING RUNS
# COUNTRY_LIST = [
#     ("France",                      :north),
#     ("Italy",                       :north),
#     ("Spain",                       :north),
#     ("United Kingdom",              :north)]

# The list is focused on Europe
COUNTRY_LIST = [
    ("Austria",                     :north),
    ("Belgium",                     :north),
    ("Bulgaria",                    :north),
    ("Denmark",                     :north),
    ("France",                      :north),
    ("Germany",                     :north),
    ("Greece",                      :north),
    ("Italy",                       :north),
    ("Netherlands",                 :north),
    ("Norway",                      :north),
    ("Poland",                      :north),
    ("Spain",                       :north),
    ("Sweden",                      :north),
    ("United_Kingdom",              :north),
#    ("United_States_of_America",     :north)
    ]


countryData = Dict( c => populateCountryData(c, h, useOptimised = false) for (c, h) in COUNTRY_LIST)

plt = plotVignette()

using FileIO, Plots, ImageIO
Plots.savefig(plt, "Vignette1.png")

plotCountriestoDisk("__beforeOptim")
saveParameters()


# Print only 3 decimals
using Printf
Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)

#-------------------------------------------------------------------------------------------------
#--
#-- First calibration
#--
updateEveryCountry(;maxtime = 20)
updateEpidemiologyOnce(maxtime = 120)

plt = plotVignette()
Plots.savefig(plt, "Vignette2.png")

plotCountriestoDisk("__first_calib")
saveParameters()


#-------------------------------------------------------------------------------------------------
#--
#-- Runs of all-country optim - disease parameters
#--
N_RUNS = 4
for run in 1:N_RUNS
    #-------------------------------------------------------------------------------------------------
    #--
    #-- Optimisition all countries at once

    #-- Build a vector of all the countries' parameters with start date constraints
    fullRange = empty([(0.0, 0.0)])
    for (c1, _) in COUNTRY_LIST
        countryRange = COUNTRY_RANGE
        countryRange[COUNTRY_PARAM_START] = approximateModelStartRange(c1)
        fullRange = vcat(fullRange, countryRange)
    end

    #-- Optimise
    best = best_candidate(bboptimize(sumCountryLossesCountries,
                                     SearchRange = fullRange;
                                     Method = :adaptive_de_rand_1_bin,
                                     MaxTime = 180,
                                     TargetFitness = 2.0,
                                     NThreads = Threads.nthreads(),
                                     TraceMode = :compact))

    #-- Store the optimised parameters
    for i in 1:COUNTRY_N
        country, _ = COUNTRY_LIST[i]

        country_start_index = (i - 1) * COUNTRY_N + 1
        country_final_index = (i - 1) * COUNTRY_N + COUNTRY_N

        global countryData[country][:params] = best[country_start_index:country_final_index]
    end

    #-------------------------------------------------------------------------------------------------
    # Update Epidemiology
    updateEpidemiologyOnce(maxtime = 60)

    #-------------------------------------------------------------------------------------------------
    # Outputs
    saveParameters()
    plotCountriestoDisk(repr(now()));

end

plotVignette()


#-------------------------------------------------------------------------------------------------
# Normal runs: Optimise country, then disease after each country
N_RUNS = 10
for run in 1:N_RUNS
    println(); println("OPTIMISING THE WORST HALF OF COUNTRIES BY AVERAGE ERRORS---------------------------")

    scores = allSingleLosses(sorted = true)
    @show scores
    for c1 in first(scores, COUNTRY_N ÷ 3)[:, 2]
        updateCountryOnce(c1; maxtime = 60)
    end

    # After having done all the countries, epidemiology again more seriously
    updateEpidemiologyOnce(; maxtime = 60)
    println()

    saveParameters()
    plotCountriestoDisk(repr(now()))
end

plotVignette()



#-------------------------------------------------------------------------------------------------
#--
#-- Global optimisation: all countries at once + disease parameters.
#--
fullRange = empty([(0.0, 0.0)])
for (c1, _) in COUNTRY_LIST
    countryRange = COUNTRY_RANGE
    countryRange[COUNTRY_PARAM_START] = approximateModelStartRange(c1)
    global fullRange = vcat(fullRange, countryRange)
end

fullRange = vcat(DISEASE_RANGE, fullRange)
result = bboptimize(fullEpidemyLoss,
                    SearchRange = fullRange;
                    Method = :adaptive_de_rand_1_bin,
                    MaxTime = 300,
                    TargetFitness = 2.0,
                    NThreads = Threads.nthreads(),
                    TraceMode = :compact)

best = best_candidate(result)
DiseaseParameters = best[1:DISEASE_N]
for i in 1:COUNTRY_LIST_N
    country, _ = COUNTRY_LIST[i]

    country_start_index = DISEASE_N + (i - 1) * COUNTRY_N + 1
    country_final_index = DISEASE_N + (i - 1) * COUNTRY_N + COUNTRY_N

    global countryData[country][:params] = best[country_start_index:country_final_index]
end

plotVignette()
saveParameters()
plotCountriestoDisk(repr(now()));
