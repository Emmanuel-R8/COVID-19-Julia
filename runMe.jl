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

using BlackBoxOptim
using PyPlot
pyplot()

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
#     ("United Kingdom",              :north),
#     ("United States of America",    :north)]

COUNTRY_LIST = [
    ("Austria",                     :north),
    ("Belgium",                     :north),
    ("Bulgaria",                    :north),
    ("Canada",                      :north),
    ("China",                       :north),
    ("Denmark",                     :north),
    ("France",                      :north),
    ("Germany",                     :north),
    ("Greece",                      :north),
    ("India",                       :north),
    ("Iran",                        :north),
    ("Italy",                       :north),
    ("Japan",                       :north),
    ("Netherlands",                 :north),
    ("Norway",                      :north),
    ("Poland",                      :north),
    ("Spain",                       :north),
    ("Sweden",                      :north),
    ("United_Kingdom",              :north),
    ("United_States_of_America",    :north)]

COUNTRY_LIST_N = length(COUNTRY_LIST)

countryData = Dict( c => populateCountryDate(c, h, useOptimised = false) for (c, h) in COUNTRY_LIST)

# Update countryData[] if the COUNTRY_LIST is updated at some intermediate step
println("Populating: ")
for (c, h) in COUNTRY_LIST
    global countryData
    if haskey(countryData, c) == false
        print(c); print(" ");
        global countryData[c] = populateCountryDate(c, h, useOptimised = false)
    end
end
println()

plotCountriestoDisk("__beforeOptim")
saveParameters()

#-------------------------------------------------------------------------------------------------
#--
#-- Global optimisation across all countries
#--

# Print only 3 decimals
using Printf
Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)



# First run to optimise the counties
println("OPTIMISING COUNTRIES---------------------------")
for (c1, _) in COUNTRY_LIST
    global countryData

    # Make a note of the disease parameters
    p = countryData[c1][:params]
    println(c1)
    print("Before                     "); @show p

    countryRange = COUNTRY_RANGE
    countryRange[COUNTRY_PARAM_START] = approximateModelStartRange(c1)

    # Determine optimal parameters for each countryw
    result = bboptimize(countryData[c1][:lossFunction],
                        SearchRange = countryRange,
                        MaxTime = 30; TraceMode = :silent)

    print("After "); @show best_candidate(result)
    println();

    countryData[c1][:params] = best_candidate(result)
end

plotVignette()
saveParameters()
plotCountriestoDisk(repr(now()))

#-------------------------------------------------------------------------------------------------
# Update Epidemiology
updateEpidemiologyOnce()
plotVignette()
saveParameters()
plotCountriestoDisk(repr(now()))


#-------------------------------------------------------------------------------------------------
# Normal runs: Optimise country, then disease after each country
for run in 1:1
    println(); println("OPTIMISING 5 WORST COUNTRIES ERRORS---------------------------")
    scores = allSingleLosses(sorted = true)
    @show scores
    for c1 in first(scores, 5)[:, 2]
        global countryData

        # Make a note of the disease parameters
        p = countryData[c1][:params]
        print(c1); println("-----------------------------")
        print("Before                     "); @show p

        # Determine optimal parameters for each countryw
        countryRange = COUNTRY_RANGE
        countryRange[COUNTRY_PARAM_START] = approximateModelStartRange(c1)

        result = bboptimize(countryData[c1][:lossFunction],
                            SearchRange = countryRange,
                            MaxTime = 20; TraceMode = :silent)

        print("After "); @show best_candidate(result)

        countryData[c1][:params] = best_candidate(result)
    end

    # After having done all the countries, epidemiology again more seriously
    println("OPTIMISING EPIDEMIOLOGY---------------------------")
    # Optimise across all countries (using "France" just to have disease parameters to optimise)q

    print("Before     "); @show DiseaseParameters
    result = bboptimize(allCountriesLoss,
                        SearchRange = DISEASE_RANGE,
                        MaxTime = 20; TraceMode = :silent)

    print("After "); @show best_candidate(result)
    println()

    global DiseaseParameters = best_candidate(result)

    saveParameters()
    plotCountriestoDisk(repr(now()))

    println()
end

plotVignette()
plotCountriestoDisk(repr(now()))
saveParameters()

allDiseaseParameters = DiseaseParameters
rename!(allDiseaseParameters, DISEASE_NAMES)

@show DiseaseParameters
# DiseaseParameters = [2.092, 3.789, 6.497, 3.672, 7.346, 0.443, 1.000, 1.004, 2.828, 0.693, 8.290, 1.000]



#------------------------------------

println("OPTIMISING COUNTRIES---------------------------")
for (c1, _) in COUNTRY_LIST

    global countryData

    # Make a note of the disease parameters
    p = countryData[c1][:params]

    # Determine optimal parameters for each country
    result = bboptimize(countryData[c1][:lossFunction],
                        SearchRange = COUNTRY_RANGE,
                        MaxTime = 20; TraceMode = :compact)

    println(); println(c1)
    print("Before         "); @show p
    print("After  "); @show best_candidate(result)

    countryData[c1][:params] = best_candidate(result)
end
