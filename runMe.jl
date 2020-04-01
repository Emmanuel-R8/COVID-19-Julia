# The [Neherlab COVID-19](https://neherlab.org/covid19/) forecast model

# Runs with 4 processes (if 4 cores CPU) - Only for optimised libraries
using Distributed
TOTAL_PROCS = 4
Distributed.addprocs(TOTAL_PROCS - Distributed.nprocs())

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
    ("Finland",                     :north),
    ("France",                      :north),
    ("Germany",                     :north),
    ("Greece",                      :north),
    ("Hungary",                     :north),
    ("Japan",                       :north),
    ("India",                       :north),
    ("Iran",                        :north),
    ("Italy",                        :north),
    ("Japan",                       :north),
    ("Netherlands",                 :north),
    ("Norway",                      :north),
    ("Poland",                      :north),
    ("Portugal",                    :north),
    ("Singapore",                   :north),
    ("Spain",                       :north),
    ("Sweden",                      :north),
    ("United_Kingdom",              :north),
    ("United_States_of_America",    :north)]

COUNTRY_LIST_N = length(COUNTRY_LIST)

countryData = Dict( c => populateCountryDate(c, h, useOptimised = false) for (c, h) in COUNTRY_LIST)

# Update countryData[] if the COUNTRY_LIST is updated at some intermediate step
print("Populating: ")
for (c, h) in COUNTRY_LIST
    global countryData
    if haskey(countryData, c) == false
        print(c); print(" ");
        global countryData[c] = populateCountryDate(c, h, useOptimised = false)
    end
end
println()

plotCountriestoDisk("beforeOptim")

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

    # Determine optimal parameters for each countryw
    result = bboptimize(countryData[c1][:lossFunction],
                        SearchRange = COUNTRY_RANGE,
                        MaxTime = 30; TraceMode = :silent)

    println(c1)
    print("Before                     "); @show p
    print("After "); @show best_candidate(result)
    println();

    countryData[c1][:params] = best_candidate(result)
end


plotCountriestoDisk("Optim1")

# Normal runs: Optimise country, then disease after each country
for run in 1:1

    println(); println("OPTIMISING COUNTRIES---------------------------")
    for (c1, _) in COUNTRY_LIST
        global countryData

        # Make a note of the disease parameters
        p = countryData[c1][:params]

        # Determine optimal parameters for each countryw
        result = bboptimize(countryData[c1][:lossFunction],
                            SearchRange = COUNTRY_RANGE,
                            MaxTime = 30; TraceMode = :silent)

        print(c1); println("-----------------------------")
        print("Before                     "); @show p
        print("After "); @show best_candidate(result)

        countryData[c1][:params] = best_candidate(result)


        # Optimise the epidemiology
        println("OPTIMISING EPIDEMIOLOGY----")
        # Optimise across all countries (using "France" just to have disease parameters to optimise)

        result = bboptimize(allCountriesLoss,
                            SearchRange = DISEASE_RANGE,
                            MaxTime = 20; TraceMode = :silent)

        print("Before     "); @show DiseaseParameters
        print("After "); @show best_candidate(result)
        println()

        global DiseaseParameters = best_candidate(result)

        println()
    end

    # After having done all the countries, epidemiology again more seriously
    println("OPTIMISING EPIDEMIOLOGY---------------------------")
    # Optimise across all countries (using "France" just to have disease parameters to optimise)

    result = bboptimize(allCountriesLoss,
                        SearchRange = DISEASE_RANGE,
                        MaxTime = 120; TraceMode = :silent)

    print("Before     "); @show DiseaseParameters
    print("After "); @show best_candidate(result)
    println()

    global DiseaseParameters = best_candidate(result)

    println()
end

plotCountriestoDisk("Optim5")

allCountryNames = DataFrame(country = [countryData[c][:name] for (c, _) in COUNTRY_LIST])
allCountryStartDates = DataFrame(startDate = [first(countryData[c][:cases].time) for (c, _) in COUNTRY_LIST])

allCountryParams = [countryData[c][:params] for (c, _) in COUNTRY_LIST]
allCountryParams = DataFrame(transpose(reduce(hcat, allCountryParams)))
rename!(allCountryParams, COUNTRY_NAMES)

allCountryParams = hcat(allCountryNames, allCountryStartDates, allCountryParams)

CSV.write("data/allCountryParameters_2020_04_01_16_51.csv", allCountryParams)

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
