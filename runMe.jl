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

DiseaseParameters = DISEASE_INIT

# If running for the first time, or no updates for a long time
# updateData()

COUNTRY_LIST = [
    ("France",                      :north),
    ("Italy",                       :north),
    ("Spain",                       :north),
    ("United Kingdom",              :north),
    ("United States of America",    :north)]

COUNTRY_LIST = [
    ("Austria",                     :north),
    ("Belgium",                     :north),
    ("Bulgaria",                    :north),
    ("Canada",                      :north),
    ("Czechia",                     :north),
    ("Denmark",                     :north),
    ("Finland",                     :north),
    ("France",                      :north),
    ("Germany",                     :north),
    ("Greece",                      :north),
    ("Hungary",                     :north),
    ("Italy",                       :north),
    ("Netherlands",                 :north),
    ("Poland",                      :north),
    ("Portugal",                    :north),
    ("Spain",                       :north),
    ("Sweden",                      :north),
    ("United Kingdom",              :north),
    ("United States of America",    :north)]

countryData = Dict( c => populateCountryDate(c, h, useOptimised = false) for (c, h) in COUNTRY_LIST)

using PyPlot
pyplot()

plotCountriestoDisk("beforeOptim")

#-------------------------------------------------------------------------------------------------
#--
#-- Global optimisation across all countries
#--

# Print only 3 decimals
using Printf
Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)

# One run is optimising the disease, then optimising the countries.
# Each run is 60 + 19*20 = about 5 minutes
for run in 1:5

    println("OPTIMISING COUNTRIES---------------------------")
    for (c1, _) in COUNTRY_LIST

        global countryData

        # Make a note of the disease parameters
        p = countryData[c1][:params]

        # Determine optimal parameters for each country
        result = bboptimize(countryData[c1][:lossFunction],
                            SearchRange = COUNTRY_RANGE,
                            MaxTime = 20; TraceMode = :silent)

        println(c1)
        print("Before                     "); @show p
        print("After "); @show best_candidate(result)
        println();

        countryData[c1][:params] = best_candidate(result)
    end

    # Optimise the epidemiology
    println("OPTIMISING EPIDEMIOLOGY---------------------------")
    # Optimise across all countries (using "France" just to have disease parameters to optimise)

    result = bboptimize(allCountriesLoss,
                        SearchRange = DISEASE_RANGE,
                        MaxTime = 60; TraceMode = :silent)

    print("Before     "); @show DiseaseParameters
    print("After "); @show best_candidate(result)
    println()

    global DiseaseParameters = best_candidate(result)

end

plotCountriestoDisk("Optim")






allCountryNames = DataFrame(country = [countryData[c][:name] for (c, _) in COUNTRY_LIST])
allCountryStartDates = DataFrame(startDate = [first(countryData[c][:cases].time) for (c, _) in COUNTRY_LIST])

allCountryParams = [countryData[c][:params] for (c, _) in COUNTRY_LIST]
allCountryParams = DataFrame(transpose(reduce(hcat, allCountryParams)))
rename!(allCountryParams, createDefaultParameters()[3])

allCountryParams = hcat(allCountryNames, allCountryStartDates, allCountryParams)

CSV.write("data/allCountryParameters.csv", allCountryParams)









#------------------------------------
# Optimise the epidemiology
println("OPTIMISING EPIDEMIOLOGY---------------------------")
# Optimise across all countries (using "France" just to have disease parameters to optimise)

result = bboptimize(allCountriesLoss,
                    SearchRange = DISEASE_RANGE,
                    MaxTime = 30; TraceMode = :compact)

best_candidate(result)
DiseaseParameters = best_candidate(result)


println("OPTIMISING COUNTRIES---------------------------")
for (c1, _) in COUNTRY_LIST

    global countryData

    # Make a note of the disease parameters
    p = countryData[c1][:params]

    # Determine optimal parameters for each country
    result = bboptimize(countryData[c1][:lossFunction],
                        SearchRange = COUNTRY_RANGE,
                        MaxTime = 15; TraceMode = :compact)

    println(); println(c1)
    print("Before         "); @show p
    print("After  "); @show best_candidate(result)

    countryData[c1][:params] = best_candidate(result)
end
