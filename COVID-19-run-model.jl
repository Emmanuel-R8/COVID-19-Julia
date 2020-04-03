#-------------------------------------------------------------------------------------------------
#--
#-- Create a range of default parameters
#--
using DataFrames, DataFramesMeta

function forecastError(actual, forecast)
    # Linear error
    # Prepare the value to never be negative, then add 1 (to avoid log errors)
    return sqrt(sum(( actual .- forecast).^2)) / length(actual)

    # Logarithmic error
    # Prepare the value to never be negative (to avoid log errors)
    # actual   = max.(actual,   0.00001)  # Should never happen
    # forecast = max.(forecast, 0.00001)
    #
    # return sqrt(sum( ( log.(actual) .- log.(forecast)).^2)) / length(actual)
end


function singleCountryLoss(country::String, countryparams)
    # finalDate = nothig to force using only the time span of actual recorded deaths
    sol = calculateSolution(country, DiseaseParameters, countryparams; finalDate = nothing)

    # Extract total deaths profile
    forecastDeaths = forecastOnActualDates(sol, country)

    return forecastError(forecastDeaths[2], forecastDeaths[3])
end


#--
#-- Calculate the sum of all the losses of all the countries to optimise disease params.
#-- Loss per country is sized as if all countries had the same 1m population
#--
function allCountriesLoss(diseaseparams)
    # The parameters passed to the individual loss is created with a mask defined
    totalError = 0.0

    for (country, _) in COUNTRY_LIST
        countryparams = countryData[country][:params]

        # finalDate = nothig to force using only the time span of actual recorded deaths
        sol = calculateSolution(country, diseaseparams, countryparams; finalDate = nothing)

        # Extract total deaths profile
        forecastDeaths = forecastOnActualDates(sol, country)

        # loss =  sum( (log.(actual) .- log.(forecast)).^ 2 ) / length(actual)
        totalError += forecastError(forecastDeaths[2], forecastDeaths[3])
    end

    return sqrt(totalError)
end



function fullEpidemyLoss(params)

    # Deconstruct the entire parameter stack
    # First are the disease parameters
    diseaseparams = params[1:DISEASE_N]

    totalLoss = 0.0

    # Then each country for which the loss is immediately calculated
    for n in 1:COUNTRY_LIST_N
        country, _ = COUNTRY_LIST[n]

        country_start_index = DISEASE_N + (n - 1) * COUNTRY_N + 1
        country_final_index = DISEASE_N + (n - 1) * COUNTRY_N + COUNTRY_N
        countryparams = params[country_start_index:country_final_index]

        # finalDate = nothig to force using only the time span of actual recorded deaths
        sol = calculateSolution(country, diseaseparams, countryparams; finalDate = nothing)

        # Extract total deaths profile
        forecastDeaths = forecastOnActualDates(sol, country)

        totalLoss += forecastError(forecastDeaths[2], forecastDeaths[3])
    end

    return sqrt(totalLoss)
end



function updateEpidemiologyOnce(;maxtime = 60)
    # Optimise the epidemiology
    println("OPTIMISING EPIDEMIOLOGY---------------------------")

    result = bboptimize(allCountriesLoss,
                        SearchRange = DISEASE_RANGE,
                        MaxTime = maxtime; TraceMode = :compact)

    best_candidate(result)
    global DiseaseParameters = best_candidate(result)
end

function updateCountryOnce(country; maxtime = 60)
    # Make a note of the disease parameters
    p = countryData[country][:params]

    # Determine optimal parameters for each countryw
    result = bboptimize(countryData[country][:lossFunction],
                        SearchRange = COUNTRY_RANGE,
                        MaxTime = 30; TraceMode = :compact)
    println(country)
    print("Before                     "); @show p
    print("After "); @show best_candidate(result)
    println();
    global countryData[country][:params] = best_candidate(result)
end


function updateCountriesAll(params)

    totalLoss = 0.0

    # Then each country for which the loss is immediately calculated
    for n in 1:COUNTRY_LIST_N
        country, _ = COUNTRY_LIST[n]

        country_start_index = (n - 1) * COUNTRY_N + 1
        country_final_index = (n - 1) * COUNTRY_N + COUNTRY_N
        countryparams = params[country_start_index:country_final_index]

        # finalDate = nothig to force using only the time span of actual recorded deaths
        sol = calculateSolution(country, DiseaseParameters, countryparams; finalDate = nothing)

        # Extract total deaths profile
        forecastDeaths = forecastOnActualDates(sol, country)

        totalLoss += forecastError(forecastDeaths[2], forecastDeaths[3])
    end

    return sqrt(totalLoss)
end
