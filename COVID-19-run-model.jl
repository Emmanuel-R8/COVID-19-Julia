#-------------------------------------------------------------------------------------------------
#--
#-- Create a range of default parameters
#--
using DataFrames, DataFramesMeta

function forecastError(country::String, actual::Array{Float64}, forecast::Array{Float64})
    # number of forcast points
    l = length(actual)

    # If logarithms are used, avoid errors.
    actual   = max.(actual,   0.0) .+ 1.0  # Should never happen
    forecast = max.(forecast, 0.0) .+ 1.0

    # Log error
    # actual = log.(actual)
    # forecast = log.(forecast)

    # Power law
    actual = sqrt.(actual)
    forecast = sqrt.(forecast)

    err_deaths = forecast .- actual

    # Penality regarding beds and icus for being negative (less forecast than actual)
    # To which extent forecast is under actual
    Δ = (max.(actual .- forecast, 0.0)).^2
    err_deaths = err_deaths .+ Δ

    # Prepare the value to never be negative (to avoid log errors)
    return sqrt( sum(err_deaths.^2) / l )
end


function forecastError(country::String, actual::Array{Int64}, forecast::Array{Float64})
    return forecastError(country, convert(Array{Float64}, actual), forecast)
end


function singleCountryLoss(country::String, diseaseparams, countryparams; finalDate = nothing)
    # finalDate = nothig to force using only the time span of actual recorded deaths
    sol = calculateSolution(country, diseaseparams, countryparams; finalDate = finalDate)

    # Extract total deaths profile
    actual = convert(Array, countryData[country][:cases][:, :deaths])
    deaths = forecastCompartmentOnActualDates(sol, "D", country)
#    beds = forecastVariableOnActualDates(sol, "BED", country)
#    icus = forecastVariableOnActualDates(sol, "ICU", country)

    return forecastError(country, actual, deaths)
    #return forecastError(country, actual, deaths, beds, icus)
end



function sumCountryLossesCountries(params)
    totalError = 0.0

    # Then each country for which the loss is immediately calculated
    for n in 1:COUNTRY_N
        country, _ = COUNTRY_LIST[n]

        country_start_index = (n - 1) * COUNTRY_N + 1
        country_final_index = (n - 1) * COUNTRY_N + COUNTRY_N
        countryparams = params[country_start_index:country_final_index]

        # loss =  sum( (log.(actual) .- log.(forecast)).^ 2 ) / length(actual)
        totalError += singleCountryLoss(country, DiseaseParameters, countryparams)
    end

    return totalError
end


#--
#-- Calculate the sum of all the losses of all the countries to optimise disease params.
#-- Loss per country is sized as if all countries had the same 1m population
#--
function sumCountryLossesDisease(diseaseparams)
    # The parameters passed to the individual loss is created with a mask defined
    totalError = 0.0

    for (country, _) in COUNTRY_LIST
        countryparams = countryData[country][:params]

        totalError += singleCountryLoss(country, diseaseparams, countryparams)
    end

    return totalError
end



function fullEpidemyLoss(params)

    # Deconstruct the entire parameter stack
    # First are the disease parameters
    diseaseparams = params[1:DISEASE_N]

    totalError = 0.0

    # Then each country for which the loss is immediately calculated
    for n in 1:COUNTRY_LIST_N
        country, _ = COUNTRY_LIST[n]

        country_start_index = DISEASE_N + (n - 1) * COUNTRY_N + 1
        country_final_index = DISEASE_N + (n - 1) * COUNTRY_N + COUNTRY_N
        countryparams = params[country_start_index:country_final_index]

        # finalDate = nothig to force using only the time span of actual recorded deaths
        sol = calculateSolution(country, diseaseparams, countryparams; finalDate = nothing)

        # Extract total deaths profile
        # Extract total deaths profile
        actual = convert(Array, countryData[country][:cases][:, :deaths])
        deaths = forecastCompartmentOnActualDates(sol, "D", country)
        beds = forecastVariableOnActualDates(sol, "BED", country)
        icus = forecastVariableOnActualDates(sol, "ICU", country)

        # loss =  sum( (log.(actual) .- log.(forecast)).^ 2 ) / length(actual)
        totalError += forecastError(country, actual, deaths, beds, icus)
    end

    return totalError
end



function updateEpidemiologyOnce(;maxtime = 60)
    # Optimise the epidemiology
    println("OPTIMISING EPIDEMIOLOGY---------------------------")

    print("Before     "); @show DiseaseParameters

    result = bboptimize(sumCountryLossesDisease,
                        SearchRange = DISEASE_RANGE;
                        Method = :adaptive_de_rand_1_bin,
                        MaxTime = maxtime,
                        TargetFitness = 2.0,
                        NThreads = Threads.nthreads(),
                        TraceMode = :compact)

    global DiseaseParameters = best_candidate(result)
    print("After      "); @show DiseaseParameters
end

function updateCountryOnce(country; maxtime = 60)
    # Make a note of the disease parameters
    println(country)
    print("Before                     ")
    @show countryData[country][:params]

    countryRange = COUNTRY_RANGE
    countryRange[COUNTRY_PARAM_START] = approximateModelStartRange(country)

    # Determine optimal parameters for each countryw
    result = bboptimize(p -> singleCountryLoss(country, DiseaseParameters, p),
                        SearchRange = countryRange;
                        Method = :adaptive_de_rand_1_bin,
                        MaxTime = maxtime,
                        TargetFitness = 2.0,
                        TraceMode = :compact)

    print("After "); @show best_candidate(result)
    println();
    global countryData[country][:params] = best_candidate(result)
end


function updateEveryCountry(; maxtime = 60)
    #-------------------------------------------------------------------------------------------------
    #--
    #-- Optimisition all countries one by one
    #--
    for (country, _) in COUNTRY_LIST
        updateCountryOnce(country; maxtime = maxtime)
    end
end
