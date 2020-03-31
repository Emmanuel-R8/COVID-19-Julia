#-------------------------------------------------------------------------------------------------
#--
#-- Create a range of default parameters
#--
using DataFrames, DataFramesMeta

#-------------------------------------------------------------------------------------------------
#--
#-- Freeze the range of parameters specific to a given country to flex the disease-specific ones
#--
function freezeCountryRange()
    # List of parameters to freeze
    listParams = ["γᵢ", "δᵤ", "mv0",                           # Always frozen at 1.0
                  "start", "infectedM", "infectiousM",
                  "mv1", "mv2", "mv3", "mv4",
                  "mv5", "mv6", "mv7", "mv8", "mv9"]

    # for each parameter
    for position in 1:PARAMS_N
        for (c, _) in COUNTRY_LIST
            value = countryData[c][:params][position]

            # If the name is to be relaxed (not in the list to be frozen)
            if findfirst(PARAMS_NAMES[position] .== listParams) == nothing
                # Sets the range at the value +/- optimisation range
                countryData[c][:range][position] = (value - PARAMS_OPTIM[position],
                                                    value + PARAMS_OPTIM[position])
            else
                # Otherwise, sets the optimisation range at that value
                countryData[c][:range][position] = (value, value)
            end
        end
    end
end

function freezeDiseaseRange()
#-- This function is the same code as freezeCountryRange, but different list

    # List of parameters to freeze
    listParams = ["γᵢ", "δᵤ", "mv0",                             # Always frozen at 1.0
                  "r₀", "tₗ", "tᵢ", "tₕ", "tᵤ",
                  "γₑ", "γⱼ", "γₖ",
                  "δₖ", "δₗ"]

    # for each parameter
    for position in 1:PARAMS_N
        for (c, _) in COUNTRY_LIST
            value = countryData[c][:params][position]

            # If the name is to be relaxed (not in the list to be frozen)
            if findfirst(PARAMS_NAMES[position] .== listParams) == nothing
                # Sets the range at the value +/- 20%
                countryData[c][:range][position] = (value - PARAMS_OPTIM[position],
                                                    value + PARAMS_OPTIM[position])
            else
                # Otherwise, sets the optimisation range at that value
                countryData[c][:range][position] = (value, value)
            end
        end
    end
end



function singleCountryLoss(country, countryparams)
    # finalDate = nothig to force using only the time span of actual recorded deaths
    sol = calculateSolution(country, DiseaseParameters, countryparams; finalDate = nothing)

    # Extract total deaths profile
    forecastDeaths = forecastOnActualDates(sol, country)

    return sqrt(sum((forecastDeaths[2] .- forecastDeaths[3]) .^ 2 ))
end


#--
#-- Calculate the sum of all the losses of all the countries to optimise disease params.
#-- Loss per country is sized as if all countries had the same 1m population
#--
function allCountriesLoss(diseaseparams)

    # The parameters passed to the individual loss is created with a mask defined
    totalLoss = 0.0

    for (country, _) in COUNTRY_LIST
        countryparams = countryData[country][:params]

        # finalDate = nothig to force using only the time span of actual recorded deaths
        sol = calculateSolution(country, diseaseparams, countryparams; finalDate = nothing)

        # Extract total deaths profile
        forecastDeaths = forecastOnActualDates(sol, country)
        loss =  sqrt(sum((forecastDeaths[2] .- forecastDeaths[3]) .^ 2 ))

        totalLoss += loss * 1000000.0 / countryData[country][:population]
    end

    return totalLoss
end



function calculateSolution(country, diseaseparams, countryparams; finalDate = nothing)

    # Deconstruct the parameters
    r₀, tₗ, tᵢ, tₕ, tᵤ, γₑ, γᵢ, γⱼ, γₖ, δₖ, δₗ, δᵤ = diseaseparams
    shiftStart, infectedM, infectiousM, mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9 = countryparams

    mitigation = [(0, mv0), (7, mv1), (14, mv2), (21, mv3), (28, mv4),
                  (35, mv5), (42, mv6), (49, mv7), (56, mv8), (63, mv9)]


    # First date should the date of the last death reported
    startDate = first(countryData[country][:cases].time)

    # Then the start of the model is shifted around the first death
    startDays = first(countryData[country][:cases].t) + shiftStart

    # Final date should the date of the last death reported unless given
    if finalDate == nothing
        endDate = last(countryData[country][:cases].time)
        endDays = last(countryData[country][:cases].t)
    else
        endDate = finalDate
        endDays = date2days(finalDate)
    end

    tSpan = (startDays, endDays)

    ## Country-specific constants
    # Those are constants which cannot be changed to improve the model.

    BED_max = countryData[country][:hospital_capacity]
    ICU_max = countryData[country][:ICU_capacity]

    Age_Pyramid = transpose( Matrix(countryData[country][:age_distribution]))
    Age_Pyramid_frac = Age_Pyramid / sum(Age_Pyramid)

    #TotalConfirmedAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :cases][1]
    #ConfirmedAtStart = TotalConfirmedAtStart .* Age_Pyramid_frac

#    TotalDeathsAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :deaths][1]
    TotalDeathsAtStart = 5
    DeathsAtStart = TotalDeathsAtStart .* Age_Pyramid_frac


    ## Parameter vector
    # Those are parameters which can be changed to improve the model
    TotalInfected = infectedM * TotalDeathsAtStart
    InfectedAtStart = TotalInfected .* Age_Pyramid_frac

    TotalInfectious = infectiousM .* TotalDeathsAtStart
    InfectiousAtStart = TotalInfectious .* Age_Pyramid_frac

    model_params = [r₀,
                   tₗ, tᵢ, tₕ, tᵤ,
                   γₑ, γᵢ, γⱼ, γₖ,
                   δₖ, δₗ, δᵤ,
                   mitigation,
                   BED_max,
                   ICU_max]


    ## Compartment vector
    # Note that values are initialised at 1 to avoid division by zero
    S0 = Age_Pyramid .- InfectedAtStart .- InfectiousAtStart .- DeathsAtStart
    E0 = InfectedAtStart
    I0 = InfectiousAtStart
    J0 = ones(nAgeGroup)
    H0 = ones(nAgeGroup)
    C0 = ones(nAgeGroup)
    R0 = ones(nAgeGroup)
    D0 = DeathsAtStart
    K0 = ones(nAgeGroup)
    L0 = ones(nAgeGroup)

    # Everybody confirmed is in hospital. Assume 1 ICU bed to stay away from zeros.
    BED = [TotalInfectious]
    ICU = [1.0]

    P0 = vcat(S0, E0, I0, J0, H0, C0, R0, D0, K0, L0, BED, ICU)

    # Differential equation solver
    model = ODEProblem(epiDynamics!, P0, tSpan, model_params)

    # Note: progress steps might be too quick to see!
    sol = solve(model, Tsit5(); progress = false)

    return sol
end
