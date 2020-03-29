#-------------------------------------------------------------------------------------------------
#--
#-- Create a range of default parameters
#--
using DataFrames, DataFramesMeta


function createDefaultParameters()
    # Create a vector containing triplets of (parameter name, initial value, range)
    params = [  ["r₀",            BaseR₀,         (2.0, 5.0)],
                ["tₗ",            4,              (3.0, 6.0)],
                ["tᵢ",            4,              (3.0, 6.0)],
                ["tₕ",            4,              (4.0, 4.0)],
                ["tᵤ",            14,             (12.0, 17.0)],
                ["γₑ",            0.8,            (0.2, 1.0)],
                ["γᵢ",            1.0,            (1.0, 1.0)],
                ["γⱼ",            1.5,            (0.2, 2.0)],
                ["γₖ",            2.0,            (0.2, 3.0)],
                ["δₖ",            1.0,            (0.9, 1.1)],
                ["δₗ",            2.0,            (0.8, 3.0)],
                ["δᵤ",            1.0,            (1.0, 1.0)],
                ["infectedM",     10.0,           (1.0, 50.0)],
                ["infectiousM",   20.0,           (1.0, 50.0)],

                # Mitigation profile with 4 points
                ["mv0",           1.0,            (1.0, 1.0)],
                ["mv1",           0.8,            (0.1, 1.5)],
                ["mv2",           0.5,            (0.1, 1.5)],
                ["mv3",           0.5,            (0.1, 1.5)],
                ["mv4",           0.5,            (0.1, 1.5)],
                ["mv5",           0.5,            (0.1, 1.5)],
                ["mv6",           0.5,            (0.1, 1.5)],
                ["mv7",           0.5,            (0.1, 1.5)],
                ["mv8",           0.5,            (0.1, 1.5)],
                ["mv9",           0.5,            (0.1, 1.5)]]

                # Fixed parameters
#                ["BED_max",       1000.0,         [1000.0, 1000.0]],
#                ["ICU_max",       1000.0,         [1000.0, 1000.0]]]

    nParams = size(params)[1]
    pNames = [params[i][1] for i in 1:nParams]
    initialParam = [params[i][2] for i in 1:nParams]
    rangeParam = Tuple{Float64, Float64}[params[i][3] for i in 1:nParams]

    return(initialParam, rangeParam, pNames)

end

#-------------------------------------------------------------------------------------------------
#--
#-- Define masks indicating variables relating to the disease and to the countries
#-- Note that 2 variables are never optimised to provide an anchor for the others
#-- They are the infectiousness of a symptomatic individual γᵢ and the fatality rate from ICU δᵤ
#--
diseaseMask = [1,
               1, 1, 1, 1,     # t parameters
               1, 0, 1, 1,     # γ parameters
               1, 1, 0,        # δ parameters
               0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

countryMask = ones(length(diseaseMask)) .- diseaseMask


#-------------------------------------------------------------------------------------------------
#--
#-- Freeze the range of parameters specific to a given country to flex the disease-specific ones
#--
function freezeCountryRange()
    # List of parameters to freeze
    listParams = ["γᵢ", "δᵤ", "mv0",                           # Always frozen
                  "infectedM", "infectiousM",
                  "mv1", "mv2", "mv3", "mv4",
                  "mv5", "mv6", "mv7", "mv8", "mv9"]

    # Get the names of the parameters
    _, _, pNames = createDefaultParameters()

    for (c, _) in COUNTRY_LIST

        # for each parameter
        for position in 1:length(pNames)

            value = countryData[c][:params][position]

            # If the name is to be relaxed (not in the list)
            if findfirst(pNames[position] .== listParams) == nothing
                # Sets the range at the value +/- 20%
                countryData[c][:range][position] = (value * (1-0.20), value * (1+0.20))
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
    listParams = ["γᵢ", "δᵤ", "mv0",                             # Always frozen
                  "r₀", "tₗ", "tᵢ", "tₕ", "tᵤ",
                  "γₑ", "γⱼ", "γₖ",
                  "δₖ", "δₗ"]

    # Get the names of the parameters
    _, _, pNames = createDefaultParameters()

    for (c, _) in COUNTRY_LIST

        # for each parameter
        for position in 1:length(pNames)

            value = countryData[c][:params][position]

            # If the name is to be relaxed (not in the list)
            if findfirst(pNames[position] .== listParams) == nothing
                # Sets the range at the value +/- 20%
                countryData[c][:range][position] = (value * (1-0.20), value * (1+0.20))
            else
                # Otherwise, sets the optimisation range at that value
                countryData[c][:range][position] = (value, value)
            end
        end
    end
end



function deathsLoss(p)
    # finalDate = nothig to force using only the time span of actual recorded deaths
    sol = calculateSolution(country, p; finalDate = nothing)

    # Extract total deaths profile
    forecastDeaths = forecastOnActualDates(sol)
    #println(sum(forecastDeaths)); println();

    return sqrt(sum( (countryData[country][:cases].deaths .- forecastDeaths) .^ 2 ))
end


#--
#-- Calculate the sum of all the losses of all the countries to optimise disease params.
#--
function allCountriesLoss(p)

    # The parameters passed to the individual loss is created with a mask defined
    totalLoss = 0.0

    for (country, _) in COUNTRY_LIST

        countryParams = countryData[country][:params]
        totalLoss += deathsLoss((diseaseMask .* p) .+ (countryMask .* countryParams))
    end

    return totalLoss
end

function calculateSolution(country, params; finalDate = nothing)

    # First date should the date of the last death reported
    startDate = first(countryData[country][:cases].time)
    startDays = first(countryData[country][:cases].t)

    # Final date should the date of the last death reported
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

    TotalConfirmedAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :cases][1]
    ConfirmedAtStart = TotalConfirmedAtStart .* Age_Pyramid_frac

    TotalDeathsAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :deaths][1]
    DeathsAtStart = TotalDeathsAtStart .* Age_Pyramid_frac


    # Deconstruct the parameters
    r₀, tₗ, tᵢ, tₕ, tᵤ, γₑ, γᵢ, γⱼ, γₖ, δₖ, δₗ, δᵤ, infectedM, infectiousM,
        mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9 = params

    mitigation = [(0, mv0), (7, mv1), (14, mv2), (21, mv3), (28, mv4),
                  (35, mv5), (42, mv6), (49, mv7), (56, mv8), (63, mv9)]

    ## Parameter vector
    # Those are parameters which can be changed to improve the model
    TotalInfected = infectedM * TotalConfirmedAtStart
    InfectedAtStart = TotalInfected .* Age_Pyramid_frac

    TotalInfectious = infectiousM .* TotalConfirmedAtStart
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
