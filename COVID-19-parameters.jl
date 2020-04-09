using Dates

# The time unit is days (as floating point)
# Day 0 is taken at 1 February 2020
const BASE_DATE = Date(2020, 2, 1)
const BASE_DAYS = 0

#%% R₀ calculations
# Default values for R_0
const BaseR₀ = 2.3954

# Seasonal forcing parameter ϵ
const ϵ = Dict(:north => 0.2, :tropical => 0.0, :south => 0.2)


#-------------------------------------------------------------------------------------------------
#--
#-- EPIDEMIOLOGY
#-

#-- Susceptibility to contagion and transition from a compartment to another
const AgeGroup = ["0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"]
const zₐ =       [0.05,   0.05,   0.10,    0.15,    0.20,    0.25,    0.30,    0.40,    0.50]
const mₐ =       [0.01,   0.03,   0.03,    0.03,    0.06,    0.10,    0.25,    0.35,    0.50]
const cₐ =       [0.05,   0.10,   0.10,    0.15,    0.20,    0.25,    0.35,    0.45,    0.55]
const fₐ =       [0.30,   0.30,   0.30,    0.30,    0.30,    0.40,    0.40,    0.50,    0.50]

const nAgeGroup = length(AgeGroup)

#-- Transition times (all in days)
# Time to infectiousness (written t\_l)
const tₗ = 5.653

# Time to infectiousness (written t\_i)
const tᵢ = 2.676

# Time in hospital bed (not ICU)
const tₕ = 2.007

# Time in ICU
const tᵤ = 13.884

# Time in asymptomatic recovery
const tᵣ = 14

#-- Contagion multiplier
# asymptomatic individuals (compartment E)
const γₑ = 0.33407

# Infected / symptomatic individuals
const γᵢ=0.25587

# Severe symptoms
const γⱼ=0.47757

# Critical symptoms
const γₖ = 5.94806

# Asymptomatic recovery
const γᵣ = 0.2

#-- Fatality Multiplier
# In ICU
const δᵤ = 1.0

# In hospital
const δₗ = 0.82117

# Out of hospital
const δₖ = 5.54388


#-------------------------------------------------------------------------------------------------
# Create parameter vectors containing decribing each parameter
#            Name             Base value      Allowed range      Disease     +/- for optim
#                                                                specific?

const DEATH_AT_MODEL_START = 1.0

const DISEASE_BASE =
         [  ["r₀",            BaseR₀,         (2.0, 3.0),        1,          1.0],
            ["tₗ",            tₗ,             (1.0, 8.0),        1,          1.0],
            ["tᵢ",            tᵢ,             (1.0, 8.0),        1,          1.0],
            ["tₕ",            tₕ,             (1.0, 8.0),        1,          1.0],
            ["tᵤ",            tᵤ,             (1.0, 25.0),       1,          2.0],
            ["tᵣ",            tᵣ,             (0.01, 40.0),      1,          2.0],
            ["γₑ",            γₑ,             (0.0, 2.0),        1,          0.3],
            ["γᵢ",            γᵢ,             (0.0, 2.0),        1,          0.0],
            ["γⱼ",            γⱼ,             (1.0, 1.0),        1,          1.0],
            ["γₖ",            γₖ,             (0.0, 10.0),       1,          2.0],
            ["γᵣ",            γᵣ,             (0.0, 10.0),       1,          2.0],
            ["δₖ",            δₖ,             (0.5, 3.0),        1,          0.5],
            ["δₗ",            δₗ,             (0.5, 9.0),        1,          1.0],
            ["δᵤ",            δᵤ,             (1.0, 1.0),        1,          0.0]]

const DISEASE_N     = size(DISEASE_BASE)[1]
const DISEASE_NAMES = [DISEASE_BASE[i][1]                        for i in 1:DISEASE_N]
const DISEASE_INIT  = [DISEASE_BASE[i][2]                        for i in 1:DISEASE_N]
const DISEASE_RANGE = Tuple{Float64, Float64}[DISEASE_BASE[i][3] for i in 1:DISEASE_N]
const DISEASE_OPTIM = [DISEASE_BASE[i][5] for i in 1:DISEASE_N]


const COUNTRY_BASE =
            # modelStart shifts between model dates which always start at 0 and real
            # dates where the model starts at a different date
            # add modelStart to go from model time to actual (eg -30 means the model starts
            # 30 days before BASE_DATE)
            # substract to go from real (eg. date of deaths) to model time
            # The range of the value is approximated by approximateModelStartRange()
         [  ["modelStart",    0.0,            (-90.0, 45.0),     0,          5.0],
            ["infectedM",     25.0,           (1.0, 5000.0),     0,          10.0],
            ["infectiousM",   25.0,           (1.0, 5000.0),     0,          10.0],

            # Mitigation profile with 4 points
            ["mv0",           1.0,            (1.0, 1.0),        0,          0.0],
            ["mv1",           0.8,            (0.1, 1.5),        0,          0.3],
            ["mv2",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv3",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv4",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv5",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv6",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv7",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv8",           0.5,            (0.1, 1.5),        0,          0.3],
            ["mv9",           0.5,            (0.1, 1.5),        0,          0.3]]

const COUNTRY_N = size(COUNTRY_BASE)[1]
const COUNTRY_NAMES = [COUNTRY_BASE[i][1]                        for i in 1:COUNTRY_N]
const COUNTRY_INIT  = [COUNTRY_BASE[i][2]                        for i in 1:COUNTRY_N]
const COUNTRY_RANGE = Tuple{Float64, Float64}[COUNTRY_BASE[i][3] for i in 1:COUNTRY_N]
const COUNTRY_OPTIM = [COUNTRY_BASE[i][5] for i in 1:COUNTRY_N]

# This is used all the time
# Get the modelStart parameter
const COUNTRY_PARAM_START = convert(Int64, findfirst("modelStart" .== COUNTRY_NAMES))
function getCountryStartDate(country::String)::Float64
   return convert(Float64, countryData[country][:params][COUNTRY_PARAM_START])
end



#-------------------------------------------------------------------------------------------------
#--
#-- Define asks indicating variables relating to the disease and to the countries
#-- Note that 2 variables are never optimised to provide an anchor for the others
#-- They are the infectiousness of a symptomatic individual γᵢ and the fatality rate from ICU δᵤ
#--
#PARAMS_diseaseMask = [PARAMS_BASE[i][4]                  for i in 1:PARAMS_N]
#PARAMS_countryMask = ones(length(PARAMS_diseaseMask)) .- PARAMS_diseaseMask


#-------------------------------------------------------------------------------------------------
#--
#-- MODELLING
#--
const COMPARTMENTS = [["S",                 "Susceptible"],
                      ["E",                 "Exposed"],
                      ["I",                 "Infectious"],
                      ["J",                 "Severe"],
                      ["H",                 "Hospitalised"],
                      ["C",                 "ICU"],
                      ["R",                 "Recovered"],
                      ["F",                 "Fully_Recovered"],
                      ["D",                 "Dead"],
                      ["K",                 "Critical"],
                      ["L",                 "Critical_Hospitalised"]]

const COMPARTMENTS_LIST = [v[1] for v in COMPARTMENTS]
const COMPARTMENTS_N    = length(COMPARTMENTS_LIST)

const VARIABLES =  [["BED",                 "Hospital_Beds"],
              ["ICU",                 "ICU_Beds"]]


#--
#-- INDIVIDUAL VARIABLES
#--
const VARIABLES_LIST = [v[1] for v in VARIABLES]
const VARIABLES_N    = length(VARIABLES_LIST)


# Where to find the index of the begiining of the compartment X
function compIndex(C:: String)
   n = findfirst(COMPARTMENTS_LIST .== C) - 1
   return n * nAgeGroup + 1
end

# Where to find the index of the variable X
function varIndex(V:: String)
   n = findfirst(VARIABLES_LIST .== V)
   return COMPARTMENTS_N * nAgeGroup + n
end
