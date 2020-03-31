using Dates

# The time unit is days (as floating point)
# Day 0 is taken at 1 March 2020
const BASE_DATE = Date(2020, 3, 1)
const BASE_DAYS = 0

#%% R₀ calculations
# Default values for R_0
const BaseR₀ = 2.7

# Seasonal forcing parameter \epsilon
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

#-- Contagion multiplier
# asymptomatic individuals (compartment E)
const γₑ = 0.50

# Infected / symptomatic individuals
const γᵢ=1.0

# Severe symptoms
const γⱼ=1.0

# Critical symptoms
const γₖ = 2.0

#-- Fatality Multiplier
# In ICU
const δᵤ = 1.0

# In hospital
const δₗ = 2.0

# Out of hospital
const δₖ = 3.0

#-- Transition times (all in days)
# Time to infectiousness (written t\_l)
const tₗ = 5.0

# Time to infectiousness (written t\_i)
const tᵢ = 3.0

# Time in hospital bed (not ICU)
const tₕ = 4.0

# Time in ICU
const tᵤ = 14.0


#-------------------------------------------------------------------------------------------------
# Create parameter vectors containing decribing each parameter
#            Name             Base value      Allowed range      Disease     +/- for optim
#                                                                specific?
DISEASE_BASE =
         [  ["r₀",            BaseR₀,         (1.0, 5.0),        1,          1.0],
            ["tₗ",            tₗ,             (3.0, 6.0),        1,          1.0],
            ["tᵢ",            tᵢ,             (3.0, 6.0),        1,          1.0],
            ["tₕ",            tₕ,             (4.0, 4.0),        1,          1.0],
            ["tᵤ",            tᵤ,             (12.0, 17.0),      1,          2.0],
            ["γₑ",            γₑ,             (0.2, 2.0),        1,          0.3],
            ["γᵢ",            γᵢ,             (1.0, 1.0),        1,          0.0],
            ["γⱼ",            γⱼ,             (0.2, 4.0),        1,          1.0],
            ["γₖ",            γₖ,             (0.2, 10.0),       1,          2.0],
            ["δₖ",            δₖ,             (0.5, 3.0),        1,          0.5],
            ["δₗ",            δₗ,             (0.5, 9.0),        1,          1.0],
            ["δᵤ",            δᵤ,             (1.0, 1.0),        1,          0.0]]

DISEASE_N     = size(DISEASE_BASE)[1]
DISEASE_NAMES = [DISEASE_BASE[i][1]                        for i in 1:DISEASE_N]
DISEASE_INIT  = [DISEASE_BASE[i][2]                        for i in 1:DISEASE_N]
DISEASE_RANGE = Tuple{Float64, Float64}[DISEASE_BASE[i][3] for i in 1:DISEASE_N]
DISEASE_OPTIM = [DISEASE_BASE[i][5] for i in 1:DISEASE_N]


COUNTRY_BASE =
         [  ["start",         0.0,            (-30.0, 30.0),     0,          5.0],
            ["infectedM",     50.0,           (1.0, 100.0),      0,          10.0],
            ["infectiousM",   20.0,           (1.0, 50.0),       0,          10.0],

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

COUNTRY_N = size(COUNTRY_BASE)[1]
COUNTRY_NAMES = [COUNTRY_BASE[i][1]                        for i in 1:COUNTRY_N]
COUNTRY_INIT  = [COUNTRY_BASE[i][2]                        for i in 1:COUNTRY_N]
COUNTRY_RANGE = Tuple{Float64, Float64}[COUNTRY_BASE[i][3] for i in 1:COUNTRY_N]
COUNTRY_OPTIM = [COUNTRY_BASE[i][5] for i in 1:COUNTRY_N]

#-------------------------------------------------------------------------------------------------
#--
#-- Define masks indicating variables relating to the disease and to the countries
#-- Note that 2 variables are never optimised to provide an anchor for the others
#-- They are the infectiousness of a symptomatic individual γᵢ and the fatality rate from ICU δᵤ
#--
#PARAMS_diseaseMask = [PARAMS_BASE[i][4]                  for i in 1:PARAMS_N]
#PARAMS_countryMask = ones(length(PARAMS_diseaseMask)) .- PARAMS_diseaseMask


#-------------------------------------------------------------------------------------------------
#--
#-- MODELLING
#--

const compartments =  ["S", "E", "I", "J", "H", "C", "R", "D", "K", "L"];

# How many compartments before D in the list of solution
const nb4D = findfirst(compartments .== "D") - 1

# index of D0
const D0Index = nb4D * nAgeGroup + 1
