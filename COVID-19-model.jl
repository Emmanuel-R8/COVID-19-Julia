#-------------------------------------------------------------------------------------------------
#--
#-- Ususal libraries
#--
using CSV, Dates
using DataFrames, DataFramesMeta
using Plots, PyPlot
using DifferentialEquations

#-------------------------------------------------------------------------------------------------
#--
#-- Date functions
#--
# The time unit is days (as floating point)
# Day 0 is taken at 1 March 2020
const BASE_DATE = Date(2020, 3, 1)
const BASE_DAYS = 0

function date2days(d)
    return convert(Float64, datetime2rata(d) - datetime2rata(BASE_DATE))
end

function days2date(d)
    return BASE_DATE + Day(d)
end

#%% R₀ calculations
# Default values for R_0
const BaseR₀ = 2.7

# Peak date
const peakDate = Dict(
    :north => date2days(Date(2020, 1, 1)),
    :tropical => date2days(Date(2020, 1, 1)),    # although no impact
    :south => date2days(Date(2020, 7, 1))
)

# Seasonal forcing parameter \epsilon
const ϵ = Dict(:north => 0.2, :tropical => 0.0, :south => 0.2)

# Gives R_0 at a given date
function R₀(d; r₀ = baseR₀, latitude = :north)
    eps = ϵ[latitude]
    peak = peakDate[latitude]

    return r₀ * (1 + eps * cos(2.0 * π * (d - peak) / 365.25))
end

#%% Epidemy mitigation
const DEFAULT_MITIGATION = [(0, 1.00), (30, 0.80), (60, 0.5), (90, 0.20)]

function getCurrentRatio(d; start = BASE_DAYS, schedule = DEFAULT_MITIGATION)
    l = length(schedule)

    # If l = 1, ratio will be the only one
    if l == 1
        return schedule[1][2]
    else
        for i in 2:l
            d1 = schedule[i-1][1]
            d2 = schedule[i  ][1]

            if d < d2
                deltaR = schedule[i][2] - schedule[i-1][2]
                return schedule[i-1][2] + deltaR * (d - d1) / (d2 - d1)
            end
        end

        # Last possible choice
        return schedule[l][2]
    end
end


#%% Commpartments
# The population will be modeled as a single vector.
# The vector will be a stack of several vectors, each of them represents a compartment.
# Each compartment vector has a size $nAgeGroup$ representing each age group.
# The compartments are: S, E, I, H, C, R, D, K, L
# We also track the hospital bed usage BED and ICU

# Population to compartments
function Pop2Comp(P)

    # To make copy/paste less prone to error
    g = 0

    S = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    E = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    I = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    J = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    H = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    C = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    R = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    D = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    K = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1
    L = P[ g*nAgeGroup + 1: (g+1)*nAgeGroup]; g += 1

    h = 1
    BED = P[ g*nAgeGroup + h: g*nAgeGroup + h]; h += 1
    ICU = P[ g*nAgeGroup + h: g*nAgeGroup + h]; h += 1

    return S, E, I, J, H, C, R, D, K, L, BED, ICU
end


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
#--
#-- MODELLING
#--

const compartments =  ["S", "E", "I", "J", "H", "C", "R", "D", "K", "L"];

# How many compartments before D in the list of solution
const nb4D = findfirst(compartments .== "D") - 1

# index of D0
const D0Index = nb4D * nAgeGroup + 1


# Helper function to never change the number of individuals in a compartment in a way that would
# make it below 0.1 (to avoid rounding errors around 0)
function ensurePositive(d,s)
    return max.(d .+ s, 0.1) .- s
end

# The dynamics of the epidemy is a function that mutates its argument with a precise signature
# Don't pay too much attetion to the print debugs/

function epiDynamics!(dP, P, params, t)

    S, E, I, J, H, C, R, D, K, L, BED, ICU = Pop2Comp(P)

    BED = BED[1]
    ICU = ICU[1]

    r₀, tₗ, tᵢ, tₕ, tᵤ, γₑ, γᵢ, γⱼ, γₖ, δₖ, δₗ, δᵤ, mitigation, BED_max, ICU_max = params
    # println(mitigation)

    ####################################
    # Arrows reflecting epidemiology - Check signs (just in case)
    EI = ones(nAgeGroup) .* E / tₗ;  EI = max.(EI, 0.0); IE = -EI
    IJ = mₐ              .* I / tᵢ;  IJ = max.(IJ, 0.0); JI = -IJ
    JK = cₐ              .* J / tₕ;  JK = max.(JK, 0.0); KJ = -JK
    HL = cₐ              .* H / tₕ;  HL = max.(HL, 0.0); LH = -HL

    # Recovery arrows
    IR = (1 .- mₐ)       .* I / tᵢ;  IR = max.(IR, 0.0); RI = -IR
    JR = (1 .- cₐ)       .* J / tₕ;  JR = max.(JR, 0.0); RJ = -JR
    HR = (1 .- cₐ)       .* H / tₕ;  HR = max.(HR, 0.0); RH = -HR
    KR = (1 .- δₖ .* fₐ) .* K / tᵤ;  KR = max.(KR, 0.0); RK = -KR
    LR = (1 .- δₗ .* fₐ) .* L / tᵤ;  LR = max.(LR, 0.0); RL = -LR
    CR = (1 .- δᵤ .* fₐ) .* C / tᵤ;  CR = max.(CR, 0.0); RC = -CR

    # Deaths
    KD = δₖ .* fₐ        .* K / tᵤ;  KD = max.(KD, 0.0); DK = -KD
    LD = δₗ .* fₐ        .* L / tᵤ;  LD = max.(LD, 0.0); DL = -LD
    CD = δᵤ .* fₐ        .* C / tᵤ;  CD = max.(CD, 0.0); DC = -CD


    ####################################
    # Bed transfers

    ####### Step 1:
    # Decrease in ICU usage after 14 days (recall that CD and CR are vectors over the age groups)
    dICU = - (sum(CD) + sum(CR))
    dICU = ensurePositive(dICU, ICU)

    # ICU beds available
    ICU_free = ICU_max - (ICU + dICU)

    # Move as many patients as possible from $L$ to $C$ in proportion of each group
    ICU_transfer = min(sum(L), ICU_free)
    LC = ICU_transfer / sum(L) .* L;                     CL = -LC

    # Overall change in ICU bed becomes
    dICU = dICU + ICU_transfer
    dICU = ensurePositive(dICU, ICU)

    # And some normal beds are freed
    dBED = -ICU_transfer
    dBED = ensurePositive(dBED, BED)
    # println(dBED); println(BED); println(floor(BED_max)); println(ICU_transfer)

    ####### Step 2:
    # Beds available
    BED_free = BED_max - (BED + dBED)

    # Move as many patients as possible from $K$ to $L$ in proportion of each group
    BED_transfer = min(sum(K), BED_free)
    KL = BED_transfer / sum(K) .* K;                     LK = -KL

    # Overall change in normal bed becomes
    dBED = dBED + BED_transfer
    dBED = ensurePositive(dBED, BED)


    ####### Step 3:
    # Beds available
    BED_free = BED_max - (BED + dBED)

    # Move as many patients as possible from $J$ to $H$ in proportion of each group
    BED_transfer = min(sum(J), BED_free)
    JH = BED_transfer / sum(J) .* J;                     HJ = -JH

    # Overall change in ICU bed becomes
    dBED = dBED + BED_transfer
    dBED = ensurePositive(dBED, BED)


    ####################################
    # Sum of all flows + Check never negative compartment

    # Susceptible
    # Calculation of β
    β = getCurrentRatio(t; start = BASE_DAYS, schedule = mitigation) .* zₐ .* r₀

    dS = -β .* (γₑ.*E + γᵢ.*I + γⱼ.*J + γₖ.*K)
    dS = min.(-0.01, dS)
    dS = ensurePositive(dS, S)
    # println(floor(sum(dS)))

    # Exposed
    dE = -dS + IE
    dE = ensurePositive(dE, E)

    # Infected.
    dI = EI + JI + RI
    dI = ensurePositive(dI, I)

    # Infected no hospital
    dJ = IJ + HJ + KJ + RJ
    dJ = ensurePositive(dJ, J)

    # Infected in hospital
    dH = JH + LH + RH
    dH = ensurePositive(dH, H)

    # Critical no hospital
    dK = JK + LK + DK + RK
    dK = ensurePositive(dK, K)

    # Critical in hospital
    dL = KL + HL + CL + DL + RL
    dL = ensurePositive(dL, L)

    # Critical in ICU
    dC = LC + DC + RC
    dC = ensurePositive(dC, C)

    # Recovery (can only increase)
    dR = IR + JR + HR + KR + LR + CR
    dR = max.(dR, 0.01)

    # Dead (can only increase)
    dD = KD + LD + CD
    dD = max.(dD, 0.01)

    # Vector change of population and update in place
    result = vcat(dS, dE, dI, dJ, dH, dC, dR, dD, dK, dL, [dBED], [dICU])
    for i = 1:length(result)
        dP[i] = result[i]
    end
end



function errorDeaths(actualD, predictedD)
    return sqrt(sum( (actualD .- predictedD) .^ 2 ))
end


function calculateTotalDeaths(sol)
    # The solutions are returned as a vector of vectors:
    #  - it is a vector of size the number of timesteps
    #  - each element of the vector is a vector of all the variables
    nSteps = length(sol.t);
    nVars  = length(sol.u[1]);

    solMat = zeros((nSteps, nVars));
    for i = 1:nSteps
        solMat[i, :] = sol.u[i]
    end;

    # Creates matrix (varialbles x time steps)
    solMat = transpose(solMat)

    # Select the rows of Dx and sums to have total deaths at each time period
    forecastDeaths = sum(solMat[D0Index:D0Index + nAgeGroup - 1, :]; dims=1)

    # Make a linear approximation of the forecast to match the actual days
    start = date2days(first(countryData[country][:cases].time))
    finish = date2days(last(countryData[country][:cases].time))

    return [ linearInterpolation(t, sol.t, forecastDeaths) for t in start:finish]
end

function my_loss_function(sol)
   tot_loss = 0.0
   if any((s.retcode != :Success for s in sol))
     tot_loss = Inf
   else

       # Calculate the total number of deaths interpolated to actual deaths
       forecastDeaths = calculateTotalDeaths(sol)

       # Calculate RMSE
       tot_loss = sqrt(sum( (ActualDeaths .- forecastDeaths) .^ 2 ))

   end
   tot_loss
end
