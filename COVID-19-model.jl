#-------------------------------------------------------------------------------------------------
#--
#-- Ususal libraries
#--
using CSV, Dates
using DataFrames, DataFramesMeta
using Plots, PyPlot
using DifferentialEquations


# The dynamics of the epidemy is a function that mutates its argument with a precise signature
# Don't pay too much attetion to the print debugs/

function epiDynamics!(dP, P, params, t)

    S, E, I, J, H, C, R, D, K, L, BED, ICU = Pop2Comp(P)

    BED = BED[1]
    ICU = ICU[1]

    r₀, tₗ, tᵢ, tₕ, tᵤ, γₑ, γᵢ, γⱼ, γₖ, δₖ, δₗ, δᵤ , mitigation, BED_max, ICU_max = params

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

function calculateSolution(country, diseaseparams, countryparams;
                           finalDate::Union{Nothing, Date} = nothing)

    # Deconstruct the parameters
    r₀, tₗ, tᵢ, tₕ, tᵤ, γₑ, γᵢ, γⱼ, γₖ, δₖ, δₗ, δᵤ = diseaseparams
    modelStart, infectedM, infectiousM, mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9 = countryparams

    mitigation = [(0, mv0), (7, mv1), (14, mv2), (21, mv3), (35, mv4),
                  (49, mv5), (63, mv6), (77, mv7), (91, mv8), (105, mv9)]


    # First date should the date of the last death reported
    startRecordDate = first(countryData[country][:cases].time)

    # The model always starts at 0 by definition.
    startModelDay = 0

    # Final date should the date of the last death reported unless given
    # The date is shifted from actual time to model time
    if finalDate == nothing
        endRecordDate = last(countryData[country][:cases].time)
    else
        endRecordDate = finalDate
    end
    endModelDay =  timeReal2Model(endRecordDate, modelStart)

    tSpan = (startModelDay, endModelDay)

    ## Country-specific constants
    # Those are constants which cannot be changed to improve the model.

    BED_max = countryData[country][:hospital_capacity]
    ICU_max = countryData[country][:ICU_capacity]

    Age_Pyramid = Array(countryData[country][:age_distribution])
    Age_Pyramid_frac = Age_Pyramid ./ sum(Age_Pyramid)

    #TotalConfirmedAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :cases][1]
    #ConfirmedAtStart = TotalConfirmedAtStart .* Age_Pyramid_frac

#    TotalDeathsAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :deaths][1]
    TotalDeathsAtStart = DEATH_AT_MODEL_START
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



# Calculate the forecast total deaths. The calculation is performed at each time steps
# of the solutions generated by the model
function calculateTotalDeaths(sol)
    # The solutions in 'soll are presented as a vector of vectors:
    #  - it is a vector of size the number of timesteps
    #  - each element of the vector is a vector of all the variables
    nSteps = length(sol.t);
    nVars  = length(sol.u[1]);

    # Creates an Array (variables x time steps)
    solMat = reduce(hcat, sol.u)

    # Select the rows of Dx and sums to have total deaths at each time period
    forecastDeaths = sum(solMat[D0Index:D0Index + nAgeGroup - 1, :]; dims=1)

    return forecastDeaths[:]
end


# Calculate the total deaths forecast by the model on the dates for which there is an actual
# record.
# calculateTotalDeaths calculates for each time step on model time
# forecastOnActualDates calculates on the time the actual dates fall on when converted
#   to model time.
function forecastOnActualDates(sol, country::String)
    # The model always starts at time 0 to deals with the mitigation ratio (which starts at 0) and
    # assumes 5 deaths at that date.
    # The 'modelStart' parameter deals with shifting the model compared to the actual record to go
    # from 'model time' to 'real time'

    # Generate the forecast deaths on the dates of the model (in model time)
    forecastDeaths = calculateTotalDeaths(sol)

    # Get the start and final dates of the deaths record (in real time)
    startRecordDate = first(countryData[country][:cases].time)
    finalRecordDate = last(countryData[country][:cases].time)

    startRecordDay = date2days(startRecordDate)
    finalRecordDay = date2days(finalRecordDate)

    # Translates the dates into model time
    startModelDay = timeReal2Model(startRecordDate, country)
    finalModelDay = timeReal2Model(finalRecordDate, country)

    # Make a linear approximation of the forecast to match the actual days

    # -- How many steps to forecast?
    l = length(countryData[country][:cases].time)

    forecast = [ linearInterpolation(t, sol.t, forecastDeaths)
                     for t in range(startModelDay, stop = finalModelDay, length = l)]

    # Select the actual deaths record on those dates
    relevantCases = countryData[country][:cases][:, :deaths]
    relevantCases = convert(Array, relevantCases)

    return [[ t for t in startRecordDay:finalRecordDay],
            relevantCases,
            forecast]
end
