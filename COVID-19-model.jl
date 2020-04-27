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
    c = 0
    S = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    E = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    I = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    J = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    H = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    C = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    R = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    F = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    D = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    K = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1
    L = P[c*nAgeGroup + 1:c*nAgeGroup + nAgeGroup]; c += 1

    r₀, tₗ, tᵢ, tₕ, tᵤ, tᵣ,
        γₑ, γᵢ, γⱼ, γₖ, γᵣ,
        δₖ, δₗ, δᵤ ,
        mitigation,
        BED_max, ICU_max, Population = params



    ####### Step 0:
    # The DE solving algorithm seems to sometimes overshoot when exploring a solution space
    # which is not acceptable with too many beds used.
    # If excess use of beds, we need to push individuals off.

    # How many ICU beds are used?
    ICU = sum(C)
    LC = CL = 0.0

    # First in ICU: move excess from C to L
    if ICU >= ICU_max
        patientMove = (ICU - ICU_max) / sum(C) .* C
        L = L .+ patientMove
        C = C .- patientMove
        CL = patientMove;                                       LC = -CL
    end

    # How many beds are used?
    BED = sum(H) + sum(L) + (ICU - ICU_max)
    JH = HJ = 0.0
    KL = LK = 0.0

    if BED >= BED_max
        # We need to move this number of patients
        patientMove =  (BED - BED_max)

        # They are moved first if severe condition (H to J);
        patientFromH = min(sum(H), patientMove); patientFromH = max(patientFromH, 0.0)

        # then moved if in critical condition (L to K)
        patientMove = patientMove - patientFromH
        patientFromL = min(sum(L), patientMove); patientFromL = max(patientFromL, 0.0)

        # Do the actual moves
        patientFromH = patientFromH / sum(H) .* H
        H = H .- patientFromH
        J = J .+ patientFromH
        HJ = patientFromH;                                      JH = -HJ

        patientFromL = patientFromL / sum(L) .* L
        L = L .- patientFromL
        K = K .+ patientFromL
        LK = patientFromL;                                      KL = -LK
    end



    ####################################
    # Arrows reflecting epidemiology - Check signs (just in case)
    # EI means from  E to I
    # EI > 0  means flow from E to I (E goes down, I goes up)
    EI = ones(nAgeGroup) .* E / tₗ;  EI = max.(EI, 0.0001);     IE = -EI
    IJ = mₐ              .* I / tᵢ;  IJ = max.(IJ, 0.0001);     JI = -IJ
    JK = cₐ              .* J / tₕ;  JK = max.(JK, 0.0001);     KJ = -JK
    HL = cₐ              .* H / tₕ;  HL = max.(HL, 0.0001);     LH = -HL

    # Asymptomatic recovery arrows
    IR = (1 .- mₐ)       .* I / tᵢ;  IR = max.(IR, 0.0001);     RI = -IR
    JR = (1 .- cₐ)       .* J / tₕ;  JR = max.(JR, 0.0001);     RJ = -JR
    HR = (1 .- cₐ)       .* H / tₕ;  HR = max.(HR, 0.0001);     RH = -HR
    KR = (1 .- (δₖ * δᵤ) .* fₐ) .* K / tᵤ;  KR = max.(KR, 0.0001);     RK = -KR
    LR = (1 .- (δₗ * δᵤ) .* fₐ) .* L / tᵤ;  LR = max.(LR, 0.0001);     RL = -LR
    CR = (1 .-       δᵤ  .* fₐ) .* C / tᵤ;  CR = max.(CR, 0.0001);     RC = -CR

    # Full recovery arrows
    RF = ones(nAgeGroup) .* R / tᵣ;  RF = max.(RF, 0.0001);     FR = -RF

    # Deaths
    KD =       δₖ .* fₐ  .* K / tᵤ;  KD = max.(KD, 0.0001);     DK = -KD
    LD =       δₗ .* fₐ  .* L / tᵤ;  LD = max.(LD, 0.0001);     DL = -LD
    CD =       δᵤ .* fₐ  .* C / tᵤ;  CD = max.(CD, 0.0001);     DC = -CD


    ####################################
    # Bed transfers

    ICU = sum(C)
    BED = sum(H) + sum(L)

    ####### Step 1:
    # Decrease in ICU usage after 14 days (recall that CD and CR are vectors over the age groups)
    dICU = -(sum(CD) + sum(CR))
    dICU = ensurePositive(dICU, ICU)
    if PRINT_DEBUG @show dICU end

    # ICU beds available
    ICU_free = max(0.0, ICU_max - (ICU + dICU))
    if PRINT_DEBUG
        @show ICU_max
        @show ICU
        @show dICU
        @show ICU_free
    end

    # Move as many patients as possible from $L$ to $C$ in proportion of each group
    ICU_transfer = min(sum(L), ICU_free)
    LC = LC .+ ICU_transfer / sum(L) .* L;                      CL = -LC
    if PRINT_DEBUG
        @show ICU_transfer
        @show sum(LC)
    end

    # Overall change in ICU bed becomes
    dICU = dICU + ICU_transfer
    dICU = ensurePositive(dICU, ICU)

    # And some normal beds are freed
    dBED = -ICU_transfer
    dBED = ensurePositive(dBED, BED)
    # println(dBED); println(BED); println(floor(BED_max)); println(ICU_transfer)

    ####### Step 2:
    # Beds available
    BED_free = max(0.0, BED_max - (BED + dBED))

    # Move as many patients as possible from $K$ to $L$ in proportion of each group
    BED_transfer = min(sum(K), BED_free)
    KL = KL .+ BED_transfer / sum(K) .* K;                      LK = -KL

    # Overall change in normal bed becomes
    dBED = dBED + BED_transfer
    dBED = ensurePositive(dBED, BED)


    ####### Step 3:
    # Beds available
    BED_free = max(0.0, BED_max - (BED + dBED))

    # Move as many patients as possible from $J$ to $H$ in proportion of each group
    BED_transfer = min(sum(J), BED_free)
    JH = JH .+ BED_transfer / sum(J) .* J;                      HJ = -JH

    # Overall change in ICU bed becomes
    dBED = dBED + BED_transfer
    dBED = ensurePositive(dBED, BED)


    ####################################
    # Sum of all flows + Check never negative compartment

    # Susceptible
    # Calculation of β
    mitigationRatio = getCurrentRatio(t; start = BASE_DAYS, schedule = mitigation)
    β = (r₀ * mitigationRatio) .* zₐ

    dS = - sum(γₑ.*E + γᵢ.*I + γⱼ.*J + γₖ.*K + γᵣ.*R)  / Population .* (S .* β)
    dS = min.(-0.0001, dS)
    dS = ensurePositive(dS, S)

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

    # Asymptomatic recovery (can only increase)
    dR = IR + JR + HR + KR + LR + CR + FR
    dR = ensurePositive(dR, R)

    # Full recovery (can only increase)
    dF = RF
    dF = max.(dF, 0.0001)

    # Dead (can only increase)
    dD = KD + LD + CD
    dD = max.(dD, 0.0001)

    # Vector change of population and update in place
    result = vcat(dS, dE, dI, dJ, dH, dC, dR, dF, dD, dK, dL)

    for i = 1:length(result)
        dP[i] = result[i]
    end

    if PRINT_DEBUG
        pf = map(sum, [dS, dE, dI, dJ, dH, dK, dL, dC, dD, dR, dF,
                       [dBED], [BED], [BED_max], [dICU], [ICU], [ICU_max],
                       dS + dE + dI + dJ + dH + dK + dL + dC + dD + dR + dF])
        @show pf
        #readline()
    end
end


function calculateSolution(country, diseaseparams, countryparams;
                           finalDate::Union{Nothing, Date} = nothing)

    # Deconstruct the parameters
    r₀, tₗ, tᵢ, tₕ, tᵤ, tᵣ,
        γₑ, γᵢ, γⱼ, γₖ, γᵣ,
        δₖ, δₗ, δᵤ = diseaseparams

    modelStart, infectedM, infectiousM, mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9 = countryparams

    mitigation = [(0, mv0), (7, mv1), (14, mv2), (21, mv3), (35, mv4),
                  (42, mv5), (49, mv6), (63, mv7), (77, mv8), (91, mv9)]


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

    Age_Pyramid = Array{Float64}(countryData[country][:age_distribution])
    Age_Pyramid_frac = Age_Pyramid ./ sum(Age_Pyramid)
    Population = sum(Age_Pyramid)

    #TotalConfirmedAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :cases][1]
    #ConfirmedAtStart = TotalConfirmedAtStart .* Age_Pyramid_frac

    # TotalDeathsAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :deaths][1]
    TotalDeathsAtStart = DEATH_AT_MODEL_START
    DeathsAtStart = TotalDeathsAtStart .* Age_Pyramid_frac


    ## Parameter vector
    # Those are parameters which can be changed to improve the model
    TotalInfected = infectedM * TotalDeathsAtStart
    InfectedAtStart = TotalInfected .* Age_Pyramid_frac

    TotalInfectious = infectiousM .* TotalDeathsAtStart
    InfectiousAtStart = TotalInfectious .* Age_Pyramid_frac

    model_params = [r₀,
                   tₗ, tᵢ, tₕ, tᵤ, tᵣ,
                   γₑ, γᵢ, γⱼ, γₖ, γᵣ,
                   δₖ, δₗ, δᵤ,
                   mitigation,
                   BED_max,
                   ICU_max,
                   Population]


    ## Compartment vector
    # Note that values are initialised at 1 to avoid division by zero
    S0 = Age_Pyramid .- InfectedAtStart .- InfectiousAtStart .- DeathsAtStart
    E0 = InfectedAtStart
    I0 = InfectiousAtStart
    J0 = 0.0001 .* ones(Float64, nAgeGroup)
    H0 = 0.0001 .* ones(Float64, nAgeGroup)
    C0 = 0.0001 .* ones(Float64, nAgeGroup)
    R0 = 0.0001 .* ones(Float64, nAgeGroup)
    F0 = 0.0001 .* ones(Float64, nAgeGroup)
    D0 = DeathsAtStart
    K0 = 0.0001 .* ones(Float64, nAgeGroup)
    L0 = 0.0001 .* ones(Float64, nAgeGroup)

    P0 = vcat(S0, E0, I0, J0, H0, C0, R0, F0, D0, K0, L0)

    # Differential equation solver
    model = ODEProblem(epiDynamics!, P0, tSpan, model_params)

    # Note: progress steps might be too quick to see!
    sol = solve(model, Tsit5(); progress = false)

    return sol
end


# Calculate the forecast total in a given compartment summed across al ages.
# The calculation is performed at each time steps
# of the solutions generated by the model
function getSummedCompartment(sol, C::String)
    # The solutions in 'sol' are presented as a vector of vectors:
    #  - it is a vector of size the number of timesteps
    #  - each element of the vector is a vector of all the variables
    # Creates an Array (variables x time steps)
    solMat = reduce(hcat, sol.u)
    l = size(solMat)[2]

    # Select the rows of Dx and sums to have total deaths at each time period
    result = sum(solMat[compIndex(C):compIndex(C) + nAgeGroup - 1, :]; dims=1)

    # Reshape to have a 1-D vector
    return reshape(result, l)
end


# Calculate the forecast total deaths. The calculation is performed at each time steps
# of the solutions generated by the model
function getVariableForecast(sol, V::String)
    # Creates an Array (variables x time steps)
    solMat = reduce(hcat, sol.u)

    # Select the rows of Dx and sums to have total deaths at each time period
    return solMat[varIndex(V), :]
end


# Calculate the total deaths forecast by the model on the dates for which there is an actual
# record.
# calculateTotalDeaths calculates for each time step on model time
# forecastOnCompartmentActualDates calculates on the time the actual dates fall on when converted
#   to model time.
function forecastCompartmentOnActualDates(sol, C::String, country::String)
    # The model always starts at time 0 to deals with the mitigation ratio (which starts at 0)
    # The 'modelStart' parameter deals with shifting the model compared to the actual record to go
    # from 'model time' to 'real time'

    # Generate the forecast deaths on the dates of the model (in model time)
    c = getSummedCompartment(sol, C)

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

    return [ linearInterpolation(t, sol.t, c) for t in range(startModelDay, stop = finalModelDay, length = l)]
end


function forecastVariableOnActualDates(sol, C::String, country::String)
    # The model always starts at time 0 to deals with the mitigation ratio (which starts at 0)
    # The 'modelStart' parameter deals with shifting the model compared to the actual record to go
    # from 'model time' to 'real time'

    # Generate the forecast deaths on the dates of the model (in model time)
    v = getVariableForecast(sol, C)

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

    return [ linearInterpolation(t, sol.t, v) for t in range(startModelDay, stop = finalModelDay, length = l)]
end
