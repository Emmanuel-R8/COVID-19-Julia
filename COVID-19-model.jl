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


function errorDeaths(actualD, predictedD)
    return sqrt(sum( (actualD .- predictedD) .^ 2 ))
end


function calculateTotalDeaths(sol)
    # The solutions are returned as a vector of vectors:
    #  - it is a vector of size the number of timesteps
    #  - each element of the vector is a vector of all the variables
    nSteps = length(sol.t);
    nVars  = length(sol.u[1]);

    # solMat = zeros((nSteps, nVars));
    # for i = 1:nSteps
    #     solMat[i, :] = sol.u[i]
    # end;
    #
    # # Creates matrix (varialbles x time steps)
    # solMat = transpose(solMat)
    solMat = reduce(hcat, sol.u)

    # Select the rows of Dx and sums to have total deaths at each time period
    forecastDeaths = sum(solMat[D0Index:D0Index + nAgeGroup - 1, :]; dims=1)

    return forecastDeaths[:]
end

function forecastOnActualDates(sol, country)

    # Calculate the forecast total deaths
    forecastDeaths = calculateTotalDeaths(sol)

    # What is the starting date of the model
    position = (findfirst("start" .== COUNTRY_NAMES))
    startDate = ceil(countryData[country][:params][position])

    # Make a linear approximation of the forecast to match the actual days
    start = max(date2days(first(countryData[country][:cases].time)), startDate)
    finish = date2days(last(countryData[country][:cases].time))

    relevantCases =  @linq countryData[country][:cases] |> where(:t .>= start) |> select(:deaths)
    relevantCases = convert(Array, relevantCases)

    return [[ t                                             for t in start:finish],
            relevantCases,
            [ linearInterpolation(t, sol.t, forecastDeaths) for t in start:finish]]
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
