
#-------------------------------------------------------------------------------------------------
filename = "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
download(filename, "data/ecdc.csv")

ecdc = CSV.read("data/ecdc.csv")

ecdc2 = @linq ecdc |>
  where(occursin.("France", :countriesAndTerritories)) |>
  orderby(:dateRep)

ecdc2.cases = cumsum(ecdc2.cases)
ecdc2.deaths = cumsum(ecdc2.deaths)

ecdc2[!, :time] = Date.(ecdc2.year, ecdc2.month, ecdc2.day)
ecdc2[!, :t] = date2days.(ecdc2.time)

ecdc2 = @where(ecdc2, :deaths .> 0)

#-------------------------------------------------------------------------------------------------
countryData["United_Kingdom"] = populateCountryDate("United_Kingdom", :north)
countryData["United_States_of_America"] = populateCountryDate("United_States_of_America", :north)

@show countryData["China"][:cases]


#-------------------------------------------------------------------------------------------------
@show countryData["France"]

countryData["Japan"] = populateCountryDate("Japan", :north; useOptimised = false)

@where(COUNTRY_DB[:cases], "Japan" .== :countriesAndTerritories)[1, :countryterritoryCode]


#-------------------------------------------------------------------------------------------------
country = "Iran"
sol = calculateSolution(country,
                        DiseaseParameters,
                        countryData[country][:params];
                        finalDate = nothing)

calculateTotalDeaths(sol)

# Extract total deaths profile
forecastDeaths = forecastOnActualDates(sol, country)

# Prepare the value to never be negative, then add 1 (to avoid log errors)
actual   = max.(forecastDeaths[2], 0.0) .+ 1.0
forecast = max.(forecastDeaths[3], 0.0) .+ 1.0

forecastError(forecastDeaths[2], forecastDeaths[3])




sol.t
countryData["Germany"][:params]

# loss =  sum( (log.(actual) .- log.(forecast)).^ 2 ) / length(actual)
loss =  sum( abs.(actual .- forecast)) / length(actual)

gr()
Plots.plot(actual)
Plots.plot!(forecast)

#-------------------------------------------------------------------------------------------------
#--
#-- Optimise starting date only

function f(d::Float64)
    global countryData[country][:params][1] = d
    sol = calculateSolution(country,
                            DiseaseParameters,
                            countryData[country][:params];
                            finalDate = nothing)

    # Extract total deaths profile
    forecastDeaths = forecastOnActualDates(sol, country)

    # Prepare the value to never be negative, then add 1 (to avoid log errors)
    actual   = max.(forecastDeaths[2], 0.0) .+ 1.0
    forecast = max.(forecastDeaths[3], 0.0) .+ 1.0

    return forecastError(actual, forecast)
end


result = bboptimize(f, SearchRange = [(-100.0, +30.0)], MaxTime = 15)




#-------------------------------------------------------------------------------------------------
allCountryParams = reshape(best[DISEASE_N + 1:end], (length(COUNTRY_LIST), COUNTRY_N))
allCountryParams = DataFrame(allCountryParams)
rename!(allCountryParams, COUNTRY_NAMES)

allCountryParams = hcat(DataFrame(Country = [c for (c, _) in COUNTRY_LIST]), allCountryParams)
CSV.write("data/allCountryParameters_FullOptim_2020_04_01_11_26.csv", allCountryParams)



function multiObjectiveEpidemyLoss(params)
    diseaseparams = params[1:DISEASE_N]

    # The parameters passed to the individual loss is created with a mask defined
    totalLoss = []

    n_country = length(COUNTRY_LIST)
    for n in 1:n_country
        country, _ = COUNTRY_LIST[n]

        country_start_index = DISEASE_N + (n-1) * COUNTRY_N + 1
        country_final_index = DISEASE_N + (n-1) * COUNTRY_N + COUNTRY_N
        countryparams = params[country_start_index:country_final_index]

        # finalDate = nothig to force using only the time span of actual recorded deaths
        sol = calculateSolution(country, diseaseparams, countryparams; finalDate = nothing)

        # Extract total deaths profile
        forecastDeaths = forecastOnActualDates(sol, country)

        # Prepare the value to never be negative, then add 1 (to avoid log errors)
        actual   = max.(forecastDeaths[2], 0.0) .+ 1.0
        forecast = max.(forecastDeaths[3], 0.0) .+ 1.0

        loss =  sqrt(sum( (log.(actual) .- log.(forecast)).^ 2 )) / length(actual)

        totalLoss = vcat(totalLoss, loss)
    end

    return totalLoss
end

result = bboptimize(multiObjectiveEpidemyLoss,
                    FitnessScheme = ParetoFitnessScheme{length(COUNTRY_LIST)}(is_minimizing=true),
                    SearchRange = fullRange,
                    MaxTime = 300,
                    TraceMode = :verbose)

countryData["China"][:params]

saveParameters()


#----------------- Save plot once
country = "China"
suffix = "temp"

clf();
ioff();
fig, ax = PyPlot.subplots();
sol = calculateSolution(country,
                        DiseaseParameters,
                        countryData[country][:params];
                        finalDate = Date(2020, 6, 1))

ax.plot(timeModel2Real.(sol.t, country),
        calculateTotalDeaths(sol),
        label = "Forecast");
ax.plot(countryData[country][:cases].t,
        countryData[country][:cases].deaths, "ro", label = "Actual", alpha = 0.3);
ax.legend(loc="lower right");
ax.set_title(country);
ax.set_xlabel("time");
ax.set_ylabel("Individuals");
ax.set_yscale("log");

PyPlot.savefig("images/country_" * country * "_" * suffix * ".png");


#-------------------------------------------------------------------------------------------------
#--
#-- Optimisition all countries at once
# Determine optimal parameters for each countryw
fullRange = empty([(0.0, 0.0)])
for (c1, _) in COUNTRY_LIST
    countryRange = COUNTRY_RANGE
    countryRange[COUNTRY_PARAM_START] = approximateModelStartRange(c1)
    global fullRange = vcat(fullRange, countryRange)
end

result = bboptimize(updateCountriesAll,
                    SearchRange = fullRange;
                    Method = :adaptive_de_rand_1_bin,
                    MaxTime = 60,
                    NThreads = Threads.nthreads(),
                    TraceMode = :compact)

best = best_candidate(result)
for i in 1:COUNTRY_LIST_N
    country, _ = COUNTRY_LIST[i]

    country_start_index = (i - 1) * COUNTRY_N + 1
    country_final_index = (i - 1) * COUNTRY_N + COUNTRY_N

    global countryData[country][:params] = best[country_start_index:country_final_index]
end


#-------------------------------------------------------------------------------------------------
updateEpidemiologyOnce()

result = bboptimize(allCountriesLoss,
                    SearchRange = DISEASE_RANGE;
                    Method = :adaptive_de_rand_1_bin,
                    NThreads = Threads.nthreads(),
                    MaxTime = 60,
                    TraceMode = :compact)

DiseaseParameters = best_candidate(result)


#-------------------------------------------------------------------------------------------------
country = "France"
sol = calculateSolution(country, DiseaseParameters, countryData[country][:params])

@profview sol = calculateSolution(country, DiseaseParameters, countryData[country][:params])
# Generate the forecast deaths on the dates of the model (in model time)
@profview forecastDeaths = calculateTotalDeaths(sol)

# Get the start and final dates of the deaths record (in real time)
startRecordDate = first(countryData[country][:cases].time)
finalRecordDate = last(countryData[country][:cases].time)

startRecordDay = date2days(startRecordDate)
finalRecordDay = date2days(finalRecordDate)

# Translates the dates into model time
startModelDay = timeReal2Model(startRecordDate, country)
finalModelDay = timeReal2Model(finalRecordDate, country)

# -- How many steps to forecast?
l = length(countryData[country][:cases].time)

length(range(startModelDay, stop = finalModelDay, length = l))

#-- Corresponding interval length
interval = (finalModelDay - startModelDay) / l

forecast = [ linearInterpolation(t, sol.t, forecastDeaths)
                 for t in startModelDay:interval:finalModelDay]


# Make a linear approximation of the forecast to match the actual days
forecast = [ linearInterpolation(t, sol.t, forecastDeaths) for t in startModelDay:finalModelDay]

l = length(sol.t)

# if x is before the very first point, return the first Y
# if  the list only has one point, same result
if (x <= sol.t[1]) | (l == 1)
    return yvar[1]
end

for i in 2:l
    x1 = xvar[i-1]
    x2 = xvar[i  ]

    if x < x2
        deltaY = yvar[i] - yvar[i-1]
        return yvar[i-1] + deltaY * (x - x1) / (x2 - x1)
    end
end

# If nothing has been returned, last possible choice
return yvar[l]


length(countryData["China"][:cases][:time])

#-------------------------------------------------------------------------------------------------
# Prints the number of cases for each country
[(c, length(countryData[c][:cases][:time])) for (c, _) in COUNTRY_LIST]


#-------------------------------------------------------------------------------------------------
updateData()

#-------------------------------------------------------------------------------------------------
for (c, _) in COUNTRY_LIST
    print(c * "   ")
    println(singleCountryLoss(c, countryData[c][:params]))
end


#-------------------------------------------------------------------------------------------------
#--
#-- CalculateSolution line by line
country = "France"

# Deconstruct the parameters
modelStart, infectedM, infectiousM,
    mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9 = countryData[country][:params]

mitigation = [(0, mv0), (7, mv1), (14, mv2), (21, mv3), (35, mv4),
              (49, mv5), (63, mv6), (77, mv7), (91, mv8), (105, mv9)]


# First date should the date of the last death reported
startRecordDate = first(countryData[country][:cases].time)
startModelDay = 0

endRecordDate = last(countryData[country][:cases].time)
endModelDay =  timeReal2Model(endRecordDate, modelStart)

tSpan = (startModelDay, endModelDay)

# Those are constants which cannot be changed to improve the model.
BED_max = countryData[country][:hospital_capacity]
ICU_max = countryData[country][:ICU_capacity]

Age_Pyramid = Array{Float64}(countryData[country][:age_distribution])
Age_Pyramid_frac = Age_Pyramid ./ sum(Age_Pyramid)
Population = sum(Age_Pyramid)

#TotalConfirmedAtStart = @where(countryData[country][:cases], :time .== startDate)[!, :cases][1]
#ConfirmedAtStart = TotalConfirmedAtStart .* Age_Pyramid_frac
TotalDeathsAtStart = DEATH_AT_MODEL_START
DeathsAtStart = TotalDeathsAtStart .* Age_Pyramid_frac


#-- Parameter vector
# Those are parameters which can be changed to improve the model
TotalInfected = infectedM * TotalDeathsAtStart
InfectedAtStart = TotalInfected .* Age_Pyramid_frac

TotalInfectious = infectiousM .* TotalDeathsAtStart
InfectiousAtStart = TotalInfectious .* Age_Pyramid_frac

model_params = vcat(DiseaseParameters, [mitigation, BED_max, ICU_max, Population])


#-- Compartment vector
# Note that values are initialised at 1 to avoid division by zero
S0 = Age_Pyramid .- InfectedAtStart .- InfectiousAtStart .- DeathsAtStart
E0 = InfectedAtStart
I0 = InfectiousAtStart
J0 = ones(nAgeGroup)
H0 = ones(nAgeGroup)
C0 = ones(nAgeGroup)
R0 = ones(nAgeGroup)
F0 = ones(nAgeGroup)
D0 = DeathsAtStart
K0 = ones(nAgeGroup)
L0 = ones(nAgeGroup)

# Everybody confirmed is in hospital. Assume 1 ICU bed to stay away from zeros.
BED = [TotalInfectious]
ICU = [1.0]

P0 = vcat(S0, E0, I0, J0, H0, C0, R0, F0, D0, K0, L0)

# Differential equation solver
model = ODEProblem(epiDynamics!, P0, tSpan, model_params)

# Note: progress steps might be too quick to see!
println("************* SEPARATOR ********************************************")
sol = solve(model, Tsit5(); progress = false)


#-------------------------------------------------------------------------------------------------
@show ensurePositive([ 1.0,  2.0], [1.0, 2.0])
@show ensurePositive([-1.0,  2.0], [1.0, 2.0])
@show ensurePositive([ 1.0, -2.0], [1.0, 2.0])
@show ensurePositive([-1.0, -2.0], [1.0, 2.0])

@show ensurePositive([ 1.0,  2.0], [-1.0, 2.0])
@show ensurePositive([-1.0,  2.0], [-1.0, 2.0])
@show ensurePositive([ 1.0, -2.0], [-1.0, 2.0])
@show ensurePositive([-1.0, -2.0], [-1.0, 2.0])

#-------------------------------------------------------------------------------------------------
# Update cases if the database is updated at somm intermediate step
for (c, _) in COUNTRY_LIST
    updateCountryCases(c)
end

#-------------------------------------------------------------------------------------------------
# Update countryData[] if the COUNTRY_LIST is updated at some intermediate step
println("Populating: ")
for (c, h) in COUNTRY_LIST
    global countryData
    if haskey(countryData, c) == false
        print(c); print(" ");
        global countryData[c] = populateCountryData(c, h, useOptimised = false)
    end
end
println()

#-------------------------------------------------------------------------------------------------
PRINT_DEBUG = true;
country = "France";
sol = calculateSolution(country, DiseaseParameters, countryData[country][:params])
PRINT_DEBUG = false;

#-------------------------------------------------------------------------------------------------
getVariableForecast(sol, "ICU")


#-------------------------------------------------------------------------------------------------
# finalDate = nothig to force using only the time span of actual recorded deaths
country = "France"
countryparams = countryData["France"][:params]
sol = calculateSolution(country, DiseaseParameters, countryparams; finalDate = nothing)

# Extract total deaths profile
actual = convert(Array{Float64}, countryData[country][:cases][:, :deaths])
deaths = forecastCompartmentOnActualDates(sol, "D", country)

@show actual, deaths
forecastError(country, actual, deaths)

#-------------------------------------------------------------------------------------------------
plotVignette()

country = "Belgium"
countryparams = countryData[country][:params]
sol = calculateSolution(country, DiseaseParameters, countryparams; finalDate = nothing)

sol.t
countryData["Belgium"][:cases].time
getSummedCompartment(sol, "D")

plotly()
p = Plots.plot(title = country)
xvar = timeModel2Real.(sol.t, country)
yvar = getSummedCompartment(sol, "S")
p = Plots.plot!(xvar, yvar, label = "S")


p = Plots.xaxis!("")
p = Plots.yaxis!("", :log10)

p = Plots.plot!(xvar, totalInCompartments, label = "Total")

xvar = countryData[country][:cases].t
yvar = countryData[country][:cases].deaths
p = Plots.scatter!(xvar, yvar, label = "", marker = :circle, markeralpha = 0.25)

Plots.plot(p)
