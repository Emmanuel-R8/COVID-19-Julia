# The [Neherlab COVID-19](https://neherlab.org/covid19/) forecast model
include("COVID-19-utils.jl")
include("COVID-19-model.jl")
include("COVID-19-run-model.jl")
include("COVID-19-data.jl")

# If running for the first time, or no updates for a long time
# updateData()

COUNTRY_LIST = [("France",                      :north),
                ("Germany",                     :north),
                ("Italy",                       :north),
                ("Spain",                       :north),
                ("Switzerland",                 :north),
                ("United Kingdom",              :north),
                ("United States of America",    :north)]

countryData = Dict( c => populateCountryDate(c, h) for (c, h) in COUNTRY_LIST)


#--------------------------------------------------------------------------------------------------
#--
#-- Find optimal parameters for every country
#--
for (c, _) in COUNTRY_LIST
    # Determine optimal parameters - 15 seconds per country
    result = bboptimize(deathsLoss, SearchRange = rangeParam, MaxTime = 15)
    global countryData[c][:params] = best_candidate(result)
end


# Select a country
country = "Spain"


#--------------------------------------------------------------------------------------------------
#--
#-- Default parameters
#--
sol = calculateSolution(country)
calculateTotalDeaths(sol)




#--------------------------------------------------------------------------------------------------
#-- Plot efault parameters
#--
pyplot();
clf();
ioff();

fig, ax = PyPlot.subplots();

ax.plot(countryData[country][:cases].time,
        calculateTotalDeaths(calculateSolution(country)),
        label = "Forecast");

ax.plot(countryData[country][:cases].time, countryData[country][:cases].deaths, "ro", label = "Actual", alpha = 0.3);

ax.legend(loc="lower right");
ax.set_xlabel("time");
ax.set_ylabel("Individuals");
ax.set_yscale("log");

gcf()



#--------------------------------------------------------------------------------------------------
#-- Optimal parameters
#--
result = bboptimize(deathsLoss, SearchRange = rangeParam, MaxTime = 15)
countryData[country][:params] = best_candidate(result)


#--------------------------------------------------------------------------------------------------
#-- Plot optimal parameters
#--
pyplot();
clf();
ioff();

fig, ax = PyPlot.subplots();

ax.plot(countryData[country][:cases].time,
                    calculateTotalDeaths(calculateSolution(country; params = countryData[country][:optimal])),
                    label = "Forecast");
ax.plot(countryData[country][:cases].time, countryData[country][:cases].deaths, "ro", label = "Actual", alpha = 0.3);

ax.legend(loc="lower right");
ax.set_xlabel("time");
ax.set_ylabel("Individuals");
ax.set_yscale("log");

gcf()



#-------------------------------------------------------------------------------------------------
#--
#-- Global optimisation across all countries
#--

# One run is optimising the disease, then optimising the countries.
for runs in 1:5

    # First freeze the countries and optimise the epidemiology
    freezeCountryRange()

    # Optimise across all countries
    result = bboptimize(allCountriesLoss, SearchRange = countryData["France"][:range], MaxTime = 30)

    #Update the parameters of each country
    for (c, _) in COUNTRY_LIST
        global countryData[c][:params] = (diseaseMask .* best_candidate(result)) .+ (countryMask .* countryData[c][:params])
    end

    #-- Freeze the parameters related to the disease
    freezeDiseaseRange()

    # Then optimise each country
    for (c, _) in COUNTRY_LIST
        # Make a note of the disease parameters
        p = countryData[c][:params]

        # Determine optimal parameters - 15 seconds per country
        result = bboptimize(deathsLoss, SearchRange = rangeParam, MaxTime = 15)
        global countryData[c][:params] = (diseaseMask .* p) .+ (countryMask .* best_candidate(result))
    end
end

allParams = [countryData[c][:params] for (c, _) in COUNTRY_LIST]

using Printf
Base.show(io::IO, f::Float64) = @printf io "%1.3f" f

@show allParams
