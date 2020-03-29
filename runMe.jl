# The [Neherlab COVID-19](https://neherlab.org/covid19/) forecast model
include("COVID-19-utils.jl")
include("COVID-19-model.jl")
include("COVID-19-run-model.jl")
include("COVID-19-data.jl")

using BlackBoxOptim

# If running for the first time, or no updates for a long time
# updateData()

COUNTRY_LIST = [
    ("Austria",                     :north),
    ("Belgium",                     :north),
    ("Bulgaria",                    :north),
    ("Canada",                      :north),
    ("Czechia",                     :north),
    ("Denmark",                     :north),
    ("Finland",                     :north),
    ("France",                      :north),
    ("Germany",                     :north),
    ("Greece",                      :north),
    ("Hungary",                     :north),
    ("Italy",                       :north),
    ("Netherlands",                 :north),
    ("Poland",                      :north),
    ("Portugal",                    :north),
    ("Italy",                       :north),
    ("Spain",                      :north),
    ("Sweden",                      :north),
    ("United Kingdom",              :north),
    ("United States of America",    :north)]

countryData = Dict( c => populateCountryDate(c, h) for (c, h) in COUNTRY_LIST)


#--------------------------------------------------------------------------------------------------
#--
#-- Find optimal parameters for every country
#--

# Select a country
country = "Spain"

for (c, _) in COUNTRY_LIST
    # Many functions rely on country being available as a global variable.
    # A for loop creates a local scope that would prevent 'country' to be
    # visible globally
    global country = c

    # Determine optimal parameters - 15 seconds per country
    result = bboptimize(deathsLoss, SearchRange = countryData[country][:range], MaxTime = 15)
    global countryData[country][:params] = best_candidate(result)
end


using PyCall
pygui(:qt5)
using PyPlot

pyplot()
PyPlot.svg(true)

country = "Italy"

#--------------------------------------------------------------------------------------------------
#-- Plot
#--
clf();
ioff();
close(fig)

sol = calculateSolution(country,
                        countryData[country][:params];
                        finalDate = Date(2020, 9, 1))

fig, ax = PyPlot.subplots();
ax.plot(sol.t, calculateTotalDeaths(sol), label = "Forecast");
ax.plot(countryData[country][:cases].t, countryData[country][:cases].deaths, "ro", label = "Actual", alpha = 0.3);
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
# Each run is 100 + 19*10 = about 5 minutes
for runs in 1:10

    # First freeze the countries and optimise the epidemiology
    freezeCountryRange()

    # Optimise across all countries (using "France" just to have disease parameters to optimise)
    result = bboptimize(allCountriesLoss,
                        SearchRange = countryData["France"][:range],
                        MaxTime = 100)

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
        result = bboptimize(deathsLoss, SearchRange = countryData[country][:range], MaxTime = 10)
        global countryData[c][:params] = (diseaseMask .* p) .+ (countryMask .* best_candidate(result))
    end
end

allCountryNames = DataFrame(country = [countryData[c][:name] for (c, _) in COUNTRY_LIST])
allCountryStartDates = DataFrame(startDate = [first(countryData[c][:cases].time) for (c, _) in COUNTRY_LIST])

allCountryParams = [countryData[c][:params] for (c, _) in COUNTRY_LIST]
allCountryParams = DataFrame(transpose(reduce(hcat, allCountryParams)))
rename!(allCountryParams, createDefaultParameters()[3])

allCountryParams = hcat(allCountryNames, allCountryStartDates, allCountryParams)

CSV.write("data/allCountryParameters.csv", allCountryParams)
