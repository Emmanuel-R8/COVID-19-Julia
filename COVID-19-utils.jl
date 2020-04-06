using Dates

#-------------------------------------------------------------------------------------------------
#--
#-- Linear interpolation of the Y value at point x depending on where it falls in/around X.
#--
function linearInterpolation(x, xvar, yvar)
    l = length(xvar)

    # if x is before the very first point, return the first Y
    # if  the list only has one point, same result
    if (x <= xvar[1]) | (l == 1)
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
end

#-------------------------------------------------------------------------------------------------
#--
#-- Date functions
#--
function date2days(d::Date)::Float64
    return convert(Float64, datetime2rata(d) - datetime2rata(BASE_DATE))
end

function days2date(d::Float64)::Date
    return BASE_DATE + Day(d)
end

# Converts a time (in days from BASE_DATE) to a time point in model time
# Example: if modelStart is 30, the model starts (model time = 0) 30 days after BASE_DATE
# Real date 45 days after BASE_DATE is time = 15 in model time.
function timeReal2Model(real_time::Float64, modelStart::Float64)
    return real_time - modelStart
end

function timeModel2Real(model_time::Float64, modelStart::Float64)
    return model_time + modelStart
end

# Equivalent functions for Dates
function timeReal2Model(real_date::Date, modelStart::Float64)
    return timeReal2Model(date2days(real_date), modelStart)
end

function timeModel2Real(model_date::Date, modelStart::Float64)
    return timeModel2Real(date2days(model_date), modelStart)
end

# Equivalent functions for countries as Strings
function timeReal2Model(real_time, country::String)
    return timeReal2Model(real_time, getCountryStartDate(country))
end

function timeModel2Real(model_time, country::String)
    return timeModel2Real(model_time, getCountryStartDate(country))
end



# Peak date
const peakDate = Dict(
    :north => date2days(Date(2020, 1, 1)),
    :tropical => date2days(Date(2020, 1, 1)),    # although no impact
    :south => date2days(Date(2020, 7, 1))
)

# Gives R_0 at a given date
function R₀(d; r₀ = baseR₀, latitude = :north)
    eps = ϵ[latitude]
    peak = peakDate[latitude]

    return r₀ * (1 + eps * cos(2.0 * π * (d - peak) / 365.25))
end

#%% Epidemy mitigation
const DEFAULT_MITIGATION = [(0, 1.0),  (7, 0.8), (14, 0.5), (21, 0.5), (35, 0.5),
                            (49, 0.5), (63, 0.5), (77, 0.5), (91, 0.5), (105, 0.5)]

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


# Helper function to never change the number of individuals in a compartment in a way that would
# make it below 0.1 (to avoid rounding errors around 0)
function ensurePositive(d,s)
    return max.(d .+ s, 0.1) .- s
end


function plotCountriestoDisk(suffix)
    for (c, _) in COUNTRY_LIST
        global country = c

        clf();
        ioff();
        fig, ax = PyPlot.subplots();
        sol = calculateSolution(country,
                                DiseaseParameters,
                                countryData[country][:params];
                                finalDate = Date(2020, 8, 1))

        ax.plot(timeModel2Real.(sol.t, country),
                calculateTotalDeaths(sol),
                label = "Forecast");

        ax.plot(countryData[country][:cases][:t],
                countryData[country][:cases].deaths, "ro", label = "Actual", alpha = 0.3);

        ax.legend(loc="lower right");
        ax.set_title(country);
        ax.set_xlabel("time");
        ax.set_ylabel("Individuals");
        ax.set_yscale("log");

        PyPlot.savefig("images/country_" * country * "_" * suffix * ".png");
    end
end


# Generates the current training loss of all countries.
function allSingleLosses(; sorted = false)
    losses = [(c, singleCountryLoss(c, countryData[c][:params])) for (c, _) in COUNTRY_LIST]
    losses = [ [e for (_, e) in losses], [c for (c, _) in losses] ]
    sort(DataFrame(losses), rev = sorted)
end



function saveParameters()
    allCountryNames = DataFrame(country = [countryData[c][:name] for (c, _) in COUNTRY_LIST])

    allCountryParams = [countryData[c][:params] for (c, _) in COUNTRY_LIST]
    allCountryParams = DataFrame(transpose(reduce(hcat, allCountryParams)))
    rename!(allCountryParams, COUNTRY_NAMES)

    allCountryParams = hcat(allCountryNames, allCountryParams)

    DF = DataFrame(DiseaseParameters')
    rename!(DF, DISEASE_NAMES)

    nowString = repr(now())
    CSV.write("data/" * nowString * "_CountryParameters.csv", allCountryParams)
    CSV.write("data/" * nowString * "_DiseaseParameters.csv", DF)
end


function plotVignette()
    plotly()

    plot_dict  = Dict()
    for (country, _) in COUNTRY_LIST
        sol = calculateSolution(country,
                                DiseaseParameters,
                                countryData[country][:params];
                                finalDate = Date(2020, 7, 1))

        plot_dict[country] = Plots.plot(title = country)
        plot_dict[country] = Plots.xaxis!("")
        plot_dict[country] = Plots.yaxis!("", :log10)

        xvar = timeModel2Real.(sol.t, country)
        yvar = calculateTotalDeaths(sol)
        plot_dict[country] = Plots.plot!(xvar, yvar, label = "")

        xvar = countryData[country][:cases][:t]
        yvar = countryData[country][:cases].deaths
        plot_dict[country] = Plots.scatter!(xvar, yvar, label = "", marker = :circle, markeralpha = 0.1)
    end

    list_plots = [plot for (_, plot) in plot_dict]
    vignette = Plots.plot(list_plots...)

    return vignette
end
