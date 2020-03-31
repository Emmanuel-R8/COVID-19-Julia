using Dates

#-------------------------------------------------------------------------------------------------
#--
#-- Linear interpolation of Y values for a new set of X values.
#--
function linearInterpolation(x, xvar, yvar)
    l = length(xvar)

    # If l = 1, ratio will be the only one
    if l == 1
        return yvar[1]
    else
        for i in 2:l
            x1 = xvar[i-1]
            x2 = xvar[i  ]

            if x < x2
                deltaY = yvar[i] - yvar[i-1]
                return yvar[i-1] + deltaY * (x - x1) / (x2 - x1)
            end
        end

        # Last possible choice
        return yvar[l]
    end
end

#-------------------------------------------------------------------------------------------------
#--
#-- Date functions
#--
function date2days(d)
    return convert(Float64, datetime2rata(d) - datetime2rata(BASE_DATE))
end

function days2date(d)
    return BASE_DATE + Day(d)
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


# Helper function to never change the number of individuals in a compartment in a way that would
# make it below 0.1 (to avoid rounding errors around 0)
function ensurePositive(d,s)
    return max.(d .+ s, 0.1) .- s
end



function plotCountriestoDisk(suffix)

    #--------------------------------------------------------------------------------------------------
    #-- Plot
    #--
    for (c, _) in COUNTRY_LIST
        global country = c

        clf();
        ioff();
        fig, ax = PyPlot.subplots();
        sol = calculateSolution(country,
                                DiseaseParameters,
                                countryData[country][:params];
                                finalDate = Date(2020, 6, 1))

        ax.plot(sol.t,
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
    end
end
