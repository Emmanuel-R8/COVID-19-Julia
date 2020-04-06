#-------------------------------------------------------------------------------------------------
#--
#-- Ususal libraries
#--
using CSV, Dates
using DataFrames, DataFramesMeta

#-------------------------------------------------------------------------------------------------
#--
#-- DATA
#--
# Update latest data from Neherlab
function updateData()

    ### COUNTRY CODES
    filename = "https://github.com/neherlab/covid19_scenarios_data/raw/master/country_codes.csv"
    download(filename, "data/country_codes.csv")

    ### LIST OF CASES
    filename = "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
    download(filename, "data/ecdc.csv")

    ### ICU beds
    filename = "https://github.com/neherlab/covid19_scenarios_data/raw/master/hospital-data/ICU_capacity.tsv"
    download(filename, "data/ICU_capacity.tsv")

    ### Hospital beds
    filename = "https://github.com/neherlab/covid19_scenarios_data/raw/master/hospital-data/hospital_capacity.csv"
    download(filename, "data/hospital_capacity.csv")
end


function loadData()
    ### Population Data
    popData = CSV.read("data/populationData.tsv", delim = "\t")

    ### COUNTRY CODES
    country_codes = select(CSV.read("data/country_codes.csv"), :name, Symbol("alpha-3"))

    ### LIST OF CASES
    ecdc = CSV.read("data/ecdc.csv")

    ### ICU beds
    ICU_capacity = select(CSV.read("data/ICU_capacity.tsv", delim = "\t"), :country, :CriticalCare)

    ### Hospital beds
    hospital_capacity = select(
        CSV.read("data/hospital_capacity.csv", types = Dict(:COUNTRY => String), limit = 1267),
        :COUNTRY,
        :YEAR,
        :VALUE,
    )
    hospital_capacity = @where(hospital_capacity, Not(ismissing.(:COUNTRY)))

    ### Age pyramids
    age_distribution = CSV.read("data/country_age_distribution.csv")

    ### Last saved parameters
    optimisedParameters = try
        CSV.read("data/allCountryParameters.csv")
    catch
        nothing
    end

    return Dict(
        :codes => country_codes,
        :popData => popData,
        :cases => ecdc,
        :ICU => ICU_capacity,
        :Beds => hospital_capacity,
        :age_distribution => age_distribution,
        :optimisedParameters => optimisedParameters
    )
end

# Let's get the Neherlab repo's data (the files can be updated with `updateData()`).
COUNTRY_DB = loadData()

#-------------------------------------------------------------------------------------------------
#--
#-- Create a dictionary populated with the information of a given country
#--
function populateCountryDate(country, hemisphere; useOptimised = true)

    # This is the start of a country specific data structure. A dictionary is good enough for that purpose.
    countryData = Dict()

    countryData = Dict([(:name, country), (:hemisphere, hemisphere)])

    countryData[:peak] = peakDate[countryData[:hemisphere]]
    countryData[:ϵ] = ϵ[countryData[:hemisphere]]
    countryData[:mitigation] = DEFAULT_MITIGATION


    # Cases from the ECDC database
    caseDB = COUNTRY_DB[:cases]

    # select the right country
    ecdc = @where(caseDB, country .== :countriesAndTerritories)

    # Add a date comlumn in the Date type and convert to days
    ecdc[!, :time] = Date.(ecdc.year, ecdc.month, ecdc.day)
    ecdc[!, :t] = date2days.(ecdc.time)

    # and order all dates in ascending order
    ecdc = @orderby(ecdc, :time)

    # Convert cases and deaths from daily deaths to cumulative amounts
    ecdc.cases = cumsum(ecdc.cases)
    ecdc.deaths = cumsum(ecdc.deaths)

    # Only retain entries where there is a strictily positive number of deaths
    ecdc = @where(ecdc, :deaths .> 0)

    # Extract country-specific information
    countryShort = ecdc[1, :countryterritoryCode]

    population = ecdc[1, :popData2018]

    ICU_capacity = try
        @where(COUNTRY_DB[:popData], occursin.(country, :name))[!, :ICUBeds][1]
    catch
        population / 10.0
    end
    ICU_capacity = convert(Float64, ICU_capacity)

    hospital_capacity = try
        @where(COUNTRY_DB[:popData], occursin.(country, :name))[!, :hospitalBeds][1]
    catch
        population / 100.0
    end
    hospital_capacity = convert(Float64, hospital_capacity)

    age_distribution = @where(COUNTRY_DB[:age_distribution], occursin.(country, :_key))[1, 2:10]

    countryData[:population] = population
    countryData[:country_code] = countryShort
    countryData[:cases] = ecdc
    countryData[:ICU_capacity] = ICU_capacity
    countryData[:hospital_capacity] = hospital_capacity
    countryData[:age_distribution] = age_distribution

    countryData[:lossFunction] = p -> singleCountryLoss(country, p)

    if useOptimised == true
        p =  @where(COUNTRY_DB[:optimisedParameters], occursin.(country, :country))[:, 3:end]
        if nrow(p) == 1
            countryData[:params] =  convert(Array, p)
        end
    else
        countryData[:params] = COUNTRY_INIT
        countryData[:range] = COUNTRY_RANGE
    end

    return countryData
end


# Update the record of cases for countries after an update
function updateCountryCases(country)

    # Cases from the ECDC database
    caseDB = COUNTRY_DB[:cases]

    # select the right country
    ecdc = @where(caseDB, country .== :countriesAndTerritories)

    # Add a date comlumn in the Date type and convert to days
    ecdc[!, :time] = Date.(ecdc.year, ecdc.month, ecdc.day)
    ecdc[!, :t] = date2days.(ecdc.time)

    # and order all dates in ascending order
    ecdc = @orderby(ecdc, :time)

    # Convert cases and deaths from daily deaths to cumulative amounts
    ecdc.cases = cumsum(ecdc.cases)
    ecdc.deaths = cumsum(ecdc.deaths)

    # Only retain entries where there is a strictily positive number of deaths
    ecdc = @where(ecdc, :deaths .> 0)

    global countryData[country][:cases] = ecdc
end


#--------------------------------------------------------------------------------------------------
#--
#-- Find the date at which cases exceed 5 for a given country
#--
function approximateModelStartRange(country::String)
    # Estimates when deaths number is about DEATH_AT_MODEL_START
    above = @where(countryData[country][:cases], :deaths .> DEATH_AT_MODEL_START)

    # Take first date
    if nrow(above) > 0
        first_date_above = first(above[:, :time])
    else
        first_date_above = BASE_DATE
    end

    # Convert to day
    above_day = date2days(first_date_above)

    return(above_day - 10.0, above_day + 2.0)
end
