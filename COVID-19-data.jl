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
    filename = "https://github.com/neherlab/covid19_scenarios_data/raw/master/case-counts/World.tsv"
    download(filename, "data/cases_world.tsv")

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
    cases = DataFrame(CSV.read("data/cases_world.tsv", delim = "\t", header = 4))

    # Add a time column in the same format as the other dataframes
    cases = hcat(DataFrame(t = date2days.(cases[:, :time])), cases)

    # Remove any row with no recorded death
    cases = cases[cases.deaths.>0, :]

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

    return Dict(
        :codes => country_codes,
        :popData => popData,
        :cases => cases,
        :ICU => ICU_capacity,
        :Beds => hospital_capacity,
        :age_distribution => age_distribution,
    )
end

# Let's get the Neherlab repo's data (the files can be updated with `updateData()`).
COUNTRY_DB = loadData();

#-------------------------------------------------------------------------------------------------
#--
#-- Create a dictionary populated with the information of a given country
#--
function populateCountryDate(country, hemisphere)

    # This is the start of a country specific data structure. A dictionary is good enough for that purpose.
    countryData = Dict()

    countryData = Dict([(:name, country), (:R₀, BaseR₀), (:hemisphere, hemisphere)])

    countryData[:peak] = peakDate[countryData[:hemisphere]]
    countryData[:ϵ] = ϵ[countryData[:hemisphere]]
    countryData[:mitigation] = DEFAULT_MITIGATION

    # Clean up and extract country-specific information
    country_codes = @where(COUNTRY_DB[:codes], occursin.(country, :name))
    countryShort = country_codes[:, Symbol("alpha-3")][1]

    cases = @where(COUNTRY_DB[:cases], occursin.(country, :location))
    sort!(cases, :time)

    ICU_capacity = try
        @where(COUNTRY_DB[:popData], occursin.(country, :name))[!, :ICUBeds][1]
    catch
        0
    end
    ICU_capacity = convert(Float64, ICU_capacity)

    hospital_capacity = try
        @where(COUNTRY_DB[:popData], occursin.(country, :name))[!, :hospitalBeds][1]
    catch
        0
    end
    hospital_capacity = convert(Float64, hospital_capacity)

    age_distribution = @where(COUNTRY_DB[:age_distribution], occursin.(country, :_key))[!, 2:10]

    countryData[:country_code] = countryShort
    countryData[:cases] = cases
    countryData[:ICU_capacity] = ICU_capacity
    countryData[:hospital_capacity] = hospital_capacity
    countryData[:age_distribution] = age_distribution

    p, r = createDefaultParameters()
    countryData[:params] = p
    countryData[:range] = r

    return countryData
end
