# Utility funtions for use with the GTFSQuerys module

"""
    duration(inTime::String)

Converts a string with a time in the format "HH:MM:SS" to a time in seconds
"""
function duration(inTime::String)
    idx = [1,2,4,5,7,8]
    mult = [10*60*60,60*60,10*60,60,10,1]
    time::Int64 = 0

    for (i,m) in zip(idx,mult)
        time = time + parse(Int64,inTime[i])*m
    end

    return time
end


# Return unique files related with a vector of parameters
function get_files(q::GTFSQuery, params::Vector{Symbol})
    files = Vector{String}()
    for p in params
        append!(files, q.p2f[p])
    end
    return unique(files)
end

"""
    get_date_filters(date_str::String, 
                     times::Vector{String}=String[];
                     df::String="yyyymmdd",
                     day_names::Vector{String},
                     time_col::Vector{String},
                     date_col::Vector{String}
                    )

    Generates a `FilterDict` containing all filters required to apply a data and time filter to the GTFS dataset. 
    Optional arguments all for customimation of the time and date columns used for the filters as well as localization of the day names.
"""
function get_date_filters(date_str::String, 
                          times::Vector{String}=String[];
                          df::String="yyyymmdd",
                          day_names::Vector{Symbol}=[:monday,:tuesday,:wednesday,:thursday,:friday,:saturday,:sunday],
                          time_col::Vector{Symbol}=[:departure_time],
                          date_col::Vector{Symbol}=[:start_date,:end_date]
                          )
    # Compute date and day of week
    d = Date(date_str,DateFormat(df))
    d_num = Dates.dayofweek(d)
    str = Dates.format(d,dateformat"yyyymmdd")

    # Add filters for date
    filters = FilterDict(
        day_names[d_num]=> x -> .==(x, "1"),
        date_col[1]     => x -> .<=(x, str),
        date_col[2]     => x -> .>=(x, str)
    )

    # Add Time filters
    if !isempty(times)
        for col in time_col
            merge!(filters,
                    FilterDict(
                        col => x -> .>=(x, times[1]),
                        col => x -> .<(x, times[2])
                        )
                  )
        end
    end
    return filters
end