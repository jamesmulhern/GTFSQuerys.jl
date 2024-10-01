module GTFSQuerys

# Required packages
using CSV, DataFrames
using Graphs, MetaGraphsNext
using JuMP, HiGHS
using Dates

# Define Exports
export GTFSQuery, FilterDict
export query, duration, get_date_filters

# Type Alias
const FilterDict = Dict{Union{Symbol, Vector{Symbol}},Function}

# Define GTFSQuery object
mutable struct GTFSQuery
    folder_path::String
    g::MetaGraph
    f2p::Dict{String,Vector{Symbol}}
    p2f::Dict{Symbol,Vector{String}}
    cache::Dict{String,DataFrame}
    default_filters::FilterDict
end

# Includes
include("utilities.jl")

#
# Constructors
# 

"""
    GTFSQuery(folder_path::String; cache_files=[""])

Creates a `GTFSQuery` object based on GTFS files stored in `folder_path`.  Optionally specify a vector of files to cache
or ""all"" to cache all files.  
"""
function GTFSQuery(folder_path::String; cache_files=[""], default_filters::FilterDict=FilterDict())
    file_names = readdir(folder_path)

    # Create dataframe with parameter infromation
    param_ref = DataFrame([name => [] for name in ["v_index", "param", "parent_file"]])
    f2p = Dict{String,Vector{Symbol}}()
    for (i, file) in pairs(file_names)

        #Create file path
        file_path = joinpath([folder_path,file])
        trim_file = split(file,".")[1]

        # Open File
        f = open(file_path)
        param_names = split(readline(f),",")
        n_params = length(param_names)
        append!(param_ref, DataFrame(v_index=fill(i, n_params), 
                                        param=Symbol.(param_names), 
                                        parent_file=fill(file, n_params)
                                      )
                )
        close(f)

        # Add row to f2p dicttionary
        f2p[file] = Symbol.(param_names)
    end

    # Create param dict
    params = unique(param_ref.param)
    p2f = Dict{Symbol,Vector{String}}()

    for p in params
        rows = param_ref[:,:param] .== p
        p2f[p] = param_ref[rows,:parent_file]
    end
    
    # Create graph
    g = MetaGraph(
        Graph();
        label_type = String,            #File name wiht extension
        vertex_data_type = Int,         #File number
        edge_data_type = Symbol,        # Shared parameter    
        graph_data = "GTFS_File_Map")

    # Add vertexs
    for (i,k) in pairs(collect(keys(f2p)))
        g[k] = i
    end

    # Add edges
    for (param,f_vector) in pairs(p2f)
        if length(f_vector) <= 1
            continue
        end

        # Add all edges between files that share a param
        for (i,f1) in pairs(f_vector), (_,f2) in pairs(f_vector[i:end])
            g[f1,f2] = param
        end
    end

    # Build cache
    if !isempty(cache_files)
        cache = Dict{String,DataFrame}()

        if cache_files == "all"
            cache_files = file_names
        end

        for f_name in cache_files

            if !in(f_name, file_names)
                throw(ArgumentError("File $(f_name) is not in GTFS dataset"))
            end

            cache[f_name] = DataFrame(CSV.File(joinpath(folder_path,f_name),types=String))

            local_params = f2p[f_name]
            # Filter data
            active_filters = filter(((k,v),) -> k in local_params, default_filters)
            for (p, func) in active_filters
                subset!(cache[f_name], p => func, skipmissing=true)
            end
        end
    else
        cache = Dict{String,DataFrame}()
    end

    # Add extra edges
    #link = :stop_id => :from_stop_id
    #g["stops.txt","transfers.txt"] = :stop_id => :from_stop_id

    return GTFSQuery(folder_path, g, f2p, p2f, cache, default_filters)
end

#
# Graph Optimization Functions
#

"""
    solve_files_IP(q::GTFSQuery, params::Vector{Symbol})

Solves an integer program to determine the optimal set of files to open to read all parameters.
"""
function solve_files_IP(q::GTFSQuery, params::Vector{Symbol})
    
    files = get_files(q, params)
    n_params = length(params)
    n_files = length(files)

    mat = zeros(n_params,n_files)
    for (i,p) in pairs(params)
        sub_files = q.p2f[p]
        for f in sub_files
            j = findfirst(files .== f)
            mat[i,j] = 1
        end
    end

    model = Model(HiGHS.Optimizer)
    # Set Model Parameters
    set_attribute(model, "output_flag", false)
    set_attribute(model, "presolve", "on")
    set_attribute(model, "time_limit", 300.0)

    @variable(model, x[1:n_files], Bin)

    b = ones(Int, n_params)

    @constraint(model, mat * x .>= b)
    @objective(model, Min, sum(x))  # TODO add weighting to files. Number of lines?
    optimize!(model)

    f_idx = findall(trunc.(Bool, value.(x)))
    return files[f_idx]
end

# Computes a sub-graph that connects all of the required files and computes a shortest path tree between all the files
# TODO can I skip the shortest paths and just build the parents dit based on the stiener tree?
function compute_links(q::GTFSQuery, files::Vector{String}, root_param::Symbol)
    if length(files) == 1
        src_f = files[1]
        parents = Dict{String,String}(src_f => "")
        extra_params = Vector{Symbol}()

        return src_f, parents, extra_params
    end

    # Find Steiner Tree connecting all files
    st = steiner_tree(q.g, [code_for(q.g, f) for f in files])

    # Find source file
    src_options = q.p2f[root_param]
    valid = [f in files for f in src_options]
    idx = findfirst(valid)
    src_f = src_options[idx]

    # Compute shortest paths
    ds = dijkstra_shortest_paths(st, code_for(q.g,src_f))
    p_idx = ds.parents

    # Determine extra paramters that are required
    extra_params = Vector{Symbol}()
    for e in edges(st)
        push!(extra_params, q.g[label_for(q.g, src(e)),label_for(q.g, dst(e))])
    end
    unique!(extra_params)

    parents = Dict{String,String}() # Child -> parent
    for (i, p_code) in pairs(p_idx)
        if p_code != 0
            child = label_for(q.g,i)
            parent = label_for(q.g,p_code)
            parents[child] = parent
            if !haskey(parents,parent)
                parents[label_for(q.g,p_code)] = ""
            end
        end
    end

    return src_f, parents, extra_params
end

#
# Data Access Functions
#
"""
    get_data(q::GTFSQuery, cur_f::String, params::Vector{Symbol}, filters::FilterDict)

Returns a view of the file dataframes that exist in the cache.  If the files does not exist in cache it loads the file from disk.
Outputs are filtered usings the filters provided and are limited to the parameters passed in the `params` argument.  
"""
function get_data(q::GTFSQuery, cur_f::String, params::Vector{Symbol}, filters::FilterDict)

    # Check if file is in cache
    if haskey(q.cache,cur_f)
        df = view(q.cache[cur_f], :, params)
    else
        # Load from CSV and store in cache
        q.cache[cur_f] = DataFrame(CSV.File(joinpath(q.folder_path,cur_f), types=String))
        df = view(q.cache[cur_f], :, params)
        #df = DataFrame(CSV.File(joinpath(q.folder_path,cur_f),types=String))
    end

    # Apply filters to data
    for (param, func) in filters
        df = subset(df, param => func, view=true, skipmissing=true)
    end
    return df
end


"""
    read_files(q::GTFSQuery, cur_f::String, parents::Dict{String,String}, params::Vector{Symbol}, filters::FilterDict)

Recusivly reads data from the GTFS data starting at `cur_f`.  Function is called recersily based on the provided `parents` dictionary.
Parameters are limited to what is provided in `params` and output is filtered by `filters`
"""
function read_files(q::GTFSQuery, 
                    cur_f::String, 
                    parents::Dict{String,String}, 
                    params::Vector{Symbol}, 
                    filters::FilterDict)

    # Get params in current file
    local_params = [p for p in params if p in q.f2p[cur_f]]

    # Read file
    df = get_data(q, cur_f, local_params, filter(((k,v),) -> all(in.(k,[local_params])), filters))

    # Find parent of current file and shared parameter
    p_file = parents[cur_f]
    if p_file != ""
        share_param = q.g[p_file,cur_f]
    else
        share_param = Symbol()
    end

    # Add data to dict
    data = Dict{String, Tuple{SubDataFrame, String, Symbol}}(cur_f=>(df,p_file,share_param))

    # Find children
    children = Vector{String}()
    for (child,parent) in pairs(parents)
        if parent == cur_f
            push!(children,child)
        end
    end

    # Recursive Calls
    if !isempty(children)
        for child in children
            sub_data = read_files(q, child, parents, params, filters)
            merge!(data,sub_data)
        end
    end

    return data
end


"""
    query(q::GTFSQuery, params::Vector{Symbol}, filters::FilterDict; default_filters::Bool=true)

Performs a query of the GTFS data.  Extracts all values of provided `params` with all filters applied.  
Use `default_filters` to activate stored default filters.
"""
function query(q::GTFSQuery, 
               params::Vector{Symbol}, 
               filters::FilterDict;
               default_filters::Bool=true)
    
    # Incoperate default filters if active
    if default_filters
        merge!(filters,q.default_filters)
    end
    
    # Get filter parameters
    #filter_params = collect(keys(filters))
    filter_params = Vector{Symbol}()
    for k in keys(filters)
        if isa(k, Symbol)
            push!(filter_params, k)
        else
            append!(filter_params, k)
        end
    end
    all_params = unique(vcat(params, filter_params))

    # Check all params exist
    # if !all(in(p,collect(keys(q.p2f))) for p in all_params)
    #     error("Not all parameters exist in GTFS dataset")
    # end

    #Check if all parameters are in the dataset
    for p in all_params
        if in(p, collect(keys(q.p2f)))
            continue
        end

        #Throw error
        throw(ArgumentError("Paramter $p is not a memeber of GTFS dataset"))
    end

    # Solve files IP
    opt_files = solve_files_IP(q, all_params)

    # Determine the best linking of files
    src_f, parents, extra_params = compute_links(q,opt_files,params[1])
    unique!(append!(all_params, extra_params))

    df_dict = read_files(q, src_f, parents, all_params, filters)

    return unique(select(merge_df_dict(df_dict),params))
end





function merge_df_dict(df_dict::Dict{String,Tuple{SubDataFrame,String,Symbol}},root_f::String="")
    
    # Find key for root file
    if isempty(root_f)
        for (k, val) in pairs(df_dict)
            if val[2] == ""
                root_f = k
            end
        end
    end

    df = df_dict[root_f][1]

    # Find children
    children = String[]
    for (k, val) in pairs(df_dict)
        if val[2] == root_f
            push!(children,k)
        end
    end

    # Loop over children nodes
    for child in children
        sub_df = merge_df_dict(df_dict, child)
        df = innerjoin(df, sub_df, on=df_dict[child][3], makeunique=true, matchmissing=:notequal)
    end

    return df
end

end # module GTFSQuerys
