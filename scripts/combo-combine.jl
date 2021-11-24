using DrWatson
@quickactivate "Hubbard"

using Hubbard
using DataFrames
using CSV
using ProgressMeter

using ArgParse
using Glob
using HDF5
using JSON

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "lattice"
            arg_type = String
            help = "shape of the lattice in the format lattice-(?,?)x(?,?)"
            required = true
    end
    parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    # latticetype, shape = parse_lattice(parsed_args["lattice"])
    lattice_str = parsed_args["lattice"]
    dense_filepath = datadir(lattice_str, "eigen", "eigen-dense-results.hdf5")
    filenames = Glob.glob("temp-eigen-dense-*.json", datadir(lattice_str, "eigen"))

    IdentifierType = NamedTuple{(:idx, :hopping, :interaction), Tuple{Int, Float64, Float64}}
    existing_results = Set(IdentifierType[])
    dense_hdf5 = h5open(dense_filepath, "cw")
    if !haskey(dense_hdf5, "eigen")
        create_group(dense_hdf5, "eigen")
    end
    h5g_eigen = dense_hdf5["eigen"]

    for c in h5g_eigen
        idx = read(attributes(c)["idx"])
        t = read(attributes(c)["hopping"])
        U = read(attributes(c)["interaction"])
        push!(existing_results, (idx=idx, hopping=t, interaction=U))
    end

    @info "Read $dense_filepath with $(length(existing_results)) parameter sets."
    for filename in filenames
        # df = CSV.read(filename, DataFrame)
        data = open(filename, "r") do io
            JSON.parse(io; dicttype=Dict, allownan=false)
        end
        for item in data
            idx::Int, hopping::Float64, interaction::Float64 = item["idx"], item["hopping"], item["interaction"]
            eigenvalues::Vector{Float64} = [item["eigenvalue"]...]
            identifier = (idx=idx, hopping=hopping, interaction=interaction)
            if identifier in existing_results
                continue
            end
            group_name = savename((idx=idx, hopping=hopping, interaction=interaction))
            group = create_group(h5g_eigen, group_name)
            attributes(group)["idx"] = idx
            attributes(group)["hopping"] = hopping
            attributes(group)["interaction"] = interaction
            group["eigenvalue"] = eigenvalues
            flush(dense_hdf5)
            push!(existing_results, identifier)
        end
        @info "Appended $filename."
    end
    @info "Total $(length(h5g_eigen)) parameter sets"
    close(dense_hdf5)
    for filename in filenames
        mv(filename, "$filename.processed")
    end
end

main()