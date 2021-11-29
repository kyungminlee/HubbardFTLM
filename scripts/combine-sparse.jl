using DrWatson
@quickactivate "HubbardFTLM"

using HubbardFTLM
using DataFrames
using ProgressMeter

using ArgParse
using Glob
using HDF5
using JSON3
using Arrow

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
    sparse_filepath = datadir(lattice_str, "eigen", "eigen-sparse-results.hdf5")
    qedpathlist = Glob.glob("eigen-sparse-results*.qed", datadir(lattice_str, "eigen", "sparsedata"))

    IdentifierType = NamedTuple{(:idx, :hopping, :interaction), Tuple{Int, Float64, Float64}}
    existing_results = Set(IdentifierType[])
    sparse_hdf5 = h5open(sparse_filepath, "cw")
    if !haskey(sparse_hdf5, "eigen-sparse") 
        create_group(sparse_hdf5, "eigen-sparse")
    end

    h5g_eigen = sparse_hdf5["eigen-sparse"]

    @info "Read $sparse_filepath with $(length(existing_results)) parameter sets."
    for qedpath in qedpathlist
        # df = CSV.read(filename, DataFrame)
        parameter = open(joinpath(qedpath, "parameter.json"), "r") do io
            JSON3.read(io; jsonlines=false, allow_inf=false)
        end

        idx::Int = parameter["idx"]
        hopping::Float64 = parameter["hopping"]
        interaction::Float64 = parameter["interaction"]
        krylovdim::Int = parameter["krylovdim"]
        seed::Int = parameter["seed"]

        eigenvalues = Float64[]
        for eigenvalue_filepath = Glob.glob("eigenvalue_*.arrow", qedpath)
            eigenvalue_table = Arrow.Table(eigenvalue_filepath)
            append!(eigenvalues, eigenvalue_table.eigenvalue)
        end

        group_name = savename((idx=idx, hopping=hopping, interaction=interaction))

        if haskey(h5g_eigen, group_name)
            HDF5.delete_object(h5g_eigen, group_name)
        end

        group = create_group(h5g_eigen, group_name)
        attributes(group)["idx"] = idx
        attributes(group)["hopping"] = hopping
        attributes(group)["interaction"] = interaction
        #attributes(group)["githash"] = githash
        #attributes(group)["timestamp"] = timestamp
        attributes(group)["type"] = "sparse"
        group["eigenvalue"] = eigenvalues
        flush(sparse_hdf5)

        @info "Appended $qedpath."
    end
    @info "Total $(length(h5g_eigen)) parameter sets"


    #=
    existing_results = Set(IdentifierType[])
    delete_list = String[]
    for c in h5g_eigen
        idx = read(attributes(c)["idx"])
        t = read(attributes(c)["hopping"])
        U = read(attributes(c)["interaction"])
        identifier = (idx=idx, hopping=t, interaction=U)
        if identifier in existing_results
            push!(delete_list, HDF5.name(c))
        end
        push!(existing_results, identifier)
    end

    @info "Delete list (not implemented yet)" delete_list
    =#

    close(sparse_hdf5)
    #for filename in filenames
    #    mv(filename, "$filename.processed")
    #end
end

main()
