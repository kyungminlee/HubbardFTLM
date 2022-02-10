using DrWatson
@quickactivate "HubbardFTLM"

using HubbardFTLM
using PyPlot
using DataFrames
using Avro
using Arrow
using Quadmath
using PyPlot
using CSV 

using ArgParse

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "lattice"
            arg_type = String
            required = true
        "--out", "-o"
            arg_type = String
            required = false
            default = ""
        "--hopping", "-t"
            arg_type = Float64
            required = true
        "--interaction", "-U"
            arg_type = Float64
            required = true
        "--temperature", "-T"
            arg_type = Float64
            nargs = '+'
            required = true
    end
    parse_args(s)
end


function partition_generator(
        df_sectors::DataFrame,
        df_dense::DataFrame,
        df_sparse::DataFrame,
        charge::Integer;
        dense_overwrite::Bool=true
    )
    function partition(temperature::Real) 
        partitions = zeros(Float64, nrow(df_sectors))
        partitioncount = zeros(Int, nrow(df_sectors))
        partitiontype = fill(' ', nrow(df_sectors))
        for row in eachrow(df_sparse)
            z = sum(row.coefficients .* exp.(-row.shifted_eigenvalues / temperature)) * df_sectors.dim[row.idx] * row.krylovdim / length(row.shifted_eigenvalues)
            partitions[row.idx] += z
            partitioncount[row.idx] += 1
            partitiontype[row.idx] = 's'
        end
        for row in eachrow(df_dense)
            if partitiontype[row.idx] != ' ' && !dense_overwrite
                continue
            end
            z = sum(exp.(-row.shifted_eigenvalues / temperature))
            partitions[row.idx] = z
            partitioncount[row.idx] = 1
            partitiontype[row.idx] = 'd'
        end
        for row in eachrow(df_sectors)
            if row.dim == 0
                @assert iszero(partitions[row.idx])
                partitioncount[row.idx] = 1
                partitiontype[row.idx] = '0'
            elseif row.charge != charge
                partitions[row.idx] = 0
                partitioncount[row.idx] = 1
                partitiontype[row.idx] = 'x'
            elseif row.idx != row.root_idx
                partitions[row.idx] = partitions[row.root_idx]
                partitioncount[row.idx] = partitioncount[row.root_idx]
                partitiontype[row.idx] = partitiontype[row.root_idx]
            end
        end
        if !all((>)(0), partitioncount)
            @show charge
            @show collect(enumerate(zip(partitioncount, partitiontype, df_sectors.dim)))
        end
        @assert ' ' ∉ partitiontype
        @assert all((>)(0), partitioncount)
        return partitions ./ partitioncount, partitiontype
    end
end


function main()
    args = parse_commandline()
    latticetype, shape = parse_lattice(args["lattice"])
    hopping = args["hopping"]
    interaction = args["interaction"]

    lattice_str = lattice_string(latticetype, shape)
    n_sites = abs(shape[1,1] * shape[2,2] - shape[2,1] * shape[1,2])

    df_sectors = let
        buf = open(datadir("sectors-$lattice_str.arrow"), "r") do io
            read(io)
        end
        df = DataFrame(Arrow.Table(buf))
        df[!, :charge] = df.nup .+ df.ndn
        df[!, :density] = df.charge ./ n_sites
        df[!, :Sz] = 0.5 * (df.nup - df.ndn)
        df[!, :Sz2] = df.Sz.^2
        df
    end

    df_dense_all = open(datadir(lattice_str, "eigen", "dense-results.avro"), "r") do io
        DataFrame(Avro.readtable(io))
    end
    df_sparse_all = open(datadir(lattice_str, "eigen", "sparse-results.avro"), "r") do io
        DataFrame(Avro.readtable(io))
    end

    charges = unique(sort(df_sectors[:, :charge]))


    charges_out = Int64[]
    susceptibilities_out = Float64[]
    temperatures_out = Float64[]

    for charge in charges
        df_dense = filter(row->row.hopping == hopping && row.interaction == interaction && df_sectors[row.idx, :charge] == charge, df_dense_all)
        df_sparse = filter(row->row.hopping == hopping && row.interaction == interaction && df_sectors[row.idx, :charge] == charge, df_sparse_all)

        minimum_eigenvalue = Inf
        if !isempty(df_dense)
            minimum_eigenvalue = min(mapreduce(minimum, min, df_dense.eigenvalues))
        end
        if !isempty(df_sparse)
            minimum_eigenvalue = min(mapreduce(minimum, min, df_sparse.eigenvalues))
        end
        isfinite(minimum_eigenvalue) || continue

        if !isempty(df_dense)
            df_dense[!, :shifted_eigenvalues] = [row.eigenvalues .- minimum_eigenvalue for row in eachrow(df_dense)]
        end
        if !isempty(df_sparse)
            df_sparse[!, :shifted_eigenvalues] = [row.eigenvalues .- minimum_eigenvalue for row in eachrow(df_sparse)]
        end

        rho_f = partition_generator(df_sectors, df_dense, df_sparse, charge)
        temperatures = args["temperature"]

        for temperature in temperatures
            rho, rho_t = rho_f(temperature)
            fil = df_sectors.charge .== charge
            rho = rho[fil]
            sdf = df_sectors[fil, :]
            push!(charges_out, charge)
            push!(temperatures_out, temperature)
            push!(susceptibilities_out, sum(sdf.Sz2 .* rho) / sum(rho) / n_sites / temperature)
        end
    end
    df_out = DataFrame((charge=charges_out, temperature=temperatures_out, susceptibility=susceptibilities_out))
    if !isempty(args["out"])
        CSV.write(args["out"], df_out)
    else
        CSV.write(stdout, df_out)
    end
end

main()
