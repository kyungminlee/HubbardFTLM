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
        df_sparse::DataFrame;
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
            elseif row.idx != row.root_idx
                partitions[row.idx] = partitions[row.root_idx]
                partitioncount[row.idx] = partitioncount[row.root_idx]
                partitiontype[row.idx] = partitiontype[row.root_idx]
            end
        end
        @assert all((>)(0), partitioncount)
        @assert ' ' ∉ partitiontype
        return partitions ./ partitioncount, partitiontype
    end
end


function main()
    args = parse_commandline()
    latticetype, shape = parse_lattice(args["lattice"])
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

    df_dense = open(datadir(lattice_str, "eigen", "dense-results.avro"), "r") do io
        DataFrame(Avro.readtable(io))
    end
    df_sparse = open(datadir(lattice_str, "eigen", "sparse-results.avro"), "r") do io
        DataFrame(Avro.readtable(io))
    end

    minimum_eigenvalue = min(mapreduce(minimum, min, df_dense.eigenvalues), mapreduce(minimum, min, df_sparse.eigenvalues))

    df_dense[!, :shifted_eigenvalues] = [row.eigenvalues .- minimum_eigenvalue for row in eachrow(df_dense)]
    df_sparse[!, :shifted_eigenvalues] = [row.eigenvalues .- minimum_eigenvalue for row in eachrow(df_sparse)]

    ρf = partition_generator(df_sectors, df_dense, df_sparse)
    temperatures = args["temperature"]

    susceptibilities = Float64[]
    for temperature in temperatures
        fil = df_sectors.charge .== 9
        ρ, ρt = ρf(temperature)
        ρ = ρ[fil]
        sdf = df_sectors[fil, :]
        push!(susceptibilities, sum(sdf.Sz2 .* ρ) / sum(ρ) / n_sites / temperature)
    end
    df_out = DataFrame((temperature=temperatures, susceptibility=susceptibilities))
    if !isempty(args["out"])
        CSV.write(args["out"], df_out)
    else
        CSV.write(stdout, df_out)
    end
end

main()