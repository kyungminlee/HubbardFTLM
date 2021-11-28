using DrWatson
@quickactivate "HubbardFTLM"

using HubbardFTLM

using LinearAlgebra
using Random

using LatticeTools
using KrylovKit
using QuantumHamiltonian
using QuantumHamiltonianParticle
using FiniteTemperatureLanczos

using DataStructures
using Formatting

using Logging
using ArgParse
using ProgressMeter

using HDF5
using JSON3
using Arrow
using DataFrames
using Serialization

using UUIDs

try
  @eval using MKL
catch
end

const SectorType = NamedTuple{(:idx, :nup, :ndn, :tii, :pii, :pic), Tuple{Int, Int, Int, Int, Int, Int}}


function process_sparse(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
    sector_indices::AbstractVector{<:Integer};
    seed::Integer=0,
    krylovdim::Integer=200,
    blocksize::Integer=10,
    nsamples::Integer=1000,
)
    @mylogmsg "Collecting sectors"
    lattice_str = lattice_string(latticetype, shape)
    sectors_filepath = datadir("sectors-$lattice_str.arrow")

    select_table = DataFrame(idx=sector_indices)
    sectors_table = DataFrame(Arrow.Table(sectors_filepath))
    @mylogmsg "Finished reading sectors"

    table = innerjoin(select_table, sectors_table, on=:idx)
    sectors = SectorType[]
    for row in eachrow(table)
        push!(sectors, (idx=row.idx, nup=row.nup, ndn=row.ndn, tii=row.tii, pii=row.pii, pic=row.pic))
    end
    return process_sparse(
        latticetype, shape, t, U, sectors;
        seed=seed, krylovdim=krylovdim, blocksize=blocksize, nsamples=nsamples
    )
end


function getbytes(rng::MersenneTwister)
    io = IOBuffer()
    serialize(io, rng)
    seekstart(io)
    buf = read(io)
    return buf
end


# function write_status!(status_group::HDF5.Group, samplecount::Integer, rng::MersenneTwister)
#     if haskey(status_group, "samplecount")
#         delete_object(status_group, "samplecount")
#     end
#     write_dataset(status_group, "samplecount", samplecount)
#     io = IOBuffer()
#     serialize(io, rng)
#     seekstart(io)
#     buf = read(io)
#     if haskey(status_group, "random_number_generator")
#         rgen = status_group["random_number_generator"]
#         rgen[:] .= zero(UInt8)
#         rgen[1:length(buf)] = buf
#     else
#         n = length(buf)
#         # m = 2^Int(ceil(log2(n+1024)))
#         m = Int(ceil((n+512)/1024))*1024
#         resize!(buf, m)
#         buf[n+1:end] .= zero(UInt8)
#         status_group["random_number_generator"] = buf
#     end
#     close(io)
# end


function process_sparse(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
    sectors::AbstractVector{SectorType}
    ;
    seed::Integer=0,
    krylovdim::Integer=200,
    blocksize::Integer=10,
    nsamples::Integer=1000,
)
    @mylogmsg "Start process_sparse"
    BR = UInt

    # Set up lattice
    if latticetype == "triangular"
        mylattice = make_triangular_lattice(shape)
    elseif latticetype == "square"
        mylattice = make_square_lattice(shape)
    else
        @mylogmsg "Unsupported lattice type $latticetype"
        exit(1)
    end

    lattice = mylattice.lattice
    ssymbed = mylattice.space_symmetry_embedding
    tsymbed = ssymbed.normal
    psymbed = ssymbed.rest

    hs, hamiltonian, c, cdag = make_hamiltonian_ladder(mylattice, t, U)
    @mylogmsg "Finished making Hamiltonian"

    pair_groups = Vector{Tuple{Int, Int}}[]
    let visited = falses(numsites(lattice.supercell), numsites(lattice.supercell))
        for (i1, (sitename1, sitefc1)) in enumerate(lattice.supercell.sites),
            (i2, (sitename2, sitefc2)) in enumerate(lattice.supercell.sites)
            visited[i1, i2] && continue
            pair_group = Set()
            for t in tsymbed.elements, p in psymbed.elements
                j1 = p(t(i1))
                j2 = p(t(i2))
                push!(pair_group, (j1, j2))
                visited[j1, j2] = true
            end
            push!(pair_groups, sort(collect(pair_group)))
        end
    end

    nup(i::Integer) = cdag(i, :up) * c(i, :up)
    ndn(i::Integer) = cdag(i, :dn) * c(i, :dn)

    correlation = []

    @mylogmsg "pair groups: $(pair_groups)"

    for pair_group in pair_groups
        push!(correlation, sum(nup(i) * nup(j) for (i, j) in pair_group) / length(pair_group))
        push!(correlation, sum(nup(i) * ndn(j) for (i, j) in pair_group) / length(pair_group))
        push!(correlation, sum(ndn(i) * nup(j) for (i, j) in pair_group) / length(pair_group))
        push!(correlation, sum(ndn(i) * ndn(j) for (i, j) in pair_group) / length(pair_group))
    end

    @assert all(isinvariant(x, y) for x in tsymbed.elements for y in correlation)
    @assert all(isinvariant(x, y) for x in psymbed.elements for y in correlation)

    correlation = [make_projector_operator(op, BR) for op in correlation]

    prev_idx = nothing
    rhsr = nothing

    prev_nupndn = nothing
    hsr = nothing

    prev_tii = nothing
    tsa = nothing
    psymbed_little = nothing

    prev_tiipiipic = nothing
    psa = nothing
    ssa = nothing

    lattice_str = lattice_string(latticetype, shape)
    githash = HubbardFTLM.getgithash()

    @mylogmsg "Starting"
    for (idx, nup, ndn, tii, pii, pic) in sectors
        @mylogmsg "Sector: idx=$idx, nup=$nup, ndn=$ndn, tii=$tii, pii=$pii, pic=$pic"

        output_groupname = savename(Dict(:idx=>idx, :t=>t, :U=>U, :krylovdim=>krylovdim, :seed=>seed))
        output_filename = savename("eigen-sparse-results", Dict(:idx=>idx, :t=>t, :U=>U, :krylovdim=>krylovdim, :seed=>seed), "qed")
        output_filepath = datadir(lattice_str, "eigen", "sparsedata", output_filename)

        @mylogmsg "File: $output_filepath"
        # isdir(dirname(output_filepath)) || mkpath(dirname(output_filepath))
        # output_file = h5open(output_filepath, "cw")
        isdir(output_filepath) || mkpath(output_filepath)
        
        if !isfile(joinpath(output_filepath, "parameter.json"))
            open(joinpath(output_filepath, "parameter.json"), "w") do io
                JSON3.write(io, (idx=idx, hopping=t, interaction=U, krylovdim=krylovdim, seed=seed))
            end
        end

        if isfile(joinpath(output_filepath, "status"))
            @mylogmsg "Reading existing random number generator"
            status = open(joinpath(output_filepath, "status"), "r") do io
                deserialize(io)
            end
            samplecount = status.samplecount
            rng = status.random_number_generator
        else
            @mylogmsg "Starting afresh"
            rng = MersenneTwister(seed)
            samplecount = 0
            status = (samplecount=samplecount, random_number_generator=rng)
            open(joinpath(output_filepath, "status"), "w") do io
                serialize(io, status)
            end
        end

        @mylogmsg "rng: $rng"
        @mylogmsg "samplecount: $samplecount"
        if samplecount >= nsamples
            @mylogmsg "samplecount $samplecount ≥ nsamples $nsamples"
            continue
        end
        GC.gc()
        @mylogmsg "Garbage collected"

        if idx == prev_idx
            @mylogmsg "Using previous sector rhsr"
            @assert !isnothing(rhsr)
        else
            if (nup, ndn) == prev_nupndn
                @mylogmsg "Using previous hsr"
                @assert !isnothing(hsr)
            else
                hsr = represent(hs, [(nup, ndn),], BR, DictIndexedVector)
                prev_nupndn = (nup, ndn)
                @mylogmsg "Finished creating hsr"
            end
            if tii == prev_tii
                @mylogmsg "Using previous tii"
                @assert !isnothing(tsa)
                @assert !isnothing(psymbed_little)
            else
                tsa = collect(get_irrep_iterator(IrrepComponent(tsymbed, tii)))
                psymbed_little = little_symmetry(tsymbed, tii, psymbed)
                prev_tii = tii
            end
            if (tii, pii, pic) == prev_tiipiipic
                @mylogmsg "Using previous psa, ssa"
                @assert !isnothing(psa)
                @assert !isnothing(ssa)
            else
                psa = collect(get_irrep_iterator(IrrepComponent(psymbed_little, pii, pic)))
                ssa = make_product_irrep(psa, tsa)
                prev_tiipiipic = (tii, pii, pic)
            end
            rhsr = symmetry_reduce(hsr, ssa)
            @mylogmsg "Finished creating rhsr"
        end
        D = dimension(rhsr)
        if D == 0
            @mylogmsg "Sector has zero dimension. Skip."
            continue
        end
        @mylogmsg "Hilbert space dimension: $D"

        @mylogmsg "Generating Hamiltonian representation"
        H = represent(rhsr, hamiltonian)
        for r in 1:blocksize
            @mylogmsg "Start computing"
            if samplecount >= nsamples
                @mylogmsg "samplecount $samplecount ≥ nsamples $nsamples"
                break
            end
            samplecount += 1
    
            v = randn(rng, ComplexF64, D)
            normalize!(v)

            iterator = LanczosIterator(x -> H*x, v)
            factorization = initialize(iterator)
            sizehint!(factorization, krylovdim)
            @mylogmsg "Finished sizehint!"
            for d in 1:krylovdim-1
                expand!(iterator, factorization)
                # @mylogmsg "Expanded factorization to $d"
            end
            @mylogmsg "Garbage collecting after Krylov"
            GC.gc()
            hmat = SymTridiagonal(factorization.αs, factorization.βs)
            eigenvalues = eigvals!(hmat)
            Arrow.write(joinpath(output_filepath, "eigenvalue_$samplecount.arrow"), (eigenvalue=eigenvalues,))

            open(joinpath(output_filepath, "status"), "w") do io
                serialize(io, (samplecount=samplecount, random_number_generator=rng))
            end
            @mylogmsg "Wrote $samplecount"
        end
    end # for (idx, nup, ndn, tii, pii, pic) in sectors
end


function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "lattice"
            arg_type = String
            help = "shape of the lattice in the format type-(?,?)x(?,?)"
            required = true
        "--hopping", "-t"
            arg_type = Float64
            default = 1.0
        "--interaction", "-U"
            arg_type = Float64
            required = true
        "--sector"
            arg_type = Int
            nargs = '+'
            required = true
        "--krylovdim"
            arg_type = Int
            range_tester = x -> x > 0
            default = 100
        "--blocksize"
            arg_type = Int
            range_tester = x -> x > 0
            default = 10
        "--nsamples"
            arg_type = Int
            range_tester = x -> x > 0
            default = 100
        "--seed"
            arg_type = Int
            default = 0
        "--force", "-f"
            action = :store_true
        "--debug", "-d"
            action = :store_true
    end
    parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    if parsed_args["debug"]
        logger = ConsoleLogger(stdout, Logging.Debug; meta_formatter=my_metafmt)
        global_logger(logger)
    else
        logger = ConsoleLogger(stdout, Logging.Info; meta_formatter=my_metafmt)
        global_logger(logger)
    end

    latticetype, shape = parse_lattice(parsed_args["lattice"])
    t = parsed_args["hopping"]
    U = parsed_args["interaction"]
    seed = parsed_args["seed"]
    krylovdim = parsed_args["krylovdim"]
    blocksize = parsed_args["blocksize"]
    nsamples = parsed_args["nsamples"]

    sectors = parsed_args["sector"]


    @mylogmsg "Starting $(first(ARGS))"
    @mylogmsg "shape: $shape"
    @mylogmsg "t: $t"
    @mylogmsg "U: $U"
    @mylogmsg "seed: $seed"
    @mylogmsg "krylovdim: $krylovdim"
    @mylogmsg "blocksize: $blocksize"
    @mylogmsg "nsamples: $nsamples"
    @mylogmsg "sectors: $sectors"

    process_sparse(latticetype, shape, t, U, sectors; seed=seed, krylovdim=krylovdim, blocksize=blocksize, nsamples=nsamples)

    @mylogmsg "Finished writing"
end


main()

