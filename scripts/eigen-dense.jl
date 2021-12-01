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
using Dates

using DataFrames
using JSON3
using Arrow

using UUIDs

using ProgressMeter

try
  @eval using MKL
catch
end

const SectorType = NamedTuple{(:idx, :nup, :ndn, :tii, :pii, :pic), Tuple{Int, Int, Int, Int, Int, Int}}


function eigen_small(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
)
    lattice_str = lattice_string(latticetype, shape)
    schedule_table = open(datadir(lattice_str, "schedule.arrow"), "r") do io
        Arrow.Table(io)
    end
    sectors::Vector{Int} = schedule_table.idx[schedule_table.type .== "small"]
    return eigen_dense(latticetype, shape, t, U, sectors)
end


function eigen_dense(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
    sector_indices::AbstractVector{<:Integer},
)
    lattice_str = lattice_string(latticetype, shape)
    sectors_table = open(datadir("sectors-$lattice_str.arrow"), "r") do io
        buf = read(io)
        DataFrame(Arrow.Table(buf))
    end
    select_table = DataFrame(idx=sector_indices)

    table = innerjoin(select_table, sectors_table, on=:idx)
    sectors = SectorType[]
    for row in eachrow(table)
        push!(sectors, (idx=row.idx, nup=row.nup, ndn=row.ndn, tii=row.tii, pii=row.pii, pic=row.pic))
    end
    return eigen_dense(latticetype, shape, t, U, sectors)
end


function eigen_dense(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
    sectors::AbstractVector{SectorType},
)
    BR = UInt

    # Set up lattice
    if latticetype == "triangular"
        mylattice = make_triangular_lattice(shape)
    elseif latticetype == "square"
        mylattice = make_square_lattice(shape)
    else
        @error "Unsupported lattice type $latticetype"
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

    @info "pair groups: $(pair_groups)"

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
    dense_filepath = datadir(lattice_str, "eigen", "eigen-dense-results.avro") # result files are avro
    temp_filepath = datadir(lattice_str, "eigen", "temp-eigen-dense-$(uuid5(uuid1(), gethostname())).jsonl")

    IdentifierType = NamedTuple{(:hopping, :interaction, :idx), Tuple{Float64, Float64, Int}}
    existing_results = Set(IdentifierType[])
    if isfile(dense_filepath)
        let table = open(dense_filepath, "r") do io
                Avro.readtable(io)
            end
            for row in eachrow(table)
                push!(existing_results, (hopping=row.hopping, interaction=row.interaction, idx=row.idx))
            end
        end
    end

    isdir(dirname(temp_filepath)) || mkpath(dirname(temp_filepath))
    output_file = open(temp_filepath, "w")
    githash = HubbardFTLM.getgithash()
    for (idx, nup, ndn, tii, pii, pic) in sectors
        @mylogmsg "Sector: idx=$idx, nup=$nup, nn=$ndn, tii=$tii, pii=$pii, pic=$pic"

        if (idx=idx, hopping=t, interaction=U) in existing_results
            @mylogmsg "Sector already calculated. Skipping $idx, $nup, $ndn, $tii, $pii, $pic."
            continue
        end

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

        @mylogmsg "Hilbert space dimension: $(dimension(rhsr))"
        @mylogmsg "Generating Hamiltonian representation"
        h = represent(rhsr, hamiltonian)

        @mylogmsg "Start computing"
        hmat = Matrix(h)
        eigenvalues = eigvals!(Hermitian(hmat))
        @mylogmsg "Finished computing"

        println(
            output_file,
            JSON3.write((
                hopping=t, interaction=U, idx=idx,
                type="dense",
                eigenvalues=eigenvalues,
                timestamp=Int64(round(Dates.datetime2unix(Dates.now()) * 1000)),
                description="{\"githash\": \"$githash\"}",
            ))
        )
        flush(output_file)
    end
    close(output_file)
end


function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "lattice"
            arg_type = String
            help = "shape of the lattice in the format lattice-(?,?)x(?,?)"
            required = true
        "--hopping", "-t"
            arg_type = Float64
            default = 1.0
        "--interaction", "-U"
            arg_type = Float64
            required = true
        "--sector"
            arg_type = Int
            nargs = '*'
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

    sectors = parsed_args["sector"]

    @mylogmsg "Starting $(first(ARGS))"
    @mylogmsg "latticetype: $latticetype"
    @mylogmsg "shape: $shape"
    @mylogmsg "t: $t"
    @mylogmsg "U: $U"
    @mylogmsg "sectors: $sectors"
    @mylogmsg "BLAS: $(BLAS.get_config())"

    if isempty(sectors)
        eigen_small(latticetype, shape, t, U)
    else
        eigen_dense(latticetype, shape, t, U, sectors)
    end
    @mylogmsg "Finished writing"
end


main()

