using DrWatson
@quickactivate "Hubbard"

using Hubbard

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

using DataFrames
using CSV

using Parquet

using UUIDs

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
    ;
    commitevery::Integer=1
)
    lattice_str = lattice_string(latticetype, shape)

    schedule_filepath = datadir("weiss", lattice_str, "schedule.csv")
    schedule_table = CSV.read(schedule_filepath, DataFrame)
    filter!(row->row.type == "small", schedule_table)
    sectors = Int[]
    for row in eachrow(schedule_table)
        push!(sectors, row.idx)
    end
    return eigen_dense(latticetype, shape, t, U, sectors; commitevery=commitevery)
end


function eigen_dense(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
    sector_indices::AbstractVector{<:Integer},
    ;
    force::Bool=false,
    commitevery::Integer=1,
)
    lattice_str = lattice_string(latticetype, shape)
    sectors_filepath = datadir("sectors-$lattice_str.csv")

    select_table = DataFrame(idx=sector_indices)
    sectors_table = CSV.read(sectors_filepath, DataFrame)

    table = innerjoin(select_table, sectors_table, on=:idx)
    sectors = SectorType[]
    for row in eachrow(table)
        push!(sectors, (idx=row.idx, nup=row.nup, ndn=row.ndn, tii=row.tii, pii=row.pii, pic=row.pic))
    end
    return eigen_dense(latticetype, shape, t, U, sectors; commitevery=commitevery, force=force)
end


function eigen_dense(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
    t::Real,
    U::Real,
    sectors::AbstractVector{SectorType},
    ;
    force::Bool=false,
    commitevery::Integer=1,
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

    @mylogmsg "Opening DB file"

    lattice_str = lattice_string(latticetype, shape)
    dense_filepath = datadir("weiss", lattice_str, "eigen-dense-results.parquet")
    temp_filepath = datadir("weiss", lattice_str, "temp-eigen-dense-$(uuid5(uuid1(), gethostname())).csv")

    ResultType = NamedTuple{(:idx, :hopping, :interaction, :eigenindex, :eigenvalue),
                         Tuple{Int, Float64, Float64, Int, Float64}}
    if !isfile(dense_filepath)
        results_table = DataFrame(ResultType[])
    else
        results_table = DataFrame(read_parquet(dense_filepath))
    end

    output_file = open(temp_filepath, "w")
    println(output_file, "idx,hopping,interaction,eigenindex,eigenvalue")

    for (idx, nup, ndn, tii, pii, pic) in sectors
        @mylogmsg "Sector: idx=$idx, nup=$nup, nn=$ndn, tii=$tii, pii=$pii, pic=$pic"

        if any(row.idx == idx && row.hopping == t && row.interaction == U for row in eachrow(results_table))
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

        for (eigenindex, eigenvalue) in enumerate(eigenvalues)
            println(output_file, "$idx,$t,$U,$eigenindex,$eigenvalue")
        end
        flush(output_file)
    end
    close(output_file)
end


function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "shape"
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
        "--force", "-f"
            action = :store_true
        "--debug", "-d"
            action = :store_true
        "--commitevery"
            arg_type = Int
            default = 1
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

    latticetype, shape = parse_lattice(parsed_args["shape"])
    t = parsed_args["hopping"]
    U = parsed_args["interaction"]

    sectors = parsed_args["sector"]
    commitevery = parsed_args["commitevery"]
    force = parsed_args["force"]

    @mylogmsg "Starting $(first(ARGS))"
    @mylogmsg "shape: $shape"
    @mylogmsg "t: $t"
    @mylogmsg "U: $U"
    @mylogmsg "sectors: $sectors"
    @mylogmsg "BLAS: $(BLAS.get_config())"
    @mylogmsg "commitevery: $commitevery"
    @mylogmsg "force: $force"

    if isempty(sectors)
        eigen_small(latticetype, shape, t, U; commitevery=commitevery)
    else
        eigen_dense(latticetype, shape, t, U, sectors; commitevery=commitevery, force=force)
    end
    @mylogmsg "Finished writing"
end


main()

