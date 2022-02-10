using DrWatson
@quickactivate "HubbardFTLM"

using HubbardFTLM

using LinearAlgebra
using Random
using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle

using Formatting
using Logging
using ArgParse
using UUIDs

using JSON3
using Arrow

function compute_sectors(latticetype::AbstractString, shape::AbstractMatrix{<:Integer}, indices::AbstractVector{<:Integer})
    BR = UInt
    lattice_str = lattice_string(latticetype, shape)
    sectors_table = open(datadir("sectors-$lattice_str.arrow"), "r") do io
        Arrow.Table(io)
    end
    if isempty(indices)
        indices = sectors_table.idx[sectors_table.idx .== sectors_table.root_idx]
    end

    # Set up lattice
    if latticetype == "triangular"
        mylattice = make_triangular_lattice(shape)
    elseif latticetype == "square"
        mylattice = make_square_lattice(shape)
    else
        @error "Unsupported lattice type $latticetype"
        exit(1)
    end
    n_sites = numsites(mylattice.lattice.supercell)

    # lattice = mylattice.lattice
    ssymbed = mylattice.space_symmetry_embedding
    tsymbed = ssymbed.normal
    psymbed = ssymbed.rest

    # Set up particle sector and Hilbert space
    ps, _, _ = electron_system()

    em = ParticleState(ps, "em", [0, 0], (0, 0))
    up = ParticleState(ps, "up", [1, 0], (1, 0))
    dn = ParticleState(ps, "dn", [0, 1], (0, 1))
    ud = ParticleState(ps, "ud", [1, 1], (1, 1))
    site = ParticleSite([em, up, dn, ud])

    hs = ParticleHilbertSpace([site for i in 1:n_sites])

    ssa_collection = Dict()
    for tii in 1:irrepcount(tsymbed)
        tsa = collect(get_irrep_iterator(IrrepComponent(tsymbed, tii)))
        psymbed_little = little_symmetry(tsymbed, tii, psymbed)
        for pii in 1:num_irreps(psymbed_little), pic in 1:irrep_dimension(psymbed_little, pii)
            psa = get_irrep_iterator(IrrepComponent(psymbed_little, pii, pic))
            ssa = collect(make_product_irrep(psa, tsa))
            ssa_collection[(tii, pii, pic)] = ssa
        end
    end

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

    jsonl_filepath = datadir("sectors-$lattice_str-$(uuid5(uuid1(), gethostname())).jsonl")
    jsonl_file = open(jsonl_filepath, "w")
    for idx_ in indices
        idx = sectors_table.idx[idx_]
        nup = sectors_table.nup[idx_]
        ndn = sectors_table.ndn[idx_]
        tii = sectors_table.tii[idx_]
        pii = sectors_table.pii[idx_]
        pic = sectors_table.pic[idx_]
        root_idx = sectors_table.root_idx[idx_]
        @assert idx == idx_

        @mylogmsg "Sector: idx=$idx, nup=$nup, nn=$ndn, tii=$tii, pii=$pii, pic=$pic"

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

        println(jsonl_file, JSON3.write((idx=idx, nup=nup, ndn=ndn, tii=tii, pii=pii, pic=pic, dim=D, root_idx=root_idx)))
        flush(jsonl_file)
    end
    close(jsonl_file)
end # function compute


function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "lattice"
            arg_type = String
            help = "shape of the lattice in the format lattice-(?,?)x(?,?)"
            required = true
        "indices"
            arg_type = Int
            nargs = '*'
            help = "if not provided, compute all root sectors"
    end
    parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    latticetype, shape = parse_lattice(parsed_args["lattice"])
    indices = parsed_args["indices"]
    compute_sectors(latticetype, shape, indices)
end


main()

