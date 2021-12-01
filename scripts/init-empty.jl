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

using DataFrames
using Arrow

function compute_sectors(latticetype::AbstractString, shape::AbstractMatrix{<:Integer})
    # BR = UInt
    lattice_str = lattice_string(latticetype, shape)

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
    
    ssym = mylattice.space_symmetry
    tsym = ssym.normal
    psym = ssym.rest

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


    ssa_collection = []
    for tii in 1:irrepcount(tsymbed)
        tsa = collect(get_irrep_iterator(IrrepComponent(tsymbed, tii)))
        psymbed_little = little_symmetry(tsymbed, tii, psymbed)
        for pii in 1:num_irreps(psymbed_little), pic in 1:irrep_dimension(psymbed_little, pii)
            psa = get_irrep_iterator(IrrepComponent(psymbed_little, pii, pic))
            ssa = collect(make_product_irrep(psa, tsa))
            push!(ssa_collection, ((tii, pii, pic), ssa))
        end
    end

    sectors = let
        T = typeof((idx=1,nup=1,ndn=1,tii=1,pii=1,pic=1,dim=1,root_idx=1))
        sectors = T[]
        sizehint!(sectors, length(ssa_collection) * length(quantum_number_sectors(hs)))
        idx = 0
        for qn in quantum_number_sectors(hs)
            nup, ndn = qn
            for ((tii, pii, pic), ssa) in ssa_collection
                idx += 1
                push!(sectors, (idx=idx, nup=nup, ndn=ndn, tii=tii, pii=pii, pic=pic, dim=-1, root_idx=-1))
            end
        end # for qn
        sectors
    end
    sector_lookup = Dict((nup, ndn, tii, pii, pic) => idx for (idx, nup, ndn, tii, pii, pic, dim, root_idx) in sectors)

    star = findstar(mylattice.lattice, psym)

    for (irow, row) in enumerate(sectors)
        (idx, nup, ndn, tii, pii, pic, dim, root_idx) = row
        @assert(irow == idx)
        if nup >= ndn
            root_nup = nup
            root_ndn = ndn
        else
            root_nup = ndn
            root_ndn = nup
        end
        root_tii = star[tii]
        root_pii = pii
        root_pic = 1
        root_idx = sector_lookup[(root_nup, root_ndn, root_tii, root_pii, root_pic)]
        sectors[irow] = (idx=idx, nup=nup, ndn=ndn, tii=tii, pii=pii, pic=pic, dim=dim, root_idx=root_idx)
    end

    @info "Writing to database"

    isdir(datadir()) || mkpath(datadir())

    Arrow.write(datadir("sectors-$lattice_str.arrow"), sectors)
end # function compute


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
    latticetype, shape = parse_lattice(parsed_args["lattice"])
    compute_sectors(latticetype, shape)
end


main()

