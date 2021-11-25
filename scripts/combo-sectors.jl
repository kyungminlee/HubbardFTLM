using DrWatson
@quickactivate "Hubbard"

using Hubbard

using LinearAlgebra
using Random
using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle

using Formatting
using Logging
using ArgParse

using DataFrames
using ProgressMeter
using Arrow

function compute_sectors(latticetype::AbstractString, shape::AbstractMatrix{<:Integer})
    BR = UInt
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

    idx = 0
    total_dimension = 0
    sectors = let
        tmp = Dict(qn => [] for qn in quantum_number_sectors(hs))
        nqns = length(quantum_number_sectors(hs))
        lck = Threads.SpinLock()
        prog = Progress(4^n_sites)
        processed_qns = 0
        Threads.@threads for qn in quantum_number_sectors(hs)
            nup, ndn = qn
            hsr = represent(hs, [qn], BR, DictIndexedVector)
            local_total_dimension = 0
            for ((tii, pii, pic), ssa) in ssa_collection
                rhsr = symmetry_reduce(hsr, [s for (s, p) in ssa], [p for (s, p) in ssa])
                local_total_dimension += dimension(rhsr)
                push!(tmp[qn], (nup=nup, ndn=ndn, tii=tii, pii=pii, pic=pic, dim=dimension(rhsr)))
            end
            lock(lck)
            total_dimension += local_total_dimension
            processed_qns += 1
            ProgressMeter.update!(prog, total_dimension; showvalues=[(:D, "$(total_dimension)/$(4^n_sites)"), (:qn, "$(processed_qns)/$(nqns)")])
            unlock(lck)
        end # for qn

        T = typeof((idx=1,nup=1,ndn=1,tii=1,pii=1,pic=1,dim=1))
        sectors = T[]
        sizehint!(sectors, length(ssa_collection) * length(quantum_number_sectors(hs)))

        idx = 0
        for qn in sort(collect(keys(tmp))), x in tmp[qn]
            idx += 1
            nup, ndn = qn
            push!(sectors, (idx=idx, nup=x.nup, ndn=x.ndn, tii=x.tii, pii=x.pii, pic=x.pic, dim=x.dim))
        end
        sectors
    end
    @info "total dimension: $(sum(x.dim for x in sectors))"
    @info "Writing to database"

    sectors_df = DataFrame(sectors)
    isdir(datadir()) || mkpath(datadir())
    Arrow.write(datadir("sectors-$lattice_str.arrow"), sectors_df)

    # h5f = h5open(datadir("sectors-$lattice_str.hdf5"), "w")
    # h5g_sectors = create_group(h5f, "sector")
    # h5g_sectors["idx"] = sectors_df.idx
    # h5g_sectors["nup"] = sectors_df.nup
    # h5g_sectors["ndn"] = sectors_df.ndn
    # h5g_sectors["tii"] = sectors_df.tii
    # h5g_sectors["pii"] = sectors_df.pii
    # h5g_sectors["pic"] = sectors_df.pic
    # h5g_sectors["dim"] = sectors_df.dim
    # attributes(h5g_sectors)["n_sites"] = n_sites
    # attributes(h5g_sectors)["lattice"] = lattice_str
    # attributes(h5g_sectors)["latticetype"] = latticetype
    # attributes(h5g_sectors)["shape"] = shape
    # close(h5f)
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

