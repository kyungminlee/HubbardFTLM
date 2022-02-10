using Random
using KrylovKit
using QuantumHamiltonian
using QuantumHamiltonianParticle
using FiniteTemperatureLanczos
using ProgressMeter

export make_hamiltonian
export make_hamiltonian_ladder
export compute_thermodynamics_dense
export compute_thermodynamics_sparse

function _prodsum(a::AbstractVector{A}, b::AbstractVector{B}) where {A, B}
    Z = promote_type(A, B)
    out = zero(Z)
    for (x, y) in zip(a, b)
        out += x * y
    end
    return out
end

function _prodsum(a::AbstractVector{A}, b::AbstractVector{B}, c::AbstractVector{C}) where {A, B, C}
    Z = promote_type(A, B, C)
    out = zero(Z)
    for (x, y, z) in zip(a, b, c)
        out += x * y * z
    end
    return out
end

function compute_thermodynamics_dense(
    H::AbstractMatrix{<:Number},
    temperatures::AbstractVector{<:AbstractFloat},
    observables::AbstractVector{<:AbstractObservable},
)
    D = size(H, 1)
    @mylogmsg "generating matrix of size $D (dense)"
    m = Matrix(H)
    # @mylogmsg "deviation from Hermiticity: $(maximum(abs2.(m - m')))"
    @mylogmsg "diagonalizing matrix of size $D (dense)"
    eigenvalues, eigenvectors = eigen!(Hermitian(m))

    @mylogmsg "diagonalized matrix of size $D (dense)"
    z = Vector{Float64}(undef, length(temperatures))
    E = Vector{Float64}(undef, length(temperatures))
    E2 = Vector{Float64}(undef, length(temperatures))
    obs_list = [Vector{eltype(obs)}(undef, length(temperatures)) for obs in observables]
    
    base_energy = minimum(eigenvalues)
    shifted_eigenvalues = eigenvalues .- base_energy
    #eigenvalues .-= base_energy

    boltzmann_T = [exp.(-shifted_eigenvalues ./ T) for T in temperatures]
    for (iT, T) in enumerate(temperatures)
        boltzmann = boltzmann_T[iT]
        z[iT] = sum(boltzmann)
        E[iT] = _prodsum(shifted_eigenvalues, boltzmann)
        E2[iT] = _prodsum(shifted_eigenvalues, shifted_eigenvalues, boltzmann)
    end

    obsU, UOUdiag = let _,
            U = eigenvectors,
            T = promote_type((eltype(obs.observable) for obs in observables)..., eltype(U))
        (
            Matrix{T}(undef, (size(first(observables).observable, 1), size(U, 2))),
            Vector{T}(undef, length(shifted_eigenvalues))
        )
    end
    for (iobs, obs) in enumerate(observables)
        U = eigenvectors
        @mylogmsg "Computing O * U matrix for observable #$(iobs)"
        mul!(obsU, obs.observable, U)
        @mylogmsg "Computing observable #$(iobs)"
        Threads.@threads for i in eachindex(shifted_eigenvalues)
            UOUdiag[i] = dot(U[:, i], obsU[:, i])
        end
        for iT in eachindex(temperatures)
            boltzmann = boltzmann_T[iT]
            obs_list[iobs][iT] = _prodsum(UOUdiag, boltzmann)
        end
    end

    @mylogmsg "Finished computing observables"
    return (matrix_type=:dense, base_energy=base_energy, z=z, E=E, E2=E2, obs=obs_list, eigenvalues=eigenvalues)
end


function compute_thermodynamics_sparse(
    H::AbstractMatrix{<:Number},
    temperatures::AbstractVector{<:AbstractFloat},
    observables::AbstractVector{<:AbstractObservable};
    rng::AbstractRNG=Random.GLOBAL_RNG,
    krylovdim::Integer=200,
    nsamples::Integer=1000,
    base_energy::Union{<:Real, <:Nothing}=nothing,
)
    # @assert isempty(observables)
    D = size(H, 1)
    krylovdim = min(D, krylovdim)
    @mylogmsg "matrix: $D (sparse / $krylovdim)"

    zs  = [zeros(Float64, nsamples) for T in temperatures]
    Es  = [zeros(Float64, nsamples) for T in temperatures]
    E2s = [zeros(Float64, nsamples) for T in temperatures]
    obss_list = [[zeros(eltype(obs), nsamples) for T in temperatures] for obs in observables]

    #@showprogress for r in 1:nsamples
    work = Vector{ComplexF64}(undef, D+2*krylovdim^2)
    
    pre_z = Vector{Float64}(undef, krylovdim)
    pre_E = Vector{Float64}(undef, krylovdim)
    pre_E2 = Vector{Float64}(undef, krylovdim)
    pre_obss = [Matrix{ComplexF64}(undef, krylovdim, krylovdim) for obs in observables]

    for r in 1:nsamples
        @mylogmsg "r=$r"
        v = randn(rng, ComplexF64, D)
        normalize!(v)
        iterator = LanczosIterator(x -> H*x, v)
        factorization = initialize(iterator)
        sizehint!(factorization, krylovdim)
        @mylogmsg "Done sizehint!"
        for d in 1:krylovdim-1
            expand!(iterator, factorization)
            @mylogmsg "Expanded factorization to $d"
            # GC.gc(false)
        end
        @mylogmsg "Garbage collecting after Krylov"
        GC.gc()

        @mylogmsg "Premeasurement"
        pm = premeasure(factorization)
        if isnothing(base_energy)
            base_energy = minimum(pm.eigen.values)
        end
        eigenvalues = copy(pm.eigen.values)
        pm.eigen.values .-= base_energy
        #pre_z  = premeasure(pm, 0)
        premeasure!(pre_z, pm)
        @mylogmsg "Finished premeasure z"
        #pre_E  = premeasure(pm, 1)
        premeasure!(pre_E, pm, 1)
        @mylogmsg "Finished premeasure E"
        premeasure!(pre_E2, pm, 2)
        #pre_E2 = premeasure(pm, 2)

        @mylogmsg "Finished premeasure E²"

        for (iobs, obs) in enumerate(observables)
            premeasure!(pre_obss[iobs], pm, obs, work)
        end
        #pre_obss = [premeasure(pm, obs, work) for obs in observables]
        @mylogmsg "Garbage collecting after premeasure obs"
        GC.gc()
        @mylogmsg "Finished garbage collecting"

        @assert base_energy isa Real

        for (iT, T) in enumerate(temperatures)
            zs[iT][r]  = measure(pre_z,  pm.eigen.values, T)
            Es[iT][r]  = measure(pre_E,  pm.eigen.values, T)
            E2s[iT][r] = measure(pre_E2, pm.eigen.values, T)
            for (iobs, pobs) in enumerate(pre_obss)
                obss_list[iobs][iT][r] = measure(pobs, pm.eigen.values, T)
            end
        end
        @mylogmsg "Garbage collecting after measure"
        GC.gc()
    end # for r
    @mylogmsg "collected $nsamples samples"
    return (matrix_type=:sparse, base_energy=base_energy, z=zs, E=Es, E2=E2s, obs=obss_list, eigenvalues=eigenvalues)
end


function make_hamiltonian(triangular, t::Real, U::Real)
    n_sites = numsites(triangular.lattice.supercell)

    ps, c, cdag = electron_system()

    em = ParticleState(ps, "em", [0, 0], (0, 0))
    up = ParticleState(ps, "up", [1, 0], (1, 0))
    dn = ParticleState(ps, "dn", [0, 1], (0, 1))
    ud = ParticleState(ps, "ud", [1, 1], (1, 1))
    site = ParticleSite([em, up, dn, ud])

    hs = ParticleHilbertSpace([site for i in 1:n_sites])

    # Construct operators (Hamiltonian and measurement)
    hopping = ParticleLadderNull(ps)
    hopping = sum(
        cdag(j, σ) * c(i, σ) + cdag(i, σ) * c(j, σ)
            for ((i, j), R) in triangular.nearest_neighbor_bonds
            for σ in [:up, :dn]
    ) * (-t)
    interaction = sum(
        cdag(i, :up) * c(i, :up) * cdag(i, :dn) * c(i, :dn)
            for i in 1:n_sites
    ) * U
    # hamiltonian = embed(hs, hopping + interaction)
    hamiltonian = make_projector_operator(hopping + interaction)
    return hs, hamiltonian
end


function make_hamiltonian_ladder(triangular, t::Real, U::Real)
    n_sites = numsites(triangular.lattice.supercell)

    ps, c, cdag = electron_system()

    em = ParticleState(ps, "em", [0, 0], (0, 0))
    up = ParticleState(ps, "up", [1, 0], (1, 0))
    dn = ParticleState(ps, "dn", [0, 1], (0, 1))
    ud = ParticleState(ps, "ud", [1, 1], (1, 1))
    site = ParticleSite([em, up, dn, ud])

    hs = ParticleHilbertSpace([site for i in 1:n_sites])

    # Construct operators (Hamiltonian and measurement)
    hopping = ParticleLadderNull(ps)
    hopping = sum(
        cdag(j, σ) * c(i, σ) + cdag(i, σ) * c(j, σ)
            for ((i, j), R) in triangular.nearest_neighbor_bonds
            for σ in [:up, :dn]
    ) * (-t)
    interaction = sum(
        cdag(i, :up) * c(i, :up) * cdag(i, :dn) * c(i, :dn)
            for i in 1:n_sites
    ) * U
    hamiltonian = make_projector_operator(hopping + interaction)
    return hs, hamiltonian, c, cdag
end
