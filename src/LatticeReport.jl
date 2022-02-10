using LatticeTools
using JSON3
using DataStructures

tostring(f::Rational{<:Integer}) = "$(f.num)/$(f.den)"
tostring(f::Vector{<:Rational{<:Integer}}) = "($(join(tostring.(f), " ")))"

function write_lattice_yaml(lattice::Lattice, point_symmetry_name::AbstractString, output_filename::AbstractString="out.yaml")
    lattice_shape = lattice.hypercube.shape_matrix
    @assert size(lattice_shape) == (2,2)
    tsym = FiniteTranslationSymmetry(lattice)
    psym = little_symmetry(tsym, PointSymmetryDatabase.find2d(point_symmetry_name))

    tsymbed = embed(lattice, tsym)
    psymbed = embed(lattice, psym)

    pair_groups = []
    let
        visited = falses(numsites(lattice.supercell), numsites(lattice.supercell))
        for (i1, (sitename1, sitefc1)) in enumerate(lattice.supercell.sites),
            (i2, (sitename2, sitefc2)) in enumerate(lattice.supercell.sites)
            visited[i1, i2] && continue
            pair_group = Set()
            for t in tsymbed.elements, p in psymbed.elements
                j1 = p(t(i1))
                j2 = p(t(i2))
                push!(pair_group, [j1, j2])
                visited[j1, j2] = true
            end
            push!(pair_groups, sort(collect(pair_group)))
        end
    end

    y_unitcell_sites = []
    for (sitename, sitefc) in lattice.unitcell.sites
        sitecc = fract2carte(lattice.unitcell, sitefc)
        push!(y_unitcell_sites, OrderedDict(
                "name" => string(sitename),
                "fractional_coordinates" => OrderedDict(
                    "whole" => sitefc.whole,
                    "fraction" => sitefc.fraction,
                ),
                "cartesian_coordinates" => sitecc
                ))
    end
    y_unitcell = OrderedDict(
        "lattice_vectors" => collect(eachcol(lattice.unitcell.latticevectors)),
        "reciprocallatticevectors" => collect(eachcol(lattice.unitcell.reciprocallatticevectors)),
        "sites" => y_unitcell_sites,
    #     "orbitals" => lattice.unitcell.orbitals,
    )

    y_supercell_sites = []
    for (sitename, sitefc) in lattice.supercell.sites
        sitecc = fract2carte(lattice.supercell, sitefc)
        push!(y_supercell_sites, OrderedDict(
                "name" => join(string.(sitename), "_"),
                "fractional_coordinates" => OrderedDict(
                    "whole" => sitefc.whole,
                    "fraction" => sitefc.fraction,
                ),
                "cartesian_coordinates" => sitecc
                ))
    end
    y_supercell = OrderedDict(
        "lattice_vectors" => collect(eachcol(lattice.supercell.latticevectors)),
        "reciprocallatticevectors" => collect(eachcol(lattice.supercell.reciprocallatticevectors)),
        "sites" => y_supercell_sites,
    #     "orbitals" => lattice.unitcell.orbitals,
    )

    y_site_pairs = pair_groups

    y_momentums = []
    let
        K = lattice.unitcell.reducedreciprocallatticevectors
        k = K * hcat(tsym.fractional_momenta...)
        for (kcol, kfcol) in zip(eachcol(k), tsym.fractional_momenta)
            push!(y_momentums, OrderedDict(
                    "cartesian_coordinates" => kcol,
                    "fractional_coordinates" => 
    #                 [
    #                     OrderedDict("numerator" => x.num, "denominator" => x.den) for x in kfcol
    #                 ],
                    [
                        "$(x.num)/$(x.den)" for x in kfcol
                    ],
                    )
                )
        end
    end    

    y_output = OrderedDict(
        "shape" => collect(eachcol(lattice_shape)),
        "unitcell" => y_unitcell,
        "supercell" => y_supercell,
        "bravais_coordinates" => lattice.bravais_coordinates,
        "momentums" => y_momentums,
        "point_group_symmetry" => psym.hermann_mauguin,
        "equivalent_site_pairs" => y_site_pairs,
    )
    open(output_filename, "w") do f
        write(f, JSON3.json(y_output))
    end
end
