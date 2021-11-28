using LinearAlgebra
using LatticeTools

function make_square_lattice(shape::AbstractMatrix{<:Integer})
    latticevectors = [1.0 0.0; 0.0 1.0]
    unitcell = make_unitcell(latticevectors, SiteType=String)
    addsite!(unitcell, "A", carte2fract(unitcell, [0.0, 0.0]))
    nnbondtypes = [ [1, 0],
                    [0, 1]]
    nnnbondtypes = [[ 1, 1],
                    [ 1,-1]]

    lattice = make_lattice(unitcell, shape)
    hypercube = lattice.hypercube
    supercell = lattice.supercell
    tsym = FiniteTranslationSymmetry(lattice)
    psym = little_symmetry(tsym, PointSymmetryDatabase.find2d("4mm"))
    ssym = SymmorphicSymmetry(tsym, psym)
    tsymbed = embed(lattice, tsym)
    psymbed = embed(lattice, psym)
    ssymbed = SymmorphicSymmetryEmbedding(tsymbed, psymbed)

    nnbonds = []
    nnnbonds = []

    for r_row in lattice.bravais_coordinates
        for colvec in nnbondtypes
            R_col, r_col = hypercube.wrap(r_row .+ colvec)
            roworb_super = ("A", r_row)
            colorb_super = ("A", r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(nnbonds, ((irow, icol), R_col))
        end
        for colvec in nnnbondtypes
            R_col, r_col = hypercube.wrap(r_row .+ colvec)
            roworb_super = ("A", r_row)
            colorb_super = ("A", r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(nnnbonds, ((irow, icol), R_col))
        end
    end

    return (unitcell=unitcell,
            lattice=lattice,
            space_symmetry=ssym,
            space_symmetry_embedding=ssymbed,
            nearest_neighbor_bonds=nnbonds,
            next_nearest_neighbor_bonds=nnnbonds
           )
end






