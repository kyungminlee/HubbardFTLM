using LinearAlgebra
using LatticeTools

#=
    .   .   .   x   .   .   .

  .   .   x   o   o   x   .   .
               \ /
    .   .   .   O - o   .   .

  .   .   .   .   .   .   .   .

=#
function make_triangular_lattice(shape::AbstractMatrix{<:Integer})
    latticevectors = [1 -0.5; 0 0.5*sqrt(3.0)];
    unitcell = make_unitcell(latticevectors, SiteType=String)
    addsite!(unitcell, "A", carte2fract(unitcell, [0.0, 0.0]))
    nnbondtypes = [ [1, 0],
                    [1, 1],
                    [0, 1] ]
    nnnbondtypes = [[ 2, 1],
                    [ 1, 2],
                    [-1, 1] ]


    nntriangletypes = [
        ([0,0], [1,0], [1,1]), # up
        ([0,0], [1,1], [0,1]), # down
    ]

    lattice = make_lattice(unitcell, shape)
    hypercube = lattice.hypercube
    supercell = lattice.supercell
    tsym = FiniteTranslationSymmetry(lattice)
    psym = little_symmetry(tsym, PointSymmetryDatabase.find2d("6mm"))
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

    nn_triangles = Tuple{
                        Vector{
                            Tuple{Tuple{Int, Int}, Vector{Int}}
                        },
                        Int}[]
    for r in lattice.bravais_coordinates
        triangle = Tuple{Tuple{Int, Int}, Vector{Int}}[]
        for i1 in 1:3
            i2 = mod(i1, 3) + 1
            r1 = nntriangletypes[1][i1]
            r2 = nntriangletypes[1][i2]
            R_row, r_row = hypercube.wrap(r .+ r1)
            R_col, r_col = hypercube.wrap(r .+ r2)
            roworb_super = ("A", r_row)
            colorb_super = ("A", r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(triangle, ((irow, icol), R_col - R_row))
        end
        push!(nn_triangles, (triangle, 1))
    end

    for r in lattice.bravais_coordinates
        triangle = Tuple{Tuple{Int, Int}, Vector{Int}}[]
        for i1 in 1:3
            i2 = mod(i1, 3) + 1
            r1 = nntriangletypes[2][i1]
            r2 = nntriangletypes[2][i2]
            R_row, r_row = hypercube.wrap(r .+ r1)
            R_col, r_col = hypercube.wrap(r .+ r2)
            roworb_super = ("A", r_row)
            colorb_super = ("A", r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(triangle, ((irow, icol), R_col - R_row))
        end
        push!(nn_triangles, (triangle, -1))
    end
            
    return (unitcell=unitcell,
            lattice=lattice,
            space_symmetry_embedding=ssymbed,
            nearest_neighbor_bonds=nnbonds,
            next_nearest_neighbor_bonds=nnnbonds,
            nearest_neighbor_triangles=nn_triangles)
end






