# HubbardFTLM

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> HubbardFTLM

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Finite Temperature Lanczos Pipeline

The scripts in the `scripts` folder support finite temperature Lanczos (+ exact) calculations for
the Hubbard model on both triangular and square lattices.
The whole computation is split into the following steps.
1. First, the symmetry sectors are identified. This is in terms of (1) number of up electrons, (2) number of down electrons,
(3) translation symmetry, and (4) little point symmetry compatible with the translation symmetry.
2. The dimension of each sector is calculated.
3. For small enough sectors, all the eigenvalues of each sector are computed by exact diagonalization.
For large sectors, the eigenvalues and the squares of the first components of the eigenvalues are computed.

The pipeline for the whole process is as follows:
```bash
julia init-empty.jl lattice_str
julia init-sectors.jl lattice_str [--sector s1 s2 ...]
julia merge-sectors.jl lattice_str
julia schedule.jl lattice_str lo hi
julia eigen-dense.jl lattice_str -t t -U U [--sector s1 s2 ...]
julia eigen-sparse.jl lattice_str -t t -U U [--sector s1 s2 ...] --krylovdim krylovdim ...
```
- `init-empty.jl` identifies all possible sectors, regardless of their dimensions (which can be zero). It also finds the equivalence between the sectors. The results are saved to `sectors-lattice_str.arrow`.
- `init-sectors.jl` computes the dimensions of the sectors, and saves the results to `sectors-....jsonl`.
- `merge-sectors.jl` merges the `jsonl` files into the `sectors-lattice_str.arrow`.
- `schedule.jl` partitions the sections into `small`, `dense`, and `sparse` depending on their dimensions.
- `eigen-dense.jl` computes the eigenvalues of the `small` and `dense` sectors.
- `eigen-sparse.jl` computes the Ritz values and the squares of the first elements of the corresponding Ritz vectors.
