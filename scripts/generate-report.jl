using DrWatson
@quickactivate "HubbardFTLM"

using HubbardFTLM
using ArgParse

function generate_report(
    latticetype::AbstractString,
    shape::AbstractMatrix{<:Integer},
)
    if latticetype == "triangular"
        mylattice = make_triangular_lattice(shape)
    elseif latticetype == "square"
        mylattice = make_square_lattice(shape)
    else
        @error "Unsupported lattice type $latticetype"
        exit(1)
    end
    mylattice = make_triangular_lattice(shape)
    lattice_str = shape_string(latticetype, shape)
    output_filepath = datadir("lattice-$(lattice_str).json")
    psymbed = mylattice.space_symmetry_embedding.rest

    write_lattice_yaml(mylattice.lattice, psymbed.hermann_mauguin, output_filepath)
end

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
    generate_report(latticetype, shape)
end

main()
