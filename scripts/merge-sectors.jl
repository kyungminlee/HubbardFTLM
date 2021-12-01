using DrWatson
@quickactivate "HubbardFTLM"

using HubbardFTLM

using ArgParse
using Glob

using DataFrames
using Arrow
using JSON3
using Tables

function merge_sectors(latticetype::AbstractString, shape::AbstractMatrix{<:Integer})
    lattice_str = lattice_string(latticetype, shape)
    sectors_filepath = datadir("sectors-$lattice_str.arrow")
    sectors_file = open(sectors_filepath, "r+")
    sectors_table = Arrow.Table(sectors_file)
    @info "Opened sector file $sectors_filepath"

    @assert all(row.idx == irow for (irow, row) in enumerate(Tables.rows(sectors_table)))
    jsonl_filepaths = Glob.glob("sectors-$lattice_str-*.jsonl", datadir())
    @info "Reading $(length(jsonl_filepaths)) jsonl files"
    for jsonl_filepath in jsonl_filepaths
        jsonl = open(jsonl_filepath, "r") do io
            JSON3.read(io, jsonlines=true)
        end
        for row in jsonl
            idx, dim = row["idx"], row["dim"]
            sectors_table.dim[idx] = dim
        end
    end
    for (irow, row) in enumerate(Tables.rows(sectors_table))
        sectors_table.dim[irow] = sectors_table.dim[row.root_idx]
    end
    missing_indices = sectors_table.idx[sectors_table.dim .< 0]
    if !isempty(missing_indices)
        @info "Missing sectors"
        @info missing_indices
    end
    close(sectors_file)
    for jsonl_filepath in jsonl_filepaths
        mv(jsonl_filepath, "$(jsonl_filepath).processed"; force=true)        
    end
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
    merge_sectors(latticetype, shape)
end


main()

