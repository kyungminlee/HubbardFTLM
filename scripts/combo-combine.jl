using DrWatson
@quickactivate "Hubbard"

using Hubbard
using DataFrames
using CSV
using ProgressMeter

using Parquet
using ArgParse
using Glob

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
    # latticetype, shape = parse_lattice(parsed_args["lattice"])
    lattice_str = parsed_args["lattice"]
    dense_filepath = datadir("weiss", lattice_str, "eigen-dense-results.parquet")
    filenames = Glob.glob("temp-eigen-*.csv", datadir("weiss", lattice_str))
    
    ResultType = NamedTuple{(:idx, :hopping, :interaction, :eigenindex, :eigenvalue),
                         Tuple{Int, Float64, Float64, Int, Float64}}
    if !isfile(dense_filepath)
        results_table = DataFrame(ResultType[])
    else
        table = read_parquet(dense_filepath)
        results_table = DataFrame(table)
        close(table)
    end
    @info "Read $dense_filepath with $(nrow(results_table)) rows."

    for filename in filenames
        df = CSV.read(filename, DataFrame)
        append!(results_table, df; cols=:setequal)
        @info "Appended $filename with $(nrow(df)) rows."
    end
    @info "Total $(nrow(results_table)) rows before removing duplicates"
    results_table = results_table[.! nonunique(results_table, [:idx, :hopping, :interaction, :eigenindex]), :]
    @info "Total $(nrow(results_table)) rows after removing duplicates"

    write_parquet("$dense_filepath.new", results_table)
    @info "Saved results to $dense_filepath.new"
    mv("$dense_filepath.new", dense_filepath; force=true)
    @info "Renamed the file to $dense_filepath"
    for filename in filenames
        mv(filename, "$filename.processed")
    end
end

main()