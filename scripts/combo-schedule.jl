using DrWatson
@quickactivate "Hubbard"

using ArgParse
using Hubbard

using DataFrames
using Arrow

function schedule_sectors(latticetype::AbstractString, shape::AbstractMatrix{<:Integer}, lo::Integer, hi::Integer)
    @assert 0 <= lo <= hi
    lattice_str = lattice_string(latticetype, shape)

    schedule = NamedTuple{(:idx, :type), Tuple{Int, String}}[]
    df = DataFrame(Arrow.Table(datadir("sectors-$lattice_str.arrow")))
    for row in eachrow(df)
        row.dim <= 0 && continue            
        if row.dim <= lo
            type = "small"
        elseif row.dim <= hi
            type="dense"
        else
            type = "sparse"
        end
        push!(schedule, (idx=row.idx, type=type))
    end
    schedule_filepath = datadir(lattice_str, "schedule.arrow")
    schedule_directory = dirname(schedule_filepath)
    isdir(schedule_directory) || mkpath(schedule_directory)
    Arrow.write(schedule_filepath, schedule)
end

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @add_arg_table! s begin
        "lattice"
            arg_type = String
            help = "shape of the lattice in the format lattice-(?,?)x(?,?)"
            required = true
        "lo"
            arg_type = Int
            required = true
            help = "0 < dim <= lo are small"
        "hi"
            arg_type = Int
            required = true
            help = "lo < dim <= hi are dense, hi < dim are sparse"
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    latticetype, shape = parse_lattice(parsed_args["lattice"])
    lo = parsed_args["lo"]
    hi = parsed_args["hi"]
    schedule_sectors(latticetype, shape, lo, hi)
end

main()
