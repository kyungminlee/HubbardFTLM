import Logging
import Dates
import Formatting

export my_metafmt
export get_rss
export parse_shape
export shape_string
export parse_lattice
export lattice_string

function get_rss()
    result = Vector{Csize_t}(undef, 1)
    ccall(:uv_resident_set_memory, Cint, (Ptr{Csize_t},), Ref(result,1))
    return result[1]
end

BN(x::Integer) = Formatting.format(Int(x), commas=true)

function my_metafmt(level, _module, group, id, file, line)
    color = Logging.default_logcolor(level)
    #prefix = (level == Logging.Warn ? "Warning" : string(level))*':'
    #prefix = string(Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS" ))* " | " * string(level) * " | " * string(id) * " | RSS: $(BN(get_rss())) | "
    date_string = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.sss")
    # prefix = "$(date_string) | $(level) | $(id) | RSS: $(BN(get_rss())) | "
    prefix = "$(date_string) | $(level) | RSS: $(BN(get_rss())) |"
    suffix = ""
    Logging.Info <= level < Logging.Warn && return color, prefix, suffix
    _module !== nothing && (suffix *= "$(_module)")
    if file !== nothing
        _module !== nothing && (suffix *= " ")
        suffix *= Base.contractuser(file)
        if line !== nothing
            suffix *= ":$(isa(line, UnitRange) ? "$(first(line))-$(last(line))" : line)"
        end
    end
    !isempty(suffix) && (suffix = "@ " * suffix)
    return color, prefix, suffix
end


global loglock = ReentrantLock()
macro mylogmsg(msg)
    return quote
        global loglock
        lock(loglock)
        @info $(esc(msg))
        flush(stdout)
        unlock(loglock)
    end
end

function parse_lattice(lattice_str::AbstractString)
    lattice_pattern = r"(\w+)-\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)x\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)"
    m = match(lattice_pattern, lattice_str)
    if isnothing(m)
        throw(ArgumentError("shape should be in format (n1a,n1b)x(n2a,n2b)"))
    end
    type = m.captures[1]
    n1a = parse(Int, m.captures[2])
    n1b = parse(Int, m.captures[3])
    n2a = parse(Int, m.captures[4])
    n2b = parse(Int, m.captures[5])
    return (type=type, shape=[n1a n2a; n1b n2b])
end

function lattice_string(type::AbstractString, shape::AbstractMatrix{<:Integer})
    n1a = shape[1,1]
    n1b = shape[2,1]
    n2a = shape[1,2]
    n2b = shape[2,2]
    return "$(strip(type))-($n1a,$n1b)x($n2a,$n2b)"
end


function parse_shape(shape_str::AbstractString)
    shape_pattern = r"\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)x\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)"
    m = match(shape_pattern, shape_str)
    if isnothing(m)
        throw(ArgumentError("shape should be in format (n1a,n1b)x(n2a,n2b)"))
    end
    n1a, n1b, n2a, n2b = [parse(Int, x) for x in m.captures]
    return [n1a n2a;
            n1b n2b]
end

function shape_string(shape::AbstractMatrix{<:Integer})
    n1a = shape[1,1]
    n1b = shape[2,1]
    n2a = shape[1,2]
    n2b = shape[2,2]
    return "($n1a,$n1b)x($n2a,$n2b)"
end
