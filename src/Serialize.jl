using Random
using DataStructures


function myserialize(rng::MersenneTwister)
    return OrderedDict(
        "seed" => rng.seed,
        "state_val" => rng.state.val,
        "vals" => rng.vals,
        "ints" => string.(rng.ints),
        "idxF" => rng.idxF,
        "idxI" => rng.idxI,
        "adv" => rng.adv,
        "adv_jump" => string(rng.adv_jump),
        "adv_vals" => rng.adv_vals,
        "adv_ints" => rng.adv_ints,
    )
end

function mydeserialize(::Type{MersenneTwister}, obj::AbstractDict)
    return MersenneTwister(
        Vector{UInt32}(obj["seed"]),
        Random.DSFMT.DSFMT_state(Vector{Int32}(obj["state_val"])),
        Vector{Float64}(obj["vals"]),
        parse.(UInt128, obj["ints"]),
        obj["idxF"],
        obj["idxI"],
        obj["adv"],
        parse(BigInt, obj["adv_jump"]),
        obj["adv_vals"],
        obj["adv_ints"],
    )
end

