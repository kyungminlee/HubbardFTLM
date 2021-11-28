function dsujoin(parents::AbstractVector{<:Integer}, a::Integer, b::Integer)
	a = dsufind(parents, a)
	b = dsufind(parents, b)
	if a < b
		parents[b] = a
		return true
	elseif b < a
		parents[a] = b
		return true
	else
		return false
	end
end

function dsufind(parents::AbstractVector{<:Integer}, v::Integer)
	p = parents[v]
	if p == v
		return v
	else
		p = dsufind(parents, p)
		parents[v] = p
		return p
	end
end

function findstar(lattice::Lattice, psym::PointSymmetry)
	unitcell = lattice.unitcell
	latticevectors = unitcell.latticevectors
	tsym = FiniteTranslationSymmetry(lattice);
	psym_elements = [op.matrix for op in elements(psym)]
	psymk_elements = [Int.(round.(transpose(latticevectors) * latticevectors * m * inv(latticevectors) * inv(transpose(latticevectors)), digits=8)) for m in psym_elements]
	momentum_list = tsym.fractional_momenta
	momentum_lookup = Dict(v=>i for (i, v) in enumerate(momentum_list))
	parents = collect(1:length(tsym.fractional_momenta))

	for (i, k) in enumerate(tsym.fractional_momenta)
		for m in psymk_elements
			k2 = mod.(m * k, 1)
			i2 = momentum_lookup[k2]
			dsujoin(parents, i, i2)
		end
	end
	for x in 1:length(tsym.fractional_momenta)
		dsufind(parents, x)
	end
	return parents
end
