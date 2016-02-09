module MultizoneTransport

import BrownianSources

abstract Zone

type UnsaturatedZone <: Zone
	v::Float64
	dispersivity::Float64
	drippoint::Float64
end

type SaturatedZone <: Zone
	x0::Array{Float64, 1}
	sigma0::Array{Float64, 1}
	v::Array{Float64, 1}
	dispersivity::Array{Float64, 1}
	drippoint::Float64#the x coordinate defining the plane/line/point where this saturated zone starts dripping down into the unsaturated zone below
	anasolfunc::Function
end

type Multizone
	uzs::Array{UnsaturatedZone, 1}#unsaturated zones
	szs::Array{SaturatedZone, 1}#saturated zones
	maxtime::Float64
	releasemass::Float64
	bss::Array{BrownianSources.BrownianSource, 1}
	sourcestrengths::Array{Function, 1}
end

function Multizone(uzs, szs, releasemass, maxtime, numsamples)
	@assert length(uzs) == length(szs)
	bss = Array(BrownianSources.BrownianSource, length(szs))
	sourcestrengths = Array(Function, length(szs))
	for i = 1:length(szs)
		bss[i] = getbrowniansource([uzs; szs[1:i - 1]], maxtime, numsamples)
		sourcestrengths[i] = t->bss[i](t)
	end
	return Multizone(uzs, szs, maxtime, releasemass, bss, sourcestrengths)
end

function getbrowniansource(zones, maxtime, numsamples)
	velocities = Array(Float64, length(zones))
	dispersivities = Array(Float64, length(zones))
	thresholds = Array(Float64, length(zones))
	i = 1
	for zone in zones
		velocities[i] = zones[i].v[1]
		dispersivities[i] = zones[i].dispersivity[1]
		thresholds[i] = zones[i].drippoint
		i += 1
	end
	return BrownianSources.BrownianSource(velocities, dispersivities, thresholds, maxtime, numsamples)
end

function getanasolargs(sz::SaturatedZone, ::Type{Val{1}})
	H = .5
	#TODO properly incorporate xb in here
	xb = NaN
	#TODO properly incorporate lambda in here
	lambda = 0.
	t0 = 0.
	return (sz.x0[1], sz.sigma0[1], sz.v[1], sqrt(2 * sz.v[1] * sz.dispersivity[1]), H, xb, lambda)
end

f(x, t, anasolargs...) = 1.

function getconcentration(mz::Multizone, x::Vector, t::Real, szindex::Int)
	#return releasemass * szs[i].anasolfunc(x, t, getanasolargs(szs[i], Val{length(x)}))
	return mz.releasemass * mz.szs[szindex].anasolfunc(x, t, getanasolargs(mz.szs[szindex], Val{length(x)})..., 0., mz.maxtime, mz.sourcestrengths[szindex])
end

end
