function getparamvec(params, prefix)
	i = Array(Int, 0)
	p = Array(Float64, 0)
	for k in keys(params)
		if startwith(k, prefix)
			push!(i, parse(Int, k[length(prefix) + 1:end]))
			push!(p, params[k])
		end
	end
	pvec = Array(Float64, length(i))
	for j = 1:length(i)
		pvec[j[i]] = p[j]
	end
	return pvec
end

function assemblemzs(params, nummzs, numsamples)
	mzs = Array(Multizone, 0)
	for i = 1:nummzs
		uzs = Array(UnsaturatedZone, 0)
		szs = Array(SaturatedZone, 0)
		j = 1
		while haskey(params, "p$(i)_uz$(j)_v")
			push!(uzs, UnsaturatedZone(params["p$(i)_uz$(j)_v"], params["p$(i)_uz$(j)_dispersivity"], params["p$(i)_uz$(j)_drippoint"]))
			x0 = getparamvec(params, "p$(i)_sz$(j)_x0_")
			sigma0 = getparamvec(params, "p$(i)_sz$(j)_sigma0_")
			v = getparamvec(params, "p$(i)_sz$(j)_v_")
			dispersivity = getparamvec(params, "p$(i)_sz$(j)_dispersivity_")
			push!(szs, SaturatedZone(x0, sigma0, v, dispersivity, params["p$(i)_sz$(j)_drippoint"], eval(params["p$(i)_sz$(j)_anasolfunc"])))
			j += 1
		end
		mzs[i] = Multizone(uzs, szs, params["p$(i)_releasemass"], params["p$(i)_maxtime"], numsamples)
	end
	return mzs
end

function madsforward(madsdata, parameters)
	obs = madsdata["Observations"]
	nummzs = madsdata["Pathways"]
	numsamples = madsdata["Time samples"]
	mzs = assemblemzs(madsdata, nummzs, numsamples)
	results = Dict()
	for obsname in obs
		thisobs = obs[obsname]
		if startswith(obsname, "mzt_")
			results[obsname] = 0.
			for i = 1:length(mzs)
				szindex = thisobs["szindices"][i]
				if szindex != false
					results[obsname] += getconcentration(mzs[i], obs[obsname]["x"], obs[obsname]["time"], obs[szindex])
				end
			end
		end
	end
	return results
end
