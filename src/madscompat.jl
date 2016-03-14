import Mads

function getparamvec(params, prefix)
	i = Array(Int, 0)
	p = Array(Float64, 0)
	for k in keys(params)
		if startswith(k, prefix)
			push!(i, parse(Int, k[length(prefix) + 1:end]))
			push!(p, params[k])
		end
	end
	pvec = Array(Float64, length(i))
	for j = 1:length(i)
		pvec[i[j]] = p[j]
	end
	return pvec
end

function assemblemzs(params, nummzs, numsamples, madsdata)
	mzs = Array(Multizone, nummzs)
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
			push!(szs, SaturatedZone(x0, sigma0, v, [params["p$(i)_sz$(j)_drippoint"], NaN, NaN], dispersivity, eval(parse(madsdata["p$(i)_sz$(j)_anasolfunc"]))))
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
	mzs = assemblemzs(parameters, nummzs, numsamples, madsdata)
	results = Dict()
	for wellkey in Mads.getwellkeys(madsdata)
		if madsdata["Wells"][wellkey]["on"]
			local wellx0::Array{Float64, 1}
			local wellx1::Array{Float64, 1}
			if haskey(madsdata["Wells"][wellkey], "z0")
				wellx0 = [madsdata["Wells"][wellkey]["x"], madsdata["Wells"][wellkey]["y"], madsdata["Wells"][wellkey]["z0"]]
				wellx1 = [madsdata["Wells"][wellkey]["x"], madsdata["Wells"][wellkey]["y"], madsdata["Wells"][wellkey]["z1"]]
				if abs( wellx1[3] - wellx0[3] ) > 0.1
					screen = true
				else
					screen = false
				end
			elseif haskey(madsdata["Wells"][wellkey], "y")
				wellx0 = [madsdata["Wells"][wellkey]["x"], madsdata["Wells"][wellkey]["y"]]
				wellx1 = [madsdata["Wells"][wellkey]["x"], madsdata["Wells"][wellkey]["y"]]
				screen = false
			else
				wellx0 = Float64[madsdata["Wells"][wellkey]["x"]]
				wellx1 = Float64[madsdata["Wells"][wellkey]["x"]]
				screen = false
			end
			for o in 1:length(madsdata["Wells"][wellkey]["obs"])
				t = madsdata["Wells"][wellkey]["obs"][o]["t"]
				conc = 0.
				for i = 1:length(mzs)
					if screen
						conc += .5 * getconcentration(mzs[i], wellx0, t, madsdata["Wells"][wellkey]["szindex"])
						conc += .5 * getconcentration(mzs[i], wellx1, t, madsdata["Wells"][wellkey]["szindex"])
					else
						conc += .5 * getconcentration(mzs[i], .5 * (wellx1 + wellx0), t, madsdata["Wells"][wellkey]["szindex"])
					end
				end
				results[string(wellkey, "_", t)] = conc
			end
		end
	end
	#=
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
	=#
	return results
end
