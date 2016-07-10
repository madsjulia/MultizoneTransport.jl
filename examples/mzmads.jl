function mzmads(madsdata)
	function mzforward(parameters)
		return MultizoneTransport.madsforward(madsdata, parameters)
	end
	return mzforward
end

