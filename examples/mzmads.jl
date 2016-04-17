import MultizoneTransport
import Anasol

function mzmads(madsdata)
	function mzforward(parameters)
		x = MultizoneTransport.madsforward(madsdata, parameters)
		return x
	end
	return mzforward
end

