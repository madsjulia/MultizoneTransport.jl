import MultizoneTransport
import Mads
import Anasol

function mzmads(madsdata, parameters)
	return MultizoneTransport.madsforward(madsdata, parameters)
end

md = Mads.loadyamlmadsfile("test.mads")
x = mzmads(md, Dict(zip(Mads.getparamkeys(md), Mads.getparamsinit(md))))
