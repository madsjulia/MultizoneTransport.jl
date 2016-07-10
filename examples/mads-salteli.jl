import Mads
import MultizoneTransport

md = Mads.loadmadsfile("test.mads")
root = Mads.getmadsrootname(md)
forward_results = Mads.forward(md)
Mads.setobservationtargets!(md, forward_results) # set observation targets to match the forward run results
Mads.plotmatches(md, forward_results) # generates 'test-match.svg'
# run(`open test-match.svg`) # works only on mac os x; open by mouse clicking if not on mac
Mads.savemadsfile(md, "test-new.mads")
nsamples = 10000000
sens_results = Mads.saltelli(md, N=nsamples)
JLD.save(root * "_saltelli_$nsamples.jld", sens_results)
Mads.plotwellSAresults(md, sens_results)