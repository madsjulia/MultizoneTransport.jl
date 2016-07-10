import Mads
import MultizoneTransport

md = Mads.loadmadsfile("test.mads")
forward_results = Mads.forward(md)
Mads.setobservationtargets!(md, forward_results) # set observation targets to match the forward run results
Mads.plotmatches(md, forward_results) # generates 'test-match.svg'
run(`open test-match.svg`) # works only on mac os x; open by mouse clicking if not on mac
Mads.savemadsfile(md, "test-new.mads")