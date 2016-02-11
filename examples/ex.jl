import MultizoneTransport
import PyPlot
import Anasol

function onedimex()
	dispersivity = 10.
	uz1 = MultizoneTransport.UnsaturatedZone(25, dispersivity, 300.)
	sz1 = MultizoneTransport.SaturatedZone([0.;], [sqrt(300.);], [50.;], [dispersivity;], 100., Anasol.long_b_b_i_cf)
	uz2 = MultizoneTransport.UnsaturatedZone(25., dispersivity, 25.)
	sz2 = MultizoneTransport.SaturatedZone([0.;], [sqrt(300.);], [50.;], [dispersivity;], 100., Anasol.long_b_b_i_cf)
	massrelease = 31.4159
	mz = MultizoneTransport.Multizone([uz1, uz2], [sz1, sz2], massrelease, 50., 10000)
	zones1 = [uz1;]
	zones2 = [uz1, sz1, uz2]
	maxtime = 50
	bs1 = MultizoneTransport.getbrowniansource(zones1, maxtime, 10000)
	bs2 = MultizoneTransport.getbrowniansource(zones2, maxtime, 10000)
	ts = linspace(0, maxtime, 1000)
	@time csupper = map(t->MultizoneTransport.getconcentration(mz, [-2.;], t, 1), ts)
	@time cslower = map(t->MultizoneTransport.getconcentration(mz, [50.;], t, 2), ts)
	pdfs1 = map(bs1, ts)
	pdfs2 = map(bs2, ts)
	PyPlot.clf()
	PyPlot.plot(ts, pdfs1, label="upper sz flux / mass")
	PyPlot.plot(ts, pdfs2, label="lower sz flux / mass")
	PyPlot.plot(ts, csupper, label="upper sz conc")
	PyPlot.plot(ts, cslower, label="lower sz conc")
	PyPlot.legend()
end

function twodimex()
	dispersivity = 10.
	uz1 = MultizoneTransport.UnsaturatedZone(25, dispersivity, 300.)
	sz1 = MultizoneTransport.SaturatedZone([0., 0., 0.], [sqrt(300.), sqrt(300.)], [50., 0.], [dispersivity, 0.1 * dispersivity], 100., Anasol.long_bb_bb_ii_cf)
	uz2 = MultizoneTransport.UnsaturatedZone(25., dispersivity, 25.)
	sz2 = MultizoneTransport.SaturatedZone([0., 0., 0.], [sqrt(300.), sqrt(300.)], [50., 0.], [dispersivity, 0.1 * dispersivity], 100., Anasol.long_bb_bb_ii_cf)
	massrelease = 31.4159e2
	mz = MultizoneTransport.Multizone([uz1, uz2], [sz1, sz2], massrelease, 50., 10000)
	zones1 = [uz1;]
	zones2 = [uz1, sz1, uz2]
	maxtime = 50
	bs1 = MultizoneTransport.getbrowniansource(zones1, maxtime, 10000)
	bs2 = MultizoneTransport.getbrowniansource(zones2, maxtime, 10000)
	ts = linspace(0, maxtime, 1000)
	@time csupper = map(t->MultizoneTransport.getconcentration(mz, [50., 0.], t, 1), ts)
	@time cslower = map(t->MultizoneTransport.getconcentration(mz, [50., 0.], t, 2), ts)
	pdfs1 = map(bs1, ts)
	pdfs2 = map(bs2, ts)
	PyPlot.clf()
	PyPlot.plot(ts, pdfs1, label="upper sz flux / mass")
	PyPlot.plot(ts, pdfs2, label="lower sz flux / mass")
	PyPlot.plot(ts, csupper, label="upper sz conc")
	PyPlot.plot(ts, cslower, label="lower sz conc")
	PyPlot.legend()
end

function threedimex()
	dispersivity = 10.
	uz1 = MultizoneTransport.UnsaturatedZone(25, dispersivity, 300.)
	sz1 = MultizoneTransport.SaturatedZone([0., 0., 0.], [sqrt(300.), sqrt(300.), sqrt(.1)], [50., 0., 0.], [dispersivity, 0.1 * dispersivity, 0.01 * dispersivity], 100., Anasol.long_bbb_bbb_iii_cf)
	uz2 = MultizoneTransport.UnsaturatedZone(25., dispersivity, 25.)
	sz2 = MultizoneTransport.SaturatedZone([0., 0., 0.], [sqrt(300.), sqrt(300.), sqrt(.1)], [50., 0., 0.], [dispersivity, 0.1 * dispersivity, 0.01 * dispersivity], 100., Anasol.long_bbb_bbb_iii_cf)
	massrelease = 31.4159e3
	mz = MultizoneTransport.Multizone([uz1, uz2], [sz1, sz2], massrelease, 50., 10000)
	zones1 = [uz1;]
	zones2 = [uz1, sz1, uz2]
	maxtime = 50
	bs1 = MultizoneTransport.getbrowniansource(zones1, maxtime, 10000)
	bs2 = MultizoneTransport.getbrowniansource(zones2, maxtime, 10000)
	ts = linspace(0, maxtime, 1000)
	@time csupper = map(t->MultizoneTransport.getconcentration(mz, [50., 0., 0.], t, 1), ts)
	@time cslower = map(t->MultizoneTransport.getconcentration(mz, [50., 0., 0.], t, 2), ts)
	pdfs1 = map(bs1, ts)
	pdfs2 = map(bs2, ts)
	PyPlot.clf()
	PyPlot.plot(ts, pdfs1, label="upper sz flux / mass")
	PyPlot.plot(ts, pdfs2, label="lower sz flux / mass")
	PyPlot.plot(ts, csupper, label="upper sz conc")
	PyPlot.plot(ts, cslower, label="lower sz conc")
	PyPlot.legend()
end

PyPlot.figure(1)
onedimex()
PyPlot.figure(2)
twodimex()
PyPlot.figure(3)
threedimex()
