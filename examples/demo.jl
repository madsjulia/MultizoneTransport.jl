import MultizoneTransport
import PyPlot
import Anasol

function demo()
	dispersivity = 10.
	drippoint = 100.
	uz1 = MultizoneTransport.UnsaturatedZone(25, dispersivity, 300.)
	#sz1 = MultizoneTransport.SaturatedZone([0.;], [sqrt(300.);], [10.;], [dispersivity;], drippoint, Anasol.long_b_d_i_cf)
	sz1 = MultizoneTransport.SaturatedZone([0., 0., 0.], [sqrt(300.), sqrt(30.), sqrt(.01)], [10., 0., 0.], [drippoint, NaN, NaN], [dispersivity, .1 * dispersivity, 0.01 * dispersivity], drippoint, Anasol.long_bbb_ddd_aii_cf)
	uz2 = MultizoneTransport.UnsaturatedZone(25., dispersivity, 100.)
	sz2 = MultizoneTransport.SaturatedZone([0., 0., 0.], [sqrt(300.), sqrt(30.), sqrt(.01)], [20, 0., 0.], [NaN, NaN, NaN], [dispersivity, .1 * dispersivity, 0.01 * dispersivity], drippoint, Anasol.long_bbb_ddd_iii_cf)
	#sz2 = MultizoneTransport.SaturatedZone([0.;], [sqrt(300.);], [50.;], [dispersivity;], drippoint, Anasol.long_b_d_i_cf)
	massrelease = 31.4159e3
	mz = MultizoneTransport.Multizone([uz1, uz2], [sz1, sz2], massrelease, 50., 10000)
	zones1 = [uz1;]
	zones2 = [uz1, sz1, uz2]
	maxtime = 50
	bs1 = MultizoneTransport.getbrowniansource(zones1, maxtime, 10000)
	bs2 = MultizoneTransport.getbrowniansource(zones2, maxtime, 10000)
	ts = collect(linspace(0, maxtime, 101))
	xs = collect(linspace(-20, drippoint, 100))
	pdfs1 = massrelease * map(bs1, ts)
	pdfs2 = massrelease * map(bs2, ts)
	modeltime = 0.
	for i = 1:length(ts)
		t = ts[i]
		print("model run time:")
		modeltime += @elapsed @time begin
			#csupper = map(x->MultizoneTransport.getconcentration(mz, [x;], t, 1), xs)
			csupper = map(x->MultizoneTransport.getconcentration(mz, [x, 0., 0.], t, 1), xs)
			#cslower = map(x->MultizoneTransport.getconcentration(mz, [x;], t, 2), xs)
			cslower = map(x->MultizoneTransport.getconcentration(mz, [x, 0., 0.], t, 2), xs)
		end
		print("plot run time:")
		@time begin
			PyPlot.clf()
			PyPlot.subplot(2, 2, 1)
			PyPlot.plot(ts[1:i], pdfs1[1:i], "k", lw=3)
			PyPlot.xlabel("time")
			PyPlot.ylabel("Mass flux into upper zone")
			PyPlot.xlim(0, maximum(ts))
			PyPlot.ylim(0, 5e3)
			PyPlot.subplot(2, 2, 3)
			PyPlot.plot(ts[1:i], pdfs2[1:i], "k", lw=3)
			PyPlot.xlabel("time")
			PyPlot.ylabel("Mass flux into lower zone")
			PyPlot.xlim(0, maximum(ts))
			PyPlot.ylim(0, 5e3)
			PyPlot.subplot(2, 2, 2)
			PyPlot.plot(xs, csupper, "k", lw=3)
			PyPlot.xlabel("x")
			PyPlot.ylabel("Concentration in upper zone")
			PyPlot.ylim(0, 5e0)
			PyPlot.subplot(2, 2, 4)
			PyPlot.plot(xs, cslower, "k", lw=3)
			PyPlot.xlabel("x")
			PyPlot.ylabel("Concentration in lower zone")
			PyPlot.ylim(0, 5e0)
			PyPlot.plt[:tight_layout]()
		end
	end
end

demo()
