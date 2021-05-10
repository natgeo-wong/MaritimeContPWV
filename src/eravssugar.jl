using ClimateERA
using Dates
using JLD2
using NPZ
using StatsBase
using PyCall

np = pyimport("numpy")

function eravsugar(
    init::AbstractDict, croot::AbstractDict;
    gpwv::AbstractRange=0:100, epwv::AbstractRange=0:100
)

    cmod,cpar,creg,etime = erainitialize(init,modID="msfc",parID="tcw",regID="SMT")

    gdata = npzread(datadir("ssugar/data_sugar_3h_45stations.npy"))
	gdata = dropdims(mean(reshape(gdata,8,:,45),dims=1),dims=1)

	npar = Array{Any,2}(undef,45,4)
	vars = np.load(datadir("ssugar/sugar_stations1.npy"))
	for iy = 1 : 4, ix = 1 : 45
		npar[ix,iy] = vars[ix][iy]
	end
	npar[:,2:4] .= parse.(Float64,npar[:,2:4])

    ngps = length(gpwv); nera = length(epwv);
    gve  = zeros(Int64,ngps-1,nera-1)
	it = 1
	rms = 0

    @info "$(Dates.now()) - Comparing GNSS ZWD vs ECMWF reanalysis TCW"

    for yr = 2008 : 2013, mo = 1 : 12

		cds,cvar = erarawread(cmod,cpar,creg,croot,Date(yr,mo))
		nt = daysinmonth(Date(yr,mo))

	    @info "$(Dates.now()) - Analysis proceeding for $yr $(monthname(mo))"

        for istn = 1 : 45

			glon = npar[istn,4]
			glat = npar[istn,3]
			ilon,ilat = regionpoint(glon,glat,creg["lon"],creg["lat"])

            etcw = cvar[ilon,ilat,:]*1
			etcw = dropdims(mean(reshape(etcw,hrindy(cmod),:),dims=1),dims=1)

			gtcw = gdata[it:(it+nt-1),istn]
			etcw = etcw[.!isnan.(gtcw)]
			gtcw = gtcw[.!isnan.(gtcw)]

			if !isempty(gtcw)
            	gve += fit(Histogram,(gtcw,etcw),(gpwv,epwv)).weights
				rms += sum((gtcw.-etcw).^2)
			end

        end

		close(cds)
		it += nt

    end

	rms = sqrt(rms/sum(gve))
    @save "$(datadir("compiled/ssugar.jld2"))" gve rms

end
