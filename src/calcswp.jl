using ClimateERA
using Dates
using Dierckx
using GeoRegions
using Logging
using NCDatasets
using NumericalIntegration

include(srcdir("common.jl"))

function t2esat(T::Real,P::Real)

    tb = T - 273.15
    if tb <= 0
    	esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
    	esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end

    r = 0.622 * esat / max(esat,P-esat)
    return r / (1+r)

end

function calcswp(
    esat::AbstractVector{<:Real},
    pre::AbstractVector{<:Real},
    psfc::Real
)

    r = cumul_integrate(pre,esat) / 9.81
    r[1] = 0
    if psfc <= 100000; pre[end] = 101235 end
    spl = Spline1D(pre,r,k=1); return spl(psfc)

end

function swp(
    init::Dict,eroot::Dict,sroot::Dict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    gres::Real=0
)

    emod,epar,ereg,etime = erainitialize(
        init;
        modID="csfc",parID="swp",regID=regID,timeID=timeID,
        gres=gres
    )

    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    pmod,ppar,_,_ = erainitialize(init,modID="dsfc",parID="p_sfc");
    smod,spar,_,_ = erainitialize(init,modID="dsfc",parID="t_sfc");
    tmod,tpar,_,_ = erainitialize(init,modID="dpre",parID="t_air");
    global_logger(ConsoleLogger(stdout,Logging.Info))

    ehr = hrindy(emod); nlon,nlat = ereg["size"]
    p = ClimateERA.erapressureload(); p = p[p.>=10]*100; np = length(p)
    p = convert.(Float32,vcat(0,p,0))
    Ta = Array{Float32,3}(undef,nlon,nlat,np)
    ps = Array{Float32,2}(undef,nlon,nlat)
    Ts = Array{Float32,2}(undef,nlon,nlat)

    esat = Vector{Float32}(undef,np+2)
    esat[1] = 0

    for dtii in datevec

        @info "$(Dates.now()) - Preallocating arrays for $(uppercase(emod["dataset"])) $(epar["name"]) data in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(dtii)) $(Dates.monthname(dtii)) ..."
        nhr = ehr * daysinmonth(dtii); swp = zeros(nlon,nlat,nhr);

        pds,pvar = erarawread(pmod,ppar,ereg,eroot,dtii);
        sds,svar = erarawread(smod,spar,ereg,eroot,dtii);

        @info "$(Dates.now()) - Calculating $(uppercase(emod["dataset"])) $(epar["name"]) data in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(dtii)) $(Dates.monthname(dtii)) ..."

        for it = 1 : nhr

            ps .= pvar[:,:,it]*1
            Ts .= svar[:,:,it]*1

            for ip = 1 : np; pre = Int16(p[ip+1]/100);
                tpar["level"] = pre; tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii);
                Ta[:,:,ip] .= tvar[:,:,it]*1; close(tds);
            end

            for ilat = 1 : nlat, ilon = 1 : nlon

                for ip = 1 : np; esat[ip+1] = t2esat(Ta[ilon,ilat,ip],p[ip+1]) end
                p[end] = ps[ilon,ilat]
                esat[end] = t2esat(Ts[ilon,ilat],p[end])
                swp[ilon,ilat,it] = calcswp(esat,p,ps[ilon,ilat])

            end

        end

        close(pds); close(sds);
        erarawsave(swp,emod,epar,ereg,dtii,sroot)

    end

    putinfo(emod,epar,ereg,etime,sroot);

end
