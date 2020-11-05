using ClimateERA
using Dates
using Dierckx
using GeoRegions
using Logging
using NCDatasets
using NumericalIntegration

include(srcdir("common.jl"))

function t2esat(T::Real)

    # e(T=273.16) = 611.657 Pa
    # L/R_v = 2.5009e6 / 461.51 = 5418.950835301511
    # 1/273.16 = 3.6608581051398447e-3

    return 611.657 * exp(5418.951 * (3.660858e-3 - 1/T))

end

function calccsf(
    rhum::AbstractVector{<:Real},
    pre::AbstractVector{<:Real},
    psfc::Real
)

    r = cumul_integrate(pre,rhum) ./ cumul_integrate(pre); r[1] = 0
    if psfc < 1000; pre[end] = 101235 end
    spl = Spline1D(pre,r,k=1); return spl(psfc)

end

function csf(
    init::Dict,eroot::Dict,sroot::Dict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    gres::Real=0
)

    emod,epar,ereg,etime = erainitialize(
        init;
        modID="csfc",parID="csf",regID=regID,timeID=timeID,
        gres=gres
    )

    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    pmod,ppar,_,_ = erainitialize(init,modID="dsfc",parID="p_sfc");
    tmod,tpar,_,_ = erainitialize(init,modID="dsfc",parID="t_sfc");
    dmod,dpar,_,_ = erainitialize(init,modID="msfc",parID="t_dew");
    rmod,rpar,_,_ = erainitialize(init,modID="mpre",parID="rhum");
    global_logger(ConsoleLogger(stdout,Logging.Info))

    ehr = hrindy(emod); nlon,nlat = ereg["size"]
    p = ClimateERA.erapressureload(); p = p[p.>=10]*100; np = length(p)
    p = convert.(Float32,vcat(0,p,0))
    rh = Array{Float32,3}(undef,nlon,nlat,np)
    ps = Array{Float32,2}(undef,nlon,nlat)
    Ts = Array{Float32,2}(undef,nlon,nlat)
    Td = Array{Float32,2}(undef,nlon,nlat)

    rhii = Vector{Float32}(undef,np+2)
    rhii[1] = 0

    for dtii in datevec

        @info "$(Dates.now()) - Preallocating arrays for $(uppercase(emod["dataset"])) $(epar["name"]) data in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(dtii)) $(Dates.monthname(dtii)) ..."
        nhr = ehr * daysinmonth(dtii); csf = zeros(nlon,nlat,nhr);

        @info "$(Dates.now()) - Calculating $(uppercase(emod["dataset"])) $(epar["name"]) data in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(dtii)) $(Dates.monthname(dtii)) ..."

        pds,pvar = erarawread(pmod,ppar,ereg,eroot,dtii);
        tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii);
        dds,dvar = erarawread(dmod,dpar,ereg,eroot,dtii);

        for it = 1 : nhr

            ps .= pvar[:,:,it]*1
            Ts .= tvar[:,:,it]*1
            Td .= dvar[:,:,it]*1

            for ip = 1 : np; pre = Int16(p[ip+1]/100);
                rpar["level"] = pre; rds,rvar = erarawread(rmod,rpar,ereg,eroot,dtii);
                rh[:,:,ip] .= rvar[:,:,it]*1; close(rds);
            end

            for ilat = 1 : nlat, ilon = 1 : nlon

                for ip = 1 : np; rhii[ip+1] = rh[ilon,ilat,ip] end
                p[end] = ps[ilon,ilat]
                rhii[end] = t2esat(Td[ilon,ilat]) / t2esat(Ts[ilon,ilat])
                csf[ilon,ilat,it] = calccsf(rhii,p,ps[ilon,ilat])

            end

        end

        close(pds); close(tds); close(dds);

    end

    putinfo(emod,epar,ereg,etime,sroot);

end
