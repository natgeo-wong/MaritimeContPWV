using ClimateERA
using Dierckx
using NumericalIntegration

function t2esat(T::Real)

    # e(T=273.16) = 611.657 Pa
    # L/R_v = 2.5009e6 / 461.51 = 5418.950835301511
    # 1/273.16 = 3.6608581051398447e-3

    return 611.657 * exp(5418.951 * (3.660858e-3 - 1/T))

end

function calcswp(
    esat::Vector{<:Real}, t::Vector{<:Real}, p::Vector{<:Real}, ps::Real
)

    esat[2:end] .= t[2:end] .* esat[2:end] ./ p[2:end]
    svec = 2.925586e-2 * cumul_integrate(p,esat)
    if ps < p[end-1]; p = vcat(p[1:(end-1)],1012.35); end
    spl = Spline1D(p,svec,k=1); return spl(ps)

end

function tairread(emod::Dict, epar::Dict, ereg::Dict, eroot::Dict, date::TimeType)

    fol = joinpath(eroot["era"],ereg["fol"],"t_air","raw",yr);
    fnc = "$(emod["prefix"])-$(ereg["fol"])-t_air-sfc-$(yr2str(date))"
    eds = Dataset(joinpath(fol,fnc))
    if haskey(eds,epar["ID"]); ID = epar["ID"]; else; ID = epar["IDnc"]; end
    return eds,eds[ID]

end

function swp(
    init::Dict,eroot::Dict,sroot::Dict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    gres::Real=0
)

    emod,epar,ereg,etime = erainitialize(
        init;
        modID="csfc",parID="swp",regID=regID,
        gres=gres
    )

    ehr = hrindy(emod); nlon,nlat = ereg["size"]
    p = ClimateERA.erapressureload()*100; np = length(p); p = vcat(0,p,0)
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    smod,spar,_,_ = erainitialize(init,modID="dsfc",parID="t_sfc");
    pmod,ppar,_,_ = erainitialize(init,modID="dsfc",parID="p_sfc");
    tmod,tpar,_,_ = erainitialize(init,modID="dpre",parID="t_air");
    global_logger(ConsoleLogger(stdout,Logging.Info))

    Ts = Array{Float32,3}(undef,nlon,nlat)
    ps = Array{Float32,3}(undef,nlon,nlat)
    Ta = Array{Float32,3}(undef,nlon,nlat,np); Taii = Vector{Float32}(undef,np+2)
    esii = Vector{Float32}(undef,np+2)

    for dtii in datevec

        @info "$(Dates.now()) - Preallocating arrays ..."
        nhr = ehr * daysinmonth(dtii); swp = zeros(nlon,nlat,nhr);

        @info "$(Dates.now()) - Calculating Saturated Vapour Path data for $(dtii) ..."
        for it = 1 : nhr

            sds,svar = erarawread(smod,spar,ereg,eroot,dtii);
            pds,pvar = erarawread(pmod,ppar,ereg,eroot,dtii);
            tds,tvar = tairread(tmod,tpar,ereg,eroot,dtii);
            Ts .= svar[:,:,it]*1;   close(sds);
            ps .= pvar[:,:,it]*1;   close(pds);
            Ta .= tvar[:,:,:,it]*1; close(tds);

            for ilat = 1 : nlat, ilon = 1 : nlon

                for ip = 1 : np; Taii[ip+1] = Ta[ilon,ilat,ip] end
                Taii[end] = Ts[ilon,ilat]; p[end] = ps[ilon,ilat];
                esii[2:end] .= t2esat(Taii[2:end])
                swp[ilon,ilat,it] = calcswp(esii,Taii,p,ps[ilon,ilat]);

            end

        end

        @info "$(Dates.now()) - Saving Saturated Vapour Path data for $(dtii) ..."
        erarawsave(swp,emod,epar,ereg,dtii,sroot)

    end

    putinfo(emod,epar,ereg,etime,sroot);

end

function csf(
    init::Dict,eroot::Dict,sroot::Dict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    gres::Real=0
)

    emod,epar,ereg,etime = erainitialize(
        init;
        modID="csfc",parID="csf",regID=regID,
        gres=gres
    )

    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    tmod,tpar,_,_ = erainitialize(init,modID="msfc",parID="tcwv");
    smod,spar,_,_ = erainitialize(init,modID="csfc",parID="swp");
    global_logger(ConsoleLogger(stdout,Logging.Info))

    for dtii in datevec

        @info "$(Dates.now()) - Calculating Column Saturation Fraction data for $(dtii) ..."
        tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii);
        sds,svar = erarawread(smod,spar,ereg,eroot,dtii);
        tcw = tvar[:,:,it]*1; close(tds);
        swp = svar[:,:,it]*1; close(sds);
        csf = tcw ./ swp

        @info "$(Dates.now()) - Saving Column Saturation Fraction data for $(dtii) ..."
        erarawsave(csf,emod,epar,ereg,dtii,sroot)

    end

    putinfo(emod,epar,ereg,etime,sroot);

end
