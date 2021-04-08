using ClimateERA
using Dates
using Dierckx
using GeoRegions
using Logging
using NCDatasets
using NumericalIntegration

include(srcdir("common.jl"))

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
    cmod,cpar,_,_ = erainitialize(init,modID="msfc",parID="tcwv");
    smod,spar,_,_ = erainitialize(init,modID="csfc",parID="swp");
    global_logger(ConsoleLogger(stdout,Logging.Info))

    for dtii in datevec

        cds,cvar = erarawread(cmod,cpar,ereg,eroot,dtii);
        sds,svar = erarawread(smod,spar,ereg,eroot,dtii);

        @info "$(Dates.now()) - Calculating $(uppercase(emod["dataset"])) $(epar["name"]) data in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(dtii)) $(Dates.monthname(dtii)) ..."

        tcwv = cvar[:] * 1
        swp  = svar[:] * 1000
        csf  = tcwv ./ swp * 100

        close(cds); close(sds);
        erarawsave(csf,emod,epar,ereg,dtii,sroot)

    end

    putinfo(emod,epar,ereg,etime,sroot);

end
