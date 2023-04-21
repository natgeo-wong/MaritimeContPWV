using ClimateERA
using GeoRegions
using JLD2
using StatsBase

include(srcdir("gnss.jl"))

function eravgnss(
    gstn::AbstractString,
    init::AbstractDict, croot::AbstractDict, proot::AbstractDict;
    gpwv::AbstractRange=0:100, epwv::AbstractRange=0:100
)

    cmod,cpar,creg,etime = erainitialize(init,modID="msfc",parID="tcw",regID="SMT")

    if init["datasetID"] == 1
        pmod,ppar,preg,_ = erainitialize(
            init,
            modID="csfc",parID="Pi_RE5",regID="SMT",
            gres=1
        )
    else
        pmod,ppar,preg,_ = erainitialize(
            init,
            modID="csfc",parID="Pi_REI",regID="SMT",
            gres=1
        )
    end

    ginf = gstationinfo(gstn,retrieveginfo());
    glon = ginf["longitude"]; glat = ginf["latitude"]

    ilon,ilat = regionpoint(glon,glat,creg["lon"],creg["lat"])
    plon,plat = regionpoint(glon,glat,preg["lon"],preg["lat"])

    ngps = length(gpwv); nera = length(epwv); nt = (etime["End"] - etime["Begin"] + 1) * 12
    gve  = zeros(Int64,ngps-1,nera-1,nt); ii = 0; nhr = hrstep(cmod)

    @info "$(Dates.now()) - Comparing GNSS ZWD vs ECMWF reanalysis TCW at Station $gstn"

    for yr = etime["Begin"] : etime["End"], mo = 1 : 12

        date = Date(yr,mo)
        gfol = datadir("gnss/$(gstn)/$(yr2str(date))")
        gfnc = joinpath("$gfol","$(gstn)-$(yrmo2str(date)).nc")

        if isfile(gfnc); ii += 1

            gds = Dataset(gfnc); gzwd = gds["zwd"][:]*1; close(gds);
            gzwd = reshape(gzwd,6*nhr,:);
            gzwd = dropdims(mean(gzwd,dims=1),dims=1);

            pds,pvar = erarawread(pmod,ppar,preg,proot,date)
            ipi  = pvar[plon,plat,:]*1; close(pds)
            gtcw = gzwd .* ipi * 1000

            cds,cvar = erarawread(cmod,cpar,creg,croot,date)
            etcw = cvar[ilon,ilat,:]*1; close(cds)

            irem = (.!ismissing.(gzwd)) .& (.!isnan.(gtcw));
            etcw = etcw[irem]; gtcw = gtcw[irem]
            gve[:,:,ii] = fit(Histogram,(gtcw,etcw),(gpwv,epwv)).weights

        end

    end

    gve = @view gve[:,:,1:ii]; gve = dropdims(sum(gve,dims=3),dims=3)
    @save "$(datadir("compiled/$(init["prefix"])/$(gstn).jld2"))" gve

end

function era5vgnss_6hourly(
    gstn::AbstractString,
    init::AbstractDict, croot::AbstractDict, proot::AbstractDict;
    gpwv::AbstractRange=0:100, epwv::AbstractRange=0:100
)

    cmod,cpar,creg,etime = erainitialize(init,modID="msfc",parID="tcw",regID="SMT")

    if init["datasetID"] != 1
        error("$(Dates.now()) - era5vgnss_6hourly() is only for ERA5 comparison")
    else
        pmod,ppar,preg,_ = erainitialize(
            init,
            modID="csfc",parID="Pi_RE5",regID="SMT",
            gres=1
        )
    end

    ginf = gstationinfo(gstn,retrieveginfo());
    glon = ginf["longitude"]; glat = ginf["latitude"]

    ilon,ilat = regionpoint(glon,glat,creg["lon"],creg["lat"])
    plon,plat = regionpoint(glon,glat,preg["lon"],preg["lat"])

    ngps = length(gpwv); nera = length(epwv); nt = (etime["End"] - etime["Begin"] + 1) * 12
    gve  = zeros(Int64,ngps-1,nera-1,nt); ii = 0; nhr = hrstep(cmod)

    @info "$(Dates.now()) - Comparing GNSS ZWD vs ECMWF reanalysis TCW at Station $gstn"

    for yr = etime["Begin"] : etime["End"], mo = 1 : 12

        date = Date(yr,mo)
        gfol = datadir("gnss/$(gstn)/$(yr2str(date))")
        gfnc = joinpath("$gfol","$(gstn)-$(yrmo2str(date)).nc")

        if isfile(gfnc); ii += 1

            gds = Dataset(gfnc); gzwd = gds["zwd"][:]*1; close(gds);
            gzwd = reshape(gzwd,36,:);
            gzwd = dropdims(mean(gzwd,dims=1),dims=1);

            pds,pvar = erarawread(pmod,ppar,preg,proot,date)
            ipi  = pvar[plon,plat,:]*1; close(pds)
            ipi  = reshape(ipi,6,:)
            ipi  = dropdims(mean(ipi,dims=1),dims=1);
            gtcw = gzwd .* ipi * 1000

            cds,cvar = erarawread(cmod,cpar,creg,croot,date)
            etcw = cvar[ilon,ilat,:]*1; close(cds)
            etcw = reshape(etcw,6,:)
            etcw = dropdims(mean(etcw,dims=1),dims=1);

            irem = (.!ismissing.(gzwd)) .& (.!isnan.(gtcw));
            etcw = etcw[irem]; gtcw = gtcw[irem]
            gve[:,:,ii] = fit(Histogram,(gtcw,etcw),(gpwv,epwv)).weights

        end

    end

    gve = @view gve[:,:,1:ii]; gve = dropdims(sum(gve,dims=3),dims=3)
    @save "$(datadir("compiled/era6/$(gstn).jld2"))" gve

end

function era5vgnss_daily(
    gstn::AbstractString,
    init::AbstractDict, croot::AbstractDict, proot::AbstractDict;
    gpwv::AbstractRange=0:100, epwv::AbstractRange=0:100
)

    cmod,cpar,creg,etime = erainitialize(init,modID="msfc",parID="tcw",regID="SMT")

    if init["datasetID"] != 1
        error("$(Dates.now()) - era5vgnss_daily() is only for ERA5 comparison")
    else
        pmod,ppar,preg,_ = erainitialize(
            init,
            modID="csfc",parID="Pi_RE5",regID="SMT",
            gres=1
        )
    end

    ginf = gstationinfo(gstn,retrieveginfo());
    glon = ginf["longitude"]; glat = ginf["latitude"]

    ilon,ilat = regionpoint(glon,glat,creg["lon"],creg["lat"])
    plon,plat = regionpoint(glon,glat,preg["lon"],preg["lat"])

    ngps = length(gpwv); nera = length(epwv); nt = (etime["End"] - etime["Begin"] + 1) * 12
    gve  = zeros(Int64,ngps-1,nera-1,nt); ii = 0; nhr = hrstep(cmod)

    @info "$(Dates.now()) - Comparing GNSS ZWD vs ECMWF reanalysis TCW at Station $gstn"

    for yr = etime["Begin"] : etime["End"], mo = 1 : 12

        date = Date(yr,mo)
        gfol = datadir("gnss/$(gstn)/$(yr2str(date))")
        gfnc = joinpath("$gfol","$(gstn)-$(yrmo2str(date)).nc")

        if isfile(gfnc); ii += 1

            gds = Dataset(gfnc); gzwd = gds["zwd"][:]*1; close(gds);
            gzwd = reshape(gzwd,144,:); ndy = size(gzwd,2)
            gtmp = zeros(ndy)
            for idy = 1 : ndy
                gii = @view gzwd[:,idy]
                gii = @view gii[.!isnan.(gii)]; nnan = length(gii)
                if nnan >= 72
                      gtmp[idy] = mean(gii)
                else; gtmp[idy] = NaN
                end
            end

            pds,pvar = erarawread(pmod,ppar,preg,proot,date)
            ipi  = pvar[plon,plat,:]*1; close(pds)
            ipi  = reshape(ipi,24,:)
            ipi  = dropdims(mean(ipi,dims=1),dims=1);
            gtcw = gtmp .* ipi * 1000

            cds,cvar = erarawread(cmod,cpar,creg,croot,date)
            etcw = cvar[ilon,ilat,:]*1; close(cds)
            etcw = reshape(etcw,24,:)
            etcw = dropdims(mean(etcw,dims=1),dims=1);

            irem = (.!ismissing.(gtmp)) .& (.!isnan.(gtmp));
            etcw = etcw[irem]; gtcw = gtcw[irem]
            gve[:,:,ii] = fit(Histogram,(gtcw,etcw),(gpwv,epwv)).weights

        end

    end

    gve = @view gve[:,:,1:ii]; gve = dropdims(sum(gve,dims=3),dims=3)
    @save "$(datadir("compiled/erady/$(gstn).jld2"))" gve

end
