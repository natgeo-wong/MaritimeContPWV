using ClimateERA
using ClimateSatellite
using GeoRegions
using JLD2
using StatsBase

include(srcdir("gnss.jl"))

function mimicvgnss(
    gstn::AbstractString,
    init::AbstractDict, eroot::AbstractDict, sroot::AbstractString;
    regID::AbstractString="GLB",
    gpwv::AbstractRange=0:100, mpwv::AbstractRange=0:100
)

    if init["datasetID"] == 1
        pmod,ppar,preg,_ = erainitialize(
            init,
            modID="csfc",parID="Pi_RE5",regID=regID,
            gres=1
        )
    else
        error("$(Dates.now()) - We use ERA5 Pi to calculate GPS PWV")
    end

    ginf = gstationinfo(gstn,retrieveginfo());
    glon = ginf["longitude"]; glat = ginf["latitude"]

    mlon,mlat = gpmlonlat(); mlon,mlat,_ = gregiongridvec(regID,mlon,mlat);

    mlon,mlat = regionpoint(glon,glat,mlon,mlat)
    plon,plat = regionpoint(glon,glat,preg["lon"],preg["lat"])

    ngps = length(gpwv); nera = length(epwv); nt = (etime["End"] - etime["Begin"] + 1) * 12
    gvm  = zeros(Int64,ngps-1,nera-1,nt); ii = 0; nhr = hrstep(cmod)

    @info "$(Dates.now()) - Comparing GNSS ZWD vs MIMIC-TPW2m at Station $gstn"

    for yr = etime["Begin"] : etime["End"], mo = 1 : 12

        date = Date(yr,mo)
        gfol = datadir("gnss/$(gstn)/$(yr2str(date))")
        gfnc = joinpath("$gfol","$(gstn)-$(yrmo2str(date)).nc")

        if isfile(gfnc); ii += 1

            gds = Dataset(gfnc); gzwd = gds["zwd"][:]*1; close(gds);
            gzwd = reshape(gzwd,6*nhr,:);
            gzwd = dropdims(mean(gzwd,dims=1),dims=1);

            pds,pvar = erarawread(pmod,ppar,preg,eroot,date)
            ipi  = pvar[plon,plat,:]*1; close(pds)
            gtcw = gzwd .* ipi * 1000

            mds,mvar = clisatrawread("mtpw2m","tpw",dtii,regID,path=sroot);
            mtpw = mvar[mlon,mlat,:]*1; close(mds)

            irem = (.!ismissing.(gtcw)) .& (.!isnan.(gtcw)) .& (.!ismissing.(mtpw)) .& (.!isnan.(mtpw));
            etcw = etcw[irem]; gtcw = gtcw[irem]
            gvm[:,:,ii] = fit(Histogram,(gtcw,mtpw),(gpwv,mpwv)).weights

        end

    end

    gvm = @view gvm[:,:,1:ii]; gvm = dropdims(sum(gvm,dims=3),dims=3)
    @save "$(datadir("compiled/$(init["prefix"])/$(gstn).jld2"))" gve

end
