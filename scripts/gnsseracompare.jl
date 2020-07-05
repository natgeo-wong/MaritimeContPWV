using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateERA
using StatsBase
using JLD2

function eravgnss(
    gstn::AbstractString,
    init::AbstractDict, eroot::AbstractDict;
    gZWD::AbstractRange=0:0.01:0.6, ePWV::AbstractRange=0:100
)

    emod,epar,ereg,etime = erainitialize(init,modID="msfc",parID="tcw",regID="SMT")
    ginf = gstationinfo(gstn,retrieveginfo()); glon = ginfo[:,2]; glat = ginfo[:,3]
    ilon,ilat = regionpoint(glon,glat,ereg["lon"],ereg["lat"])
    nzwd = length(gZWD); npwv = length(ePWV); nt = (etime["End"] - etime["Begin"] + 1) * 12
    gve  = zeros(Int64,nzwd,npwv,nt); ii = 0; nhr = hrstep(emod)

    @info "$(Dates.now()) -  Comparing GNSS ZWD vs ECMWF reanalysis TCW at Station $gstn"

    for yr = etime["Begin"] : etime["End"], mo = 1 : 12

        date = Date(yr,mo)
        gfol = datadir("gnss/gstn/$(yr2str(date))")
        gfnc = joinpath("$fol","$(gstn)-$(yrmo2str(date)).nc")

        if isfile(gfnc); ii = 1

            gds = Dataset(gfnc); gzwd = gds["zwd"][:]*1; close(gds);
            gzwd = reshape(gzwd,6*nhr,:);
            gzwd = dropdims(mean(gzwd,dims=1),dims=1)

            eds,evar = erarawread(emod,epar,ereg,eroot,date)
            etcw = evar[ilon,ilat,:]*1; close(eds)

            irem = (gzwd .!= missing) .& (gzwd .!= NaN);
            etcw = etcw[irem]; gzwd = gzwd[irem]
            gve[:,:,ii] = fit(Histogram,(gzwd,etcw),(gZWD,ePWV)).weights

        end

    end

    gve = @view gve[:,:,1:ii]; gve = dropdims(sum(gve,3),dims=3)
    @save "$(datadir("compiled/$(gstn).jld2"))" gve

end

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/")
mkpath(datadir("compiled"))

gstns = retrieveginfo()[:,1]; eravgnss.(gstns,init,eroot)
