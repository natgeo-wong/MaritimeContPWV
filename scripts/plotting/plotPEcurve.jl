using DrWatson
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plotPEcurvegeneral(dID::AbstractString; fthr::Integer=25)

    @load "$(datadir("compiled/tcwVprcp/$(dID).jld2"))" prcp freq
    if uppercase(dID) == "ERA5"
          prcp = prcp * 1000 # Units is in m (hourly, so no adjust for time)
    else; prcp = prcp * 3600 # Units is in mm s-1
    end

    freq = dropdims(sum(freq,dims=(1,2)),dims=(1,2))
    prcp = dropdims(sum(prcp,dims=(1,2)),dims=(1,2)) ./ freq
    prcp[freq.<25] .= NaN

    axs[1].plot(10:0.5:90,prcp,lw=1,label="$(uppercase(dID))",legend="ul")

end

function plotPEcurvelandsea(
      dID::AbstractString, jj::Integer;
      lthr::Real=0.75, othr::Real=0.25, fthr::Integer=25
)

    @load "$(datadir("compiled/tcwVprcp/$(dID).jld2"))" prcp freq
    if uppercase(dID) == "ERA5"
          prcp = prcp * 1000 # Units is in m (hourly, so no adjust for time)
    else; prcp = prcp * 3600 # Units is in mm s-1
    end

    freqall = dropdims(sum(freq,dims=(1,2)),dims=(1,2))
    prcpall = dropdims(sum(prcp,dims=(1,2)),dims=(1,2)) ./ freqall
    prcpall[freqall.<fthr] .= NaN

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
    lsmsea = lsm .< othr; lsmlnd = lsm .> lthr; lsmcst = (lsm .<= lthr) .& (lsm .>= othr)

    freqsea = freq .* lsmsea; prcpsea = prcp .* lsmsea
    freqlnd = freq .* lsmlnd; prcplnd = prcp .* lsmlnd
    freqcst = freq .* lsmcst; prcpcst = prcp .* lsmcst

    freqsea = dropdims(sum(freqsea,dims=(1,2)),dims=(1,2))
    prcpsea = dropdims(sum(prcpsea,dims=(1,2)),dims=(1,2)) ./ freqsea
    prcpsea[freqsea.<fthr] .= NaN

    freqlnd = dropdims(sum(freqlnd,dims=(1,2)),dims=(1,2))
    prcplnd = dropdims(sum(prcplnd,dims=(1,2)),dims=(1,2)) ./ freqlnd
    prcplnd[freqlnd.<fthr] .= NaN

    freqcst = dropdims(sum(freqcst,dims=(1,2)),dims=(1,2)); freqcst[freqcst.<25] .= 0
    prcpcst = dropdims(sum(prcpcst,dims=(1,2)),dims=(1,2)) ./ freqcst
    prcpcst[freqcst.<fthr] .= NaN

    axs[jj].plot(10:0.5:90,prcpall,lw=1,label="All",legend="ul",c="k")
    axs[jj].plot(10:0.5:90,prcpsea,lw=1,label="Sea",legend="ul",c="b")
    axs[jj].plot(10:0.5:90,prcplnd,lw=1,label="Land",legend="ul",c="g")
    # axs[jj].plot(10:0.5:90,prcpcst,lw=1,label="Coast",legend="ul",c="r")

    axs[jj].format(
        xlim=(10,90),xlabel="Precipitable Water / mm",
        ylim=(1e-4,100),ylabel=L"Precipitation Rate / mm hr$^{-1}$",yscale="log",
        title="Precipitation Data Source: $(uppercase(dID))",abc=true
    )

end

arr = [[0,1,1,0],[2,2,3,3]]
pplt.close(); f,axs = pplt.subplots(arr,aspect=2,axwidth=3,sharex=3,sharey=3);

plotPEcurvegeneral("gpm");
plotPEcurvegeneral("era5");
axs[1].format(
    xlim=(10,90),
    ylim=(1e-4,100),ylabel=L"Precipitation Rate / mm hr$^{-1}$",yscale="log",
    title="Summary",abc=true,
    suptitle="P-E curve (2001-2018)"
)

plotPEcurvelandsea("gpm",2,othr=0.5,lthr=0.5);
plotPEcurvelandsea("era5",3,othr=0.5,lthr=0.5);

f.savefig(plotsdir("PEcurve.png"),transparent=false,dpi=200)
