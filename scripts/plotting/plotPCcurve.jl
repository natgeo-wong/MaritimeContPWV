using DrWatson
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plotPEcurvegeneral(dID::AbstractString, axs, month::Integer=0; density::Integer=25)

    if iszero(month)
          @load "$(datadir("compiled/csfVprcp/$(dID).jld2"))" prcp freq csf
    else; @load "$(datadir("compiled/csfVprcp/$(dID)-$month.jld2"))" prcp freq csf
    end

    if uppercase(dID) == "ERA5"
          prcp = prcp * 1000 # Units is in m (hourly, so no adjust for time)
    else; prcp = prcp * 3600 # Units is in mm s-1
    end

    freq = dropdims(sum(freq,dims=(1,2)),dims=(1,2))
    prcp = dropdims(sum(prcp,dims=(1,2)),dims=(1,2)) ./ freq
    prcp[freq.<25] .= NaN

    axs[1].plot(csf,prcp,lw=1,label="$(uppercase(dID))",legend="ul")

end

function plotPEcurvelandsea(
      dID::AbstractString, jj::Integer,axs, month::Integer=0;
      coast::Real=0.5, density::Real=0.05
)

    if iszero(month)
          @load "$(datadir("compiled/csfVprcp/$(dID).jld2"))" prcp freq csf
    else; @load "$(datadir("compiled/csfVprcp/$(dID)-$month.jld2"))" prcp freq csf
    end

    if uppercase(dID) == "ERA5"
          prcp = prcp * 1000 # Units is in m (hourly, so no adjust for time)
    else; prcp = prcp * 3600 # Units is in mm s-1
    end

    prcps = prcp ./ freq

    nbins = length(csf)
    freqall = dropdims(sum(freq,dims=(1,2)),dims=(1,2))
    prcpall = dropdims(sum(prcp,dims=(1,2)),dims=(1,2)) ./ freqall
    freqall = freqall / sum(freqall) * nbins
    prcpall[freqall.<density] .= NaN

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
    lsmsea = lsm .< coast; lsmlnd = lsm .> coast

    freqsea = freq .* lsmsea; prcpsea = prcp .* lsmsea
    freqlnd = freq .* lsmlnd; prcplnd = prcp .* lsmlnd

    lsmsea = dropdims(lsmsea,dims=3)
    lsmlnd = dropdims(lsmlnd,dims=3)

    prcpsea25 = zeros(nbins); prcpsea50 = zeros(nbins); prcpsea75 = zeros(nbins)
    prcplnd25 = zeros(nbins); prcplnd50 = zeros(nbins); prcplnd75 = zeros(nbins)

    for icsf = 1 : nbins

        prcpseaii = prcps[:,:,icsf]; prcpseaii = prcpseaii[lsmsea]
        prcpseaii = prcpseaii[.!isnan.(prcpseaii)]

        if !isempty(prcpseaii)
            prcpsea25[icsf] = percentile(prcpseaii,5)
            prcpsea75[icsf] = percentile(prcpseaii,95)
        end

        prcplndii = prcps[:,:,icsf]; prcplndii = prcplndii[lsmlnd]
        prcplndii = prcplndii[.!isnan.(prcplndii)]

        if !isempty(prcplndii)
            prcplnd25[icsf] = percentile(prcplndii,5)
            prcplnd75[icsf] = percentile(prcplndii,95)
        end

    end

    freqsea = dropdims(sum(freqsea,dims=(1,2)),dims=(1,2))
    prcpsea = dropdims(sum(prcpsea,dims=(1,2)),dims=(1,2)) ./ freqsea
    freqsea = freqsea / sum(freqsea) * nbins
    prcpsea[freqsea.<density] .= NaN
    prcpsea25[freqsea.<density] .= NaN
    prcpsea75[freqsea.<density] .= NaN

    freqlnd = dropdims(sum(freqlnd,dims=(1,2)),dims=(1,2))
    prcplnd = dropdims(sum(prcplnd,dims=(1,2)),dims=(1,2)) ./ freqlnd
    freqlnd = freqlnd / sum(freqlnd) * nbins
    prcplnd[freqlnd.<density] .= NaN
    prcplnd25[freqlnd.<density] .= NaN
    prcplnd75[freqlnd.<density] .= NaN

    axs[jj].plot(csf,prcpall,lw=1,label="All",legend="ul",c="k")

    axs[jj].plot(csf,prcpsea,fadedata=[prcpsea25,prcpsea75],lw=1,c="b")
    axs[jj].plot(csf,prcplnd,fadedata=[prcplnd25,prcplnd75],lw=1,c="g")

    axs[jj].plot(csf,prcpsea,lw=1,label="Sea",legend="ul",c="b")
    axs[jj].plot(csf,prcplnd,lw=1,label="Land",legend="ul",c="g")

    axs[jj].format(
        xlim=(0,nbins),xlabel="Column Saturation Fraction",
        ylim=(1e-4,30),ylabel=L"Precipitation Rate / mm hr$^{-1}$",yscale="log",
        title="Precipitation Data Source: $(uppercase(dID))",abc=true
    )

end


# pplt.close();
# sb = [[0,1,1,0],[2,2,3,3]]; f,axs = pplt.subplots(sb,aspect=2,axwidth=3,sharex=3,sharey=3)
#                             # f,axs = pplt.subplots(ncols=2,aspect=2,axwidth=4);
#
# plotPEcurvegeneral("gpm");
# plotPEcurvegeneral("era5");
# axs[1].format(
#     xlim=(0,120),
#     ylim=(1e-4,30),ylabel=L"Precipitation Rate / mm hr$^{-1}$",yscale="log",
#     title="Summary",abc=true,
#     suptitle="P-C curve"
# )
#
# plotPEcurvelandsea("gpm",2,density=0.05);
# plotPEcurvelandsea("era5",3,density=0.05);
#
# f.savefig(plotsdir("PCcurve.png"),transparent=false,dpi=200)

for mo in 1 : 12

    pplt.close();
    sb = [[0,1,1,0],[2,2,3,3]]; f,axs = pplt.subplots(sb,aspect=2,axwidth=3,sharex=3,sharey=3)
                                # f,axs = pplt.subplots(ncols=2,aspect=2,axwidth=4);

    plotPEcurvegeneral("gpm",axs,mo);
    plotPEcurvegeneral("era5",axs,mo);
    axs[1].format(
        xlim=(0,120),
        ylim=(1e-4,30),ylabel=L"Precipitation Rate / mm hr$^{-1}$",yscale="log",
        title="Summary",abc=true,
        suptitle="P-C curve ($(monthname(mo)))"
    )

    plotPEcurvelandsea("gpm",2,axs,mo,density=0.05);
    plotPEcurvelandsea("era5",3,axs,mo,density=0.05);

    f.savefig(plotsdir("PCcurve-$mo.png"),transparent=false,dpi=200)

end
