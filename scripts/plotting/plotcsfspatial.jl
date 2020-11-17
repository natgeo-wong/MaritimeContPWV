using DrWatson
@quickactivate "MaritimeContPWV"
using NCDatasets

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plotPEcurvelandsea(
      dID::AbstractString, jj::Integer, axs; month::Integer=0,
      csfthr::Real, coast::Real=0.5, density::Real=0.05
)

    if iszero(month)
          @load "$(datadir("compiled/csfVprcp/$(dID).jld2"))" prcp freq csf
    else; @load "$(datadir("compiled/csfVprcp/$(dID)-$month.jld2"))" prcp freq csf
    end

    if uppercase(dID) == "ERA5"
          prcp = prcp * 1000 # Units is in m (hourly, so no adjust for time)
    else; prcp = prcp * 3600 # Units is in mm s-1
    end

    nbins = length(csf)
    iicsf = argmin(abs.(csf.-csfthr))
    prcps = prcp ./ freq

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc"))
    lon = mds["longitude"][:]
    lat = mds["latitude"][:]
    lsm = mds["lsm"][:]*1
    close(mds)

    lsmsea = lsm .< coast; freqsea = freq .* lsmsea; prcpsea = prcp .* lsmsea
    freqsea = dropdims(sum(freqsea,dims=(1,2)),dims=(1,2))
    prcpsea = dropdims(sum(prcpsea,dims=(1,2)),dims=(1,2)) ./ freqsea
    freqsea = freqsea / sum(freqsea) * nbins

    coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
    x = coord[:,1]; y = coord[:,2];

    axs[jj].contourf(
        lon,lat,prcps[:,:,iicsf]'/prcpsea[iicsf],
        norm="segmented",cmap="drywet",extend="both",
        levels=[0.05,0.1,0.2,0.5,0.67,1,1.5,2,5,10,20]
    )
    axs[jj].plot(x,y,c="k",lw=0.5)
    axs[jj].format(
        xlim=(90,165),ylim=(-15,20),urtitle="CSF: $csfthr %",
        xlocator=90:15:165
    )

end

pplt.close(); f,axs = pplt.subplots(nrows=5,ncols=2,aspect=15/7,axwidth=3)

c = axs[1].contourf(
    [1,2],[1,2],ones(2,2)*NaN,
    norm="segmented",cmap="drywet",extend="both",
    levels=[0.05,0.1,0.2,0.5,0.67,1,1.5,2,5,10,20]
)

plotPEcurvelandsea("gpm",1,axs,csfthr=40)
plotPEcurvelandsea("gpm",3,axs,csfthr=50)
plotPEcurvelandsea("gpm",5,axs,csfthr=60)
plotPEcurvelandsea("gpm",7,axs,csfthr=70)
plotPEcurvelandsea("gpm",9,axs,csfthr=85)

plotPEcurvelandsea("era5",2,axs,csfthr=40)
plotPEcurvelandsea("era5",4,axs,csfthr=50)
plotPEcurvelandsea("era5",6,axs,csfthr=60)
plotPEcurvelandsea("era5",8,axs,csfthr=70)
plotPEcurvelandsea("era5",10,axs,csfthr=85)

axs[1].format(title="GPM",suptitle="Rainfall Rate / Ratio against Domain Ocean Mean")
axs[2].format(title="ERA5")

f.colorbar(c,loc="r")
f.savefig(plotsdir("csfspatial-ratio.png"),transparent=false,dpi=200)

# for mo = 1 : 12
#
#     pplt.close(); f,axs = pplt.subplots(nrows=5,ncols=2,aspect=15/7,axwidth=3)
#
#     c = axs[1].contourf(
#         [1,2],[1,2],ones(2,2)*NaN,
#         norm="segmented",cmap="drywet",extend="both",
#         levels=[0.05,0.1,0.2,0.5,0.67,1,1.5,2,5,10,20]
#     )
#
#     plotPEcurvelandsea("gpm",1,axs,month=mo,csfthr=40)
#     plotPEcurvelandsea("gpm",3,axs,month=mo,csfthr=50)
#     plotPEcurvelandsea("gpm",5,axs,month=mo,csfthr=60)
#     plotPEcurvelandsea("gpm",7,axs,month=mo,csfthr=70)
#     plotPEcurvelandsea("gpm",9,axs,month=mo,csfthr=85)
#
#     plotPEcurvelandsea("era5",2,axs,month=mo,csfthr=40)
#     plotPEcurvelandsea("era5",4,axs,month=mo,csfthr=50)
#     plotPEcurvelandsea("era5",6,axs,month=mo,csfthr=60)
#     plotPEcurvelandsea("era5",8,axs,month=mo,csfthr=70)
#     plotPEcurvelandsea("era5",10,axs,month=mo,csfthr=85)
#
#     axs[1].format(title="GPM",suptitle="Rainfall Rate ($(monthname(mo))) / Ratio against Domain Ocean Mean")
#     axs[2].format(title="ERA5")
#
#     f.colorbar(c,loc="r")
#     f.savefig(plotsdir("csfspatial-ratio-$mo.png"),transparent=false,dpi=200)
#
# end
