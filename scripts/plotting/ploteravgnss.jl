using DrWatson
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

include(srcdir("gnss.jl"))

function ploteravsgnss(stn::AbstractString,dset::AbstractString)

    pplt.close(); f,axs = pplt.subplots(axwidth=3)
    @load "$(datadir("compiled/$dset/$(stn).jld2"))" gve
    c = axs[1].contourf(
        0.5:99.5,0.5:99.5,gve' ./ sum(gve) .* 100,
        norm="segmented",levels=[0,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5],
    )
    axs[1].plot([0,100],[0,100],c="k",lw=0.5)
    axs[1].format(
        xlim=(0,100),xlabel="GNSS-derived Precipitable Water / mm",
        ylim=(0,100),ylabel="$(uppercase(dset))-retrieved Precipitable Water / mm",
        suptitle="$(stn)"
    )
    f.colorbar(c,loc="r")
    f.savefig(plotsdir("eravsgnss-$dset/$(stn).png"),transparent=false,dpi=200)

end

function ploteravsgnssall(stn::AbstractString,dset::AbstractString,jj::Integer)

    @load "$(datadir("compiled/$(stn).jld2"))" gve
    axs[jj].contourf(
        0.5:99.5,0.5:99.5,gve' ./ sum(gve) .* 100,
        levels=0:0.25:2.5
        #norm="segmented",levels=[0,0.1,0.2,0.5,1,2,5,10,20,50,100],
    )
    axs[jj].format(xlim=(0,100),ylim=(0,100))

end

function ploteravsgnsscompiled(dset::AbstractString,jj::Integer;title::AbstractString)

    gstns = retrieveginfo()[:,1];
    gvec = zeros(100,100)

    for gstn in gstns
        @load "$(datadir("compiled/$dset/$(gstn).jld2"))" gve
        gvec = gvec + gve
    end

    sgvec = sum(gvec); gvec[gvec.==0] .= NaN
    nbins = prod(size(gvec))

    c = axs[jj].pcolormesh(
        0:100,0:100,gvec' ./ sgvec .* nbins,cmap="gnbu",
        #levels=0:0.1:1,extend="max"
        norm="segmented",extend="both",
        levels=15:15:150,
        cmap_kw=Dict("left"=>0.02)
    )
    axs[jj].plot([0,100],[0,100],c="k",lw=0.2)
    axs[jj].format(
        xlim=(0,100),xlabel="GNSS-derived Precipitable Water / mm",
        ylim=(0,100),ylabel="Reanalysis Precipitable Water / mm",
        suptitle="All GNSS observations",
        ultitle="$title\nNo. of Obs.: $(@sprintf("%d",sum(sgvec)))",grid="on",
        abc=true
    )

    return c

end

function ploteravsgnssnormal(dset::AbstractString,jj::Integer;title::AbstractString)

    gstns = retrieveginfo()[:,1];
    gvec = zeros(100,100)

    for gstn in gstns
        @load "$(datadir("compiled/$dset/$(gstn).jld2"))" gve
        gvec = gvec + gve
    end

    mgvec = maximum(gvec); sgvec = sum(gvec); gvec[gvec.==0] .= NaN

    c = axs[jj].pcolormesh(
        0:100,0:100,gvec' ./ mgvec,cmap="gnbu",
        levels=[0,0.01,0.1,0.2,0.5,0.9,1],
        cmap_kw=Dict("left"=>0.1)
    )
    axs[jj].plot([0,100],[0,100],c="k",lw=0.2)
    axs[jj].format(
        xlim=(0,100),xlabel="GPS-derived Precipitable Water / mm",
        ylim=(0,100),ylabel="Reanalysis Precipitable Water / mm",
        suptitle="All GNSS observations",
        ultitle="$title\nNo. of Obs.: $(@sprintf("%d",sum(sgvec)))",grid="on",
        abc=true
    )

    return c

end

gstns = retrieveginfo()[:,1]; jj = 0; dset = "erai"
mkpath(plotsdir("eravsgnss-$dset"))

# pplt.close(); f,axs = pplt.subplots(nrows=8,ncols=8,axwidth=4)
#
# for gstn in gstns; global jj += 1
#     ploteravsgnssall(gstn,dset,jj)
# end
#
# axs[1].contourf(
#     (1:2)/10,1:2,[[NaN,NaN],[NaN,NaN]],
#     levels=0:0.25:2.5
#     #norm="segmented",levels=[0,0.1,0.2,0.5,1,2,5,10,20,50,100],
# )
#
# f.savefig(plotsdir("eravgnss-$dset.png"),transparent=false,dpi=200)

# for gstn in gstns
#     ploteravsgnss(gstn,dset)
# end

pplt.close(); f,axs = pplt.subplots(ncols=3,axwidth=2)
c = ploteravsgnsscompiled("era5",1,title="ERA5 (Hourly)")
# c = ploteravsgnsscompiled("era6",2,title="ERA5 (6-Hourly)")
c = ploteravsgnsscompiled("erady",2,title="ERA5 (Daily)")
f.colorbar(c,loc="r",label="Density")
f.savefig(plotsdir("eravgnsscompiled2.png"),transparent=false,dpi=200)

pplt.close(); f,axs = pplt.subplots(ncols=2,axwidth=2)
c = ploteravsgnssnormal("era5",1,title="ERA5 (Hourly)")
c = ploteravsgnssnormal("erady",2,title="ERA5 (Daily)")
f.colorbar(c,loc="r",label="Normalized Frequency")
f.savefig(plotsdir("eravgnssnormal2.png"),transparent=false,dpi=200)
