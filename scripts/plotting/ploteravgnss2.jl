using DrWatson
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

include(srcdir("gnss.jl"))

function ploteravsgnsscompiled(dset::AbstractString,ax::PyObject;title::AbstractString)

    gstns = retrieveginfo()[:,1];
    gvec = zeros(100,100)

    for gstn in gstns
        @load "$(datadir("compiled/$dset/$(gstn).jld2"))" gve
        gvec = gvec + gve
    end

    sgvec = sum(gvec); gvec[gvec.==0] .= NaN
    nbins = prod(size(gvec))

    c = ax.pcolormesh(
        0:100,0:100,gvec' ./ sgvec .* nbins,cmap="gnbu",
        #levels=0:0.1:1,extend="max"
        norm="segmented",extend="both",
        levels=(1:10)*10,
        cmap_kw=Dict("left"=>0.02)
    )
    ax.plot([0,100],[0,100],c="k",lw=0.2)
    ax.format(
        ultitle="$title\nNo. of Obs.: $(@sprintf("%d",sum(sgvec)))",grid="on",
        abc=true
    )

    return c

end

pplt.close(); f,axs = pplt.subplots(ncols=3,axwidth=2,sharey=0,sharex=0,wspace=[1.2,0])


ods = NCDataset(datadir("era5-GLBx0.25-oro.nc"));
mds = NCDataset(datadir("era5-GLBx0.25-lsm.nc"));
oro_e5 = ods["z"][:,:,1]/9.81/1000;
lsm_e5 = mds["lsm"][:,:,1]*1
oro_e5[oro_e5.<0] .= 0
lon_e5 = ods["longitude"][:]*1
lat_e5 = ods["latitude"][:]*1
close(ods); close(mds)

ginfo = retrieveginfo()
glon = ginfo[:,2]; glat = ginfo[:,3]

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

c = axs[1].contourf(lon_e5,lat_e5,oro_e5',cmap="speed_r",levels=0.2:0.2:2,extend="both",cmap_kw=Dict("left"=>0.2))
axs[1].plot(x,y,c="w",lw=0.5)
axs[1].scatter(glon,glat,s=5,c="r")
axs[1].format(
    xlim=(95,107),ylim=(-6,6),xlocator=95:2:107,
    xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
    abc=true
)
axs[1].colorbar(c,loc="r",title="Surface Height / km")

c = ploteravsgnsscompiled("era5",axs[2],title="Hourly-Averaged")
axs[2].format(
    xlim=(0,100),xlabel="SuGAr Precipitable Water / mm",
    ylim=(0,100),ylabel="ERA5 Precipitable Water / mm",
    abc=true
)
c = ploteravsgnsscompiled("erady",axs[3],title="Daily-Averaged")
axs[3].format(
    xlim=(0,100),xlabel="SuGAr Precipitable Water / mm",
    ylim=(0,100),ylocator=0:20:100,yticklabels=[],
    abc=true
)
f.colorbar(c,loc="r",label="Density")
f.savefig(plotsdir("eravgnsscompiled2.png"),transparent=false,dpi=200)
