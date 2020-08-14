using DrWatson
using DelimitedFiles
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

ds  = NCDataset(datadir("era5-SEAx0.25-lsm.nc"));
lsm = ds["lsm"][:,:,1]*1
lon = ds["longitude"][:]*1
lat = ds["latitude"][:]*1
close(ds)

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); proj = pplt.Proj("eqc"); f,axs = pplt.subplots(
    nrows=1,proj=proj,axwidth=10
)

axs[1].contourf(
    lon,lat,lsm',cmap="drywet_r",colorbar="r",
    norm="segmented",levels=[0,0.5,1],
)
axs[1].plot(x,y,c="k",lw=0.5)
axs[1].format(lonlim=(90,165),latlim=(-15,20))

f.savefig(plotsdir("lsm.png"),transparent=false,dpi=200)
