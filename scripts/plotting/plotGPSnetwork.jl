using DrWatson
using DelimitedFiles
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

ods = NCDataset(datadir("era5-GLBx0.25-oro.nc"));
mds = NCDataset(datadir("era5-GLBx0.25-lsm.nc"));
oro_e5 = ods["z"][:,:,1]/9.81/1000;
lsm_e5 = mds["lsm"][:,:,1]*1
oro_e5[oro_e5.<0] .= 0
lon_e5 = ods["longitude"][:]*1
lat_e5 = ods["latitude"][:]*1
close(ods); close(mds)

ods = NCDataset(datadir("erai-GLBx0.25-oro.nc"));
mds = NCDataset(datadir("erai-GLBx0.25-lsm.nc"));
oro_ei = ods["z"][:,:,1]/9.81/1000;
lsm_ei = mds["lsm"][:,:,1]*1
oro_ei[oro_ei.<0] .= 0
lon_ei = ods["longitude"][:]*1
lat_ei = ods["latitude"][:]*1
close(ods); close(mds)

ginfo = retrieveginfo()
glon = ginfo[:,2]; glat = ginfo[:,3]

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); f,axs = pplt.subplots(
    ncols=2,axwidth=3
)

axs[1].contourf(
    lon_e5,lat_e5,oro_e5',cmap="speed_r",levels=0:0.2:2
)
axs[1].plot(x,y,c="w",lw=0.5)
axs[1].scatter(glon,glat,s=5,c="r")
axs[1].format(xlim=(95,107),ylim=(-6,6),abc=true,rtitle="ERA5")

c = axs[2].contourf(
    lon_ei,lat_ei,oro_ei',cmap="speed_r",levels=0:0.2:2
)
axs[2].plot(x,y,c="w",lw=0.5)
axs[2].scatter(glon,glat,s=5,c="r")
axs[2].format(xlim=(95,107),ylim=(-6,6),abc=true,rtitle="ERA-Interim")

f.colorbar(c,loc="r",title="Surface Height / km")
f.savefig(plotsdir("oro.png"),transparent=false,dpi=200)
