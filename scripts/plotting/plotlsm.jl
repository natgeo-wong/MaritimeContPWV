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

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); proj = pplt.Proj("eqc"); f,axs = pplt.subplots(proj=proj,nrows=2,axwidth=5)

axs[1].contourf(lon_e5,lat_e5,oro_e5',cmap="speed_r",levels=0:0.2:2)
axs[1].plot(x,y,c="k",lw=0.5)
axs[1].format(lonlim=(90,165),latlim=(-15,20),abc=true,rtitle="ERA5")

c = axs[2].contourf(lon_ei,lat_ei,oro_ei',cmap="speed_r",levels=0:0.2:2)
axs[2].plot(x,y,c="k",lw=0.5)
axs[2].format(lonlim=(90,165),latlim=(-15,20),abc=true,rtitle="ERA-Interim")

f.colorbar(c,loc="r",title="Surface Height / km")
f.savefig(plotsdir("lsm.png"),transparent=false,dpi=200)
