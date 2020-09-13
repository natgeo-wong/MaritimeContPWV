using DrWatson
using DelimitedFiles
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

cds = NCDataset(datadir("compiled/erac5-SEAx0.25-csf-sfc.nc"));
lon = cds["longitude"][:]*1
lat = cds["latitude"][:]*1
csf = cds["average"]*1
#close(cds);

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); proj = pplt.Proj("eqc"); f,axs = pplt.subplots(proj=proj,nrows=1,axwidth=5)

c = axs[1].contourf(lon,lat,csf',cmap="speed_r",levels=0.25:0.05:0.75)
axs[1].plot(x,y,c="k",lw=0.5)
axs[1].format(lonlim=(90,165),latlim=(-15,20))

f.colorbar(c,loc="r",title="Surface Height / km")
f.savefig(plotsdir("csf.png"),transparent=false,dpi=200)
