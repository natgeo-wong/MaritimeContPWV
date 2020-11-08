using DrWatson
using DelimitedFiles
@quickactivate "MaritimeContPWV"
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

cds = NCDataset(datadir("era5-SEAx0.25-csf-sfc-197901.nc"));
lon = cds["longitude"][:]*1
lat = cds["latitude"][:]*1
csf = cds["csf"][:]*1
#close(cds);

csfm = dropdims(mean(csf,dims=3),dims=3); nt = size(csf,3)

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

for it = 1 : nt

    pplt.close(); f,axs = pplt.subplots(nrows=1,axwidth=5,aspect=15/7)

    c = axs[1].contourf(lon,lat,csf[:,:,it]',cmap="speed_r",levels=0:10:100)
    axs[1].plot(x,y,c="k",lw=0.5)
    axs[1].format(xlim=(90,165),ylim=(-15,20))

    f.colorbar(c,loc="r",title="Surface Height / km")
    f.savefig(plotsdir("test/csf-$it.png"),transparent=false,dpi=200)

end
