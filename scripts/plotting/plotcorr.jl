using DrWatson
@quickactivate "MaritimeContPWV"
using NCDatasets

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

ds  = NCDataset(datadir("compiled/eravmimic-SEAx0.25-corr.nc"));
rho = ds["rho"][:]
lon = ds["longitude"][:]
lat = ds["latitude"][:]

close(ds)
coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); f,axs = pplt.subplots(
    nrows=1,aspect=15/7,axwidth=6
)

axs[1].contourf(
    lon,lat,rho',cmap="viridis",colorbar="r",levels=0:0.1:1,
)
axs[1].plot(x,y,c="k",lw=0.5)
axs[1].format(xlim=(90,165),ylim=(-15,20),suptitle="ERA5 vs MIMIC Correlation")

f.savefig(plotsdir("eravmimic-SEAx0.25.png"),transparent=false,dpi=200)

ds = NCDataset(datadir("compiled/eravmimic-SEAx0.25.nc"))
freqall = ds["bin_frq"][:]; pwv = ds["pwv"][:]; close(ds)

mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
thr = 0.5; lsmsea = lsm .< thr; lsmlnd = lsm .> thr;

freqsea = dropdims(sum(freqall .* lsmsea,dims=(1,2)),dims=(1,2));
freqlnd = dropdims(sum(freqall .* lsmlnd,dims=(1,2)),dims=(1,2));
freqall = dropdims(sum(freqall,dims=(1,2)),dims=(1,2));

freqall = freqall ./ sum(freqall); freqall[freqall.==0] .= NaN;
freqlnd = freqlnd ./ sum(freqlnd); freqlnd[freqlnd.==0] .= NaN;
freqsea = freqsea ./ sum(freqsea); freqsea[freqsea.==0] .= NaN;

pplt.close(); f,axs = pplt.subplots(ncols=3,axwidth=2);

axs[1].pcolormesh(
    0:100,0:100,freqall' .* 100^2,cmap="gnbu",
    #levels=0:0.1:1,extend="max"
    norm="segmented",extend="max",
    levels=[0,1,10,20,50,90,100],
)
axs[1].plot([0,100],[0,100],c="k",lw=0.2)
axs[1].format(
    xlim=(0,100),xlabel="MIMIC-TPW2m Precipitable Water / mm",
    ylim=(0,100),ylabel="Reanalysis Precipitable Water / mm",
    suptitle="All GNSS observations",rtitle="All",
    abc=true,grid="on"
)

axs[2].pcolormesh(
    0:100,0:100,freqlnd' .* 100^2,cmap="gnbu",
    #levels=0:0.1:1,extend="max"
    norm="segmented",extend="max",
    levels=[0,1,10,20,50,90,100],
)
axs[2].plot([0,100],[0,100],c="k",lw=0.2)
axs[2].format(
    xlim=(0,100),xlabel="MIMIC-TPW2m Precipitable Water / mm",
    ylim=(0,100),ylabel="Reanalysis Precipitable Water / mm",
    suptitle="All GNSS observations",rtitle="Land",
    abc=true,grid="on"
)

c = axs[3].pcolormesh(
    0:100,0:100,freqsea' .* 100^2,cmap="gnbu",
    #levels=0:0.1:1,extend="max"
    norm="segmented",extend="max",
    levels=[0,1,10,20,50,90,100],
)
axs[3].plot([0,100],[0,100],c="k",lw=0.2)
axs[3].format(
    xlim=(0,100),xlabel="MIMIC-TPW2m Precipitable Water / mm",
    ylim=(0,100),ylabel="Reanalysis Precipitable Water / mm",
    suptitle="All MIMIC-TPW2m observations",rtitle="Sea",
    abc=true,grid="on"
)

f.colorbar(c,loc="r",label="Density")
f.savefig(plotsdir("eravmimic.png"),transparent=false,dpi=200)
