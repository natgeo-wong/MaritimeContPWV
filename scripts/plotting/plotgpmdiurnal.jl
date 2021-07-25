using DrWatson
@quickactivate "MaritimeContPWV"

using JLD2

using ImageShow, FileIO
using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

@load datadir("test.jld2") rlon rlat gpm θmat
coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); f,axs = pplt.subplots(nrows=1,axwidth=2.5,aspect=1)

c = axs[1].contourf(rlon,rlat,θmat',cmap="romaO",levels=0:24)
axs[1].plot(x,y,c="k",lw=0.5)
axs[1].format(xlim=(95,107),ylim=(-6,6))
axs[1].colorbar(c,loc="r")

f.savefig(plotsdir("diurnalphase.png"),transparent=false,dpi=200)
load(plotsdir("diurnalphase.png"))
