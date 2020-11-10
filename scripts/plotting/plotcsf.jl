using DrWatson
using DelimitedFiles
@quickactivate "MaritimeContPWV"
using ClimateERA
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/")
emod,epar,ereg,etime = erainitialize(init;modID="csfc",parID="csf",regID="SEA")
cds,cvar = eracmpread("average",emod,epar,ereg,eroot); avg = cvar[:]*1; close(cds)
cds,cvar = eracmpread("variability_diurnal",emod,epar,ereg,eroot)
dhr = cvar[:]*1
lon = cds["longitude"][:]
lat = cds["latitude"][:]
close(cds)

coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
x = coord[:,1]; y = coord[:,2];

pplt.close(); f,axs = pplt.subplots(nrows=2,axwidth=5,aspect=15/7)

c = axs[1].contourf(lon,lat,avg',cmap="Blues",levels=0:10:100)
axs[1].plot(x,y,c="k",lw=0.5)
axs[1].format(rtitle="Average",suptitle="Column Saturation Fraction / %")
axs[1].colorbar(c,loc="r")

c = axs[2].contourf(lon,lat,dhr',cmap="Blues",levels=0:10)
axs[2].plot(x,y,c="k",lw=0.5)
axs[2].format(rtitle="Diurnal Variability")
axs[2].colorbar(c,loc="r")

for ax in axs
    ax.format(xlim=(90,165),ylim=(-15,20),xlocator=(6:11)*15)
end

f.savefig(plotsdir("csf.png"),transparent=false,dpi=200)
