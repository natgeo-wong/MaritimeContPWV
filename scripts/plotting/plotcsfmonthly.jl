using DrWatson
@quickactivate "MaritimeContPWV"
using Dates
using DelimitedFiles
using ClimateERA
using JLD2

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plotcsfmonthly(
    init::Dict,eroot::Dict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    gres::Real=0
)

    emod,epar,ereg,etime = erainitialize(
        init;
        modID="csfc",parID="csf",regID=regID,timeID=timeID,
        gres=gres
    )
    nlon,nlat = ereg["size"]; elon = ereg["lon"]; elat = ereg["lat"]
    nt = etime["End"] - etime["Begin"] + 1

    csfavg = zeros(nlon,nlat,12)
    csfdhr = zeros(nlon,nlat,24,12)

    for yr = etime["Begin"] : etime["End"]
        cds,cvar = eraanaread("domain_monthly_mean_climatology",emod,epar,ereg,eroot,Date(yr));
        csfavg += cvar[:]*1; close(cds)
        cds,cvar = eraanaread("domain_monthly_mean_hourly",emod,epar,ereg,eroot,Date(yr));
        csfdhr += cvar[:]*1; close(cds)
    end

    csfdhr = dropdims(maximum(csfdhr,dims=3) .- minimum(csfdhr,dims=3),dims=3)
    csfavg = csfavg / nt
    csfdhr = csfdhr / nt

    coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
    x = coord[:,1]; y = coord[:,2];

    for mo = 1 : 12

        pplt.close(); f,axs = pplt.subplots(ncols=2,axwidth=3,aspect=15/7)

        c = axs[1].contourf(elon,elat,csfavg[:,:,mo]',cmap="Blues",levels=25:5:75,extend="both")
        axs[1].plot(x,y,c="k",lw=0.5)
        axs[1].format(rtitle="Average",suptitle="Column Saturation Fraction ($(monthname(mo))) / %")
        axs[1].colorbar(c,loc="r")

        c = axs[2].contourf(elon,elat,csfdhr[:,:,mo]',cmap="Blues",levels=0:10,extend="max")
        axs[2].plot(x,y,c="k",lw=0.5)
        axs[2].format(rtitle="Diurnal Variability")
        axs[2].colorbar(c,loc="r")

        for ax in axs
            ax.format(xlim=(90,165),ylim=(-15,20),xlocator=(6:11)*15)
        end

        f.savefig(plotsdir("csf-month-$mo.png"),transparent=false,dpi=200)

    end

end

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/");
plotcsfmonthly(init,eroot,regID="SEA")
