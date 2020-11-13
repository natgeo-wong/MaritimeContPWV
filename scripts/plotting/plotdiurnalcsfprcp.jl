using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using ClimateSatellite
using Dates
using DelimitedFiles
using GeoRegions
using Logging
using NCDatasets
using Printf
using Statistics

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function retrievecsfprcp(
    init::Dict,eroot::Dict,sroot::AbstractString;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    gres::Real=0
)

    emod,epar,ereg,etime = erainitialize(
        init;
        modID="csfc",parID="csf",regID=regID,timeID=timeID,
        gres=gres
    )

    nlon,nlat = ereg["size"]; elon = ereg["lon"]; elat = ereg["lat"]

    @info "$(Dates.now()) - Extracting relevant closest-coordinate points of GPM precipitation for each of the ERA5 total column water grid points ..."

    lon,lat = gpmlonlat(); rlon,rlat,_ = gregiongridvec(regID,lon,lat);
    nrlon = length(rlon); nrlat = length(rlat)

    csf = zeros(nlon,nlat,24)
    gpm = zeros(nrlon,nrlat,48)

    for yr = etime["Begin"] : etime["End"]

        @info "$(Dates.now()) - Extracting GPM Precipitation Rate and ERA5 Column Saturation Fraction data in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(yr) ..."

        cds,cvar = eraanaread("domain_yearly_mean_hourly",emod,epar,ereg,eroot,Date(yr));
        csf += cvar[:]*1; close(cds)

        rds,rvar = clisatanaread(
            "gpmimerg","prcp_rate","domain_yearly_mean_hourly",
            Date(yr),regID,path=sroot
        );
        gpm += rvar[:]*3600; close(rds);

    end

    nt = etime["End"] - etime["Begin"] + 1
    gpm = dropdims(mean(reshape(gpm,nrlon,nrlat,2,:),dims=3),dims=3)

    return elon,elat,rlon,rlat,csf/nt,gpm/nt

end

function plotcsfprcp(
    elon::Vector{<:Real}, elat::Vector{<:Real}, csf::Array{<:Real,3},
    rlon::Vector{<:Real}, rlat::Vector{<:Real}, prcp::Array{<:Real,3}
)

    coord = readdlm(plotsdir("SEA.txt"),comments=true,comment_char='#')
    x = coord[:,1]; y = coord[:,2];
    mcsf = dropdims(mean(csf,dims=3),dims=3)

    for hr = 1 : 24

        pplt.close(); f,axs = pplt.subplots(nrows=1,axwidth=5,aspect=15/7)

        c = axs[1].contourf(
            elon,elat,(csf[:,:,hr].-mcsf)',
            cmap="drywet",levels=vcat(-5:-1:-1,1:5),extend="both"
        )
        axs[1].plot(x,y,c="k",lw=0.5)
        axs[1].contour(rlon,rlat,prcp[:,:,hr]',levels=[0,0.5],linewidth=0.5,color="r")
        axs[1].format(xlim=(90,165),ylim=(-15,20),coast=true,xlocator=90:15:165)
        axs[1].colorbar(c,loc="b",title="Column Relative Humidity")

        if !isdir(plotsdir("csfanim")); mkpath(plotsdir("csfanim")) end
        f.savefig(plotsdir(
            "csfanim/csfprcp-$(@sprintf("%02d",hr)).png"),
            transparent=false,dpi=200
        )

    end

    rcsf = dropdims(maximum(csf,dims=3),dims=3) .- dropdims(minimum(csf,dims=3),dims=3)

    pplt.close(); f,axs = pplt.subplots(nrows=1,axwidth=5,aspect=15/7)

    c = axs[1].contourf(
        elon,elat,rcsf',
        cmap="Blues",levels=0:10,extend="max"
    )
    axs[1].plot(x,y,c="k",lw=0.5)
    axs[1].format(xlim=(90,165),ylim=(-15,20),coast=true,xlocator=90:15:165)
    axs[1].colorbar(c,loc="b",title="Column Relative Humidity")

    if !isdir(plotsdir("csfanim")); mkpath(plotsdir("csfanim")) end
    f.savefig(plotsdir(
        "csfanim/csf-range.png"),
        transparent=false,dpi=200
    )

    rprcp = dropdims(maximum(prcp,dims=3),dims=3) .- dropdims(minimum(prcp,dims=3),dims=3)

    pplt.close(); f,axs = pplt.subplots(nrows=1,axwidth=5,aspect=15/7)

    c = axs[1].contourf(
        rlon,rlat,rprcp',
        cmap="Blues",levels=0:10,extend="max"
    )
    axs[1].plot(x,y,c="k",lw=0.5)
    axs[1].format(xlim=(90,165),ylim=(-15,20),coast=true,xlocator=90:15:165)
    axs[1].colorbar(c,loc="b",title="Column Relative Humidity")

    if !isdir(plotsdir("csfanim")); mkpath(plotsdir("csfanim")) end
    f.savefig(plotsdir(
        "csfanim/prcp-range.png"),
        transparent=false,dpi=200
    )

end

droot = "/n/kuangdss01/lab/"
init,eroot = erastartup(aID=2,dID=1,path=droot);
elon,elat,rlon,rlat,csf,prcp = retrievecsfprcp(
    init,eroot,droot,
    regID="SEA",timeID=[2001,2018]
)
plotcsfprcp(elon,elat,csf,rlon,rlat,prcp)
