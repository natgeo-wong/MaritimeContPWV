using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using ClimateSatellite
using Dates
using GeoRegions
using Logging
using NCDatasets

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

    for hr = 1 : 24

        pplt.close(); proj = pplt.Proj("eqc");
        f,axs = pplt.subplots(proj=proj,nrows=1,axwidth=5,aspect=2)

        c = axs[1].contourf(elon,elat,csf[:,:,hr]',cmap="Blues",levels=0.25:0.05:0.75)
        axs[1].contour(rlon,rlat,prcp[:,:,hr]',levels=[0,0.5],linewidth=0.5,color="r")
        axs[1].format(lonlim=(90,165),latlim=(-15,20),coast=true)
        axs[1].colorbar(c,loc="b",title="Column Relative Humidity")

        if !isdir(plotsdir("csfanim")); mkpath(plotsdir("csfanim")) end
        f.savefig(plotsdir(
            "csfanim/csfprcp-$(@sprintf("%02d",hr)).png"),
            transparent=false,dpi=200
        )

    end

end

droot = "/n/kuangdss01/lab/"
init,eroot = erastartup(aID=2,dID=1,path=droot);
elon,elat,rlon,rlat,csf,prcp = retrievecsfprcp(
    init,eroot,droot,
    regID="SEA",timeID=[2001,2018]
)
plotcsfprcp(elon,elat,csf,rlon,rlat,prcp)
