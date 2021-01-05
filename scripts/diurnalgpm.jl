using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateSatellite
using Dates
using GeoRegions
using Interpolations
using JLD2
using Logging
using LsqFit
using NCDatasets
using Statistics

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

mfit(t,p) = p[1] * cos.((t.-p[2])*pi/12)

function retrieveprcp(
    sroot::AbstractString;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0
)

    lon,lat = gpmlonlat(); rlon,rlat,_ = gregiongridvec(regID,lon,lat);
    nrlon = length(rlon); nrlat = length(rlat)

    gpm  = zeros(nrlon,nrlat,48)
    θmat = zeros(nrlon,nrlat)

    nt = timeID[2] - timeID[1] + 1
    for yr = timeID[1] : timeID[2]

        @info "$(now()) - Extracting GPM Precipitation Rate data in $regID domain for $(yr) ..."

        rds,rvar = clisatanaread(
            "gpmimerg","prcp_rate","domain_yearly_mean_hourly",
            Date(yr),regID,path=sroot
        );
        gpm += rvar[:]*3600; close(rds);

    end

    gpm .= gpm / nt
    vart = zeros(49)
    p0 = [0.5, 0.5]

    for ilat = 1 : nrlat, ilon = 1 : nrlon

        t = (0:48)/2 .+ rlon[ilon]/15
        vart[1:48] = gpm[ilon,ilat,:]
        vart[end]  = gpm[ilon,ilat,1]
        itp = interpolate(vart,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,t)
        etp = extrapolate(stp,Periodic())

        gpm[ilon,ilat,:] .= etp[0:0.5:23.5]

        fit = curve_fit(mfit, 0:0.5:23.5, (@view gpm[ilon,ilat,:]), p0)
        if fit.param[1] < 0
              θmat[ilon,ilat] = mod(fit.param[2]+24,24)
        else; θmat[ilon,ilat] = mod(fit.param[2],24)
        end

    end

    @save datadir("test.jld2") rlon rlat gpm θmat

end

droot = "/n/kuangdss01/lab/"
retrieveprcp(droot,regID="SEA",timeID=[2001,2018])
