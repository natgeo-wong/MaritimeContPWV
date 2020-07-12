using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using Dates
using JLD2
using NCDatasets
using Statistics

function pecompile(
    dID::AbstractString;
    regID::AbstractString, timeID::Union{Integer,Vector}=0,
    hres::Real=0
)

    ttype = typeof(timeID)
    if ttype <: Array
          datevec = collect(Date(minimum(timeID),1):Month(1):Date(maximum(timeID),12));
    else; datevec = collect(Date(timeID,1):Month(1):Date(timeID,12))
    end

    step = @sprintf("%.2f",eraregionstep(regID,hres))
    fol  = datadir("$dID/$(yr2str(datevec[1]))")
    fnc  = joinpath(fol,"$(dID)-$(regID)x$(step)-tcwvVprcp-$(yrmo2str(datevec[1])).nc");
    ds   = Dataset(fnc); nlon = ds.dim["longitude"]; nlat = ds.dim["latitude"]
    ntcw = ds.dim["tcwv"];
    prcp = zeros(Float64,nlon,nlat,ntcw);
    freq = zeros(Int64,nlon,nlat,ntcw)
    lon  = ds["longitude"][:]; lat = ds["latitude"][:]; tcw = ds["tcwv"][:]
    close(ds)

    for dtii in datevec

        fol  = datadir("$dID/$(yr2str(dtii))")
        fnc  = joinpath(fol,"$(dID)-$(regID)x$(step)-tcwvVprcp-$(yrmo2str(dtii)).nc")
        ds   = Dataset(fnc); prcp_avg = ds["prcp_avg"][:]; bin_frq = ds["bin_frq"][:]
        prcp += prcp_avg .* bin_frq; freq += bin_frq; close(ds)

    end

    @save "$(datadir("compiled/tcwVprcp/$(dID).jld2"))" lon lat tcw prcp freq

end

mkpath(datadir("compiled/tcwVprcp"))
pecompile("era",regID="SEA",timeID=[1980,2019])
pecompile("gpm",regID="SEA",timeID=[2001,2018])
