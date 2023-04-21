using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using GeoRegions
using Dates
using JLD2
using Logging
using NCDatasets
using Printf
using Statistics

include(srcdir("common.jl"))

function pccompile(
    dID::AbstractString;
    regID::AbstractString, timeID::Union{Integer,Vector}=0,
    hres::Real=0
)

    ttype = typeof(timeID)
    if ttype <: Array
          datevec = collect(Date(minimum(timeID),1):Month(1):Date(maximum(timeID),12));
    else; datevec = collect(Date(timeID,1):Month(1):Date(timeID,12))
    end

    @info "$(Dates.now()) - Extracting dimensional data ..."
    step = @sprintf("%.2f",eraregionstep(regID,hres))
    fol  = datadir("$dID/$(yr2str(datevec[1]))")
    fnc  = joinpath(fol,"$(dID)-$(regID)x$(step)-csfVprcp-$(yrmo2str(datevec[1])).nc");
    ds   = Dataset(fnc); nlon = ds.dim["longitude"]; nlat = ds.dim["latitude"]
    ncsf = ds.dim["csf"];

    @info "$(Dates.now()) - Preallocating arrays ..."
    prcp = zeros(Float64,nlon,nlat,ncsf);
    freq = zeros(Int64,nlon,nlat,ncsf)
    lon  = ds["longitude"][:]; lat = ds["latitude"][:]; csf = ds["csf"][:]
    close(ds)

    for dtii in datevec

        @info "$(Dates.now()) - Extracting P-E curve data for $(gregionfullname(regID)) (Horizontal Resolution: $(step)) during $(year(dtii)) $(Dates.monthname(dtii)) ..."
        fol  = datadir("$dID/$(yr2str(dtii))")
        fnc  = joinpath(fol,"$(dID)-$(regID)x$(step)-csfVprcp-$(yrmo2str(dtii)).nc")
        ds   = Dataset(fnc); prcp_avg = ds["prcp_avg"][:]; bin_frq = ds["bin_frq"][:]

        @info "$(Dates.now()) - Adding P-E curve information for $(year(dtii)) $(Dates.monthname(dtii)) ..."
        prcp_avg[isnan.(prcp_avg)] .= 0;
        prcp += prcp_avg .* bin_frq; freq += bin_frq; close(ds)

    end

    @info "$(Dates.now()) - Saving compiled P-C curve information ..."
    @save "$(datadir("compiled/csfVprcp/$(dID).jld2"))" lon lat csf prcp freq

end

function pccompilemonth(
    dID::AbstractString;
    month::Integer,
    regID::AbstractString, timeID::Vector{<:Integer}=0,
    hres::Real=0
)

    ttype = typeof(timeID)
    datevec = collect(Date(minimum(timeID),month):Year(1):Date(maximum(timeID),month))

    @info "$(Dates.now()) - Extracting dimensional data ..."
    step = @sprintf("%.2f",eraregionstep(regID,hres))
    fol  = datadir("$dID/$(yr2str(datevec[1]))")
    fnc  = joinpath(fol,"$(dID)-$(regID)x$(step)-csfVprcp-$(yrmo2str(datevec[1])).nc");
    ds   = Dataset(fnc); nlon = ds.dim["longitude"]; nlat = ds.dim["latitude"]
    ncsf = ds.dim["csf"];

    @info "$(Dates.now()) - Preallocating arrays ..."
    prcp = zeros(Float64,nlon,nlat,ncsf);
    freq = zeros(Int64,nlon,nlat,ncsf)
    lon  = ds["longitude"][:]; lat = ds["latitude"][:]; csf = ds["csf"][:]
    close(ds)

    for dtii in datevec

        @info "$(Dates.now()) - Extracting P-E curve data for $(gregionfullname(regID)) (Horizontal Resolution: $(step)) during $(year(dtii)) $(Dates.monthname(dtii)) ..."
        fol  = datadir("$dID/$(yr2str(dtii))")
        fnc  = joinpath(fol,"$(dID)-$(regID)x$(step)-csfVprcp-$(yrmo2str(dtii)).nc")
        ds   = Dataset(fnc); prcp_avg = ds["prcp_avg"][:]; bin_frq = ds["bin_frq"][:]

        @info "$(Dates.now()) - Adding P-E curve information for $(year(dtii)) $(Dates.monthname(dtii)) ..."
        prcp_avg[isnan.(prcp_avg)] .= 0;
        prcp += prcp_avg .* bin_frq; freq += bin_frq; close(ds)

    end

    @info "$(Dates.now()) - Saving compiled P-C curve information ..."
    @save "$(datadir("compiled/csfVprcp/$(dID)-$month.jld2"))" lon lat csf prcp freq

end

mkpath(datadir("compiled/csfVprcp"))
pccompile("era",regID="SEA",timeID=[1980,2019])
pccompile("gpm",regID="SEA",timeID=[2001,2018])

for mo in 1 : 12
    pccompilemonth("era",month=mo,regID="SEA",timeID=[1980,2019])
    pccompilemonth("gpm",month=mo,regID="SEA",timeID=[2001,2018])
end
