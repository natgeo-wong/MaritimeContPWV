using ClimateERA
using ClimateSatellite
using Dates
using GeoRegions
using Logging
using NCDatasets
using StatsBase

include(srcdir("common.jl"))

function tcwfreq(
    init::AbstractDict, eroot::AbstractDict;
    regID::AbstractString="GLB",
    timeID::Union{Integer,Vector}=0,
    pwv::AbstractRange=0:100
)

    emod,epar,ereg,etime = erainitialize(
        init,
        modID="msfc",parID="tcwv",regID=regID,timeID=timeID
    )

    nlon,nlat = ereg["size"]; npwv = length(pwv)
    elon = ereg["lon"]; elat = ereg["lat"]
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    @info "$(Dates.now()) - Preallocating data arrays to find frequency ..."
    evmmon = zeros(Int32,nlon,nlat,npwv-1)
    evmcom = zeros(Int32,nlon,nlat,npwv-1)

    for dtii in datevec

        @info "$(Dates.now()) - Extracting ERA5 total column water data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) during $(year(dtii)) $(Dates.monthname(dtii)) ..."

        tds,tvar = erarawread(emod,epar,ereg,eroot,dtii); tcw = tvar[:]*1; close(tds)

        for ilat = 1 : nlat, ilon = 1 : nlon

            etcwii = @view tcw[ilon,ilat,:]
            evmmon[ilon,ilat,:,:] .= fit(Histogram,etcwii,pwv).weights
            evmcom[ilon,ilat,:,:] += fit(Histogram,etcwii,pwv).weights

        end

        tcwfreqsave(evmmon,pwv,ereg,dtii)

    end

    tcwfreqsave(evmcom,pwv,ereg)

end

function tcwfreqsave(
    evmfreq::Array{<:Real,3}, pwv::AbstractRange, ereg::Dict, date::TimeType
)

    @info "$(Dates.now()) - Saving binned frequencies for ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) ..."

    fol = datadir("compiled/$(yr2str(date))"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"tcwvfreqsave-$(ereg["fol"])-$(yrmo2str(date)).nc");
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    ds.dim["longitude"] = ereg["size"][1]
    ds.dim["latitude"]  = ereg["size"][2]
    ds.dim["pwv"]       = length(pwv)
    ds.dim["bin"]       = length(pwv) - 1

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    ncpwv = defVar(ds,"pwv",Float32,("pwv",),attrib = Dict(
        "long_name" => "total_column_water_vapour",
        "full_name" => "Total Column Water Vapour",
        "units"     => "kg m^{-2}"
    ))

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","bin"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    ncpwv[:] = collect(pwv); ncbfrq[:] = evmfreq;

    close(ds)

    @info "$(Dates.now()) - Binned frequencies for ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) has been saved into $(fnc)."

end


function tcwfreqsave(
    evmfreq::Array{<:Real,3}, pwv::AbstractRange, ereg::Dict
)

    @info "$(Dates.now()) - Saving binned frequencies for ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates ..."

    fol = datadir("compiled"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"tcwvfreqsave-$(ereg["fol"]).nc");
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    ds.dim["longitude"] = ereg["size"][1]
    ds.dim["latitude"]  = ereg["size"][2]
    ds.dim["pwv"]       = length(pwv)
    ds.dim["bin"]       = length(pwv) - 1

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    ncpwv = defVar(ds,"pwv",Float32,("pwv",),attrib = Dict(
        "long_name" => "total_column_water_vapour",
        "full_name" => "Total Column Water Vapour",
        "units"     => "kg m^{-2}"
    ))

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","bin"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    ncpwv[:] = collect(pwv); ncbfrq[:] = evmfreq;

    close(ds)

    @info "$(Dates.now()) - Binned frequencies for ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates has been saved into $(fnc)."

end
