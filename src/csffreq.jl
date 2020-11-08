using ClimateERA
using ClimateSatellite
using Dates
using GeoRegions
using Logging
using NCDatasets
using StatsBase

include(srcdir("common.jl"))

function csffreq(
    init::AbstractDict, eroot::AbstractDict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    nbins::Integer=100
)

    emod,epar,ereg,etime = erainitialize(
        init,
        modID="csfc",parID="csf",regID=regID,timeID=timeID
    )

    nlon,nlat = ereg["size"]; elon = ereg["lon"]; elat = ereg["lat"]
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    @info "$(Dates.now()) - Preallocating data arrays to find frequency ..."
    csfvec = (0:nbins)/nbins
    csfspt = zeros(Int32,nlon,nlat,nbins)
    csfcum = zeros(Int32,nlon,nlat,nbins)

    if !isdir(datadir("compiled/csffreq")); mkpath(datadir("compiled/csffreq")) end

    for dtii in datevec

        @info "$(Dates.now()) - Extracting ERA5 column saturated fraction data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) during $(year(dtii)) $(Dates.monthname(dtii)) ..."

        cds,cvar = erarawread(emod,epar,ereg,eroot,dtii); csf = cvar[:]*1; close(cds)

        for ilat = 1 : nlat, ilon = 1 : nlon

            csfii = @view csf[ilon,ilat,:]
            csfspt[ilon,ilat,:,:] .= fit(Histogram,csfii,csfvec).weights
            csfcum[ilon,ilat,:,:] += fit(Histogram,csfii,csfvec).weights

        end

        csffreqsave(csfspt,csfvec,ereg,dtii)

    end

    csffreqsave(csfcum,csfvec,ereg)

end

function csffreqsave(
    csfspt::Array{<:Real,3}, csfvec::AbstractRange, ereg::Dict, date::TimeType
)

    @info "$(Dates.now()) - Saving binned frequencies for ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) ..."

    fol = datadir("compiled/csffreq/$(yr2str(date))"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"csffreqsave-$(ereg["fol"])-$(yrmo2str(date)).nc");
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    ds.dim["longitude"] = ereg["size"][1]
    ds.dim["latitude"]  = ereg["size"][2]
    ds.dim["csf"]       = length(csfvec)
    ds.dim["bin"]       = length(csfvec) - 1

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nccsf = defVar(ds,"csf",Float32,("csf",),attrib = Dict(
        "long_name" => "fraction_column_saturation",
        "full_name" => "Column Saturation Fraction",
        "units"     => "kg m^{-2}"
    ))

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","bin"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    nccsf[:] = collect(csfvec); ncbfrq[:] = csfspt;

    close(ds)

    @info "$(Dates.now()) - Binned frequencies for ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) has been saved into $(fnc)."

end


function csffreqsave(
    csfcum::Array{<:Real,3}, csfvec::AbstractRange, ereg::Dict
)

    @info "$(Dates.now()) - Saving binned frequencies for ERA5 Column Saturation Fraction in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates ..."

    fol = datadir("compiled"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"csffreqsave-$(ereg["fol"]).nc");
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    ds.dim["longitude"] = ereg["size"][1]
    ds.dim["latitude"]  = ereg["size"][2]
    ds.dim["csf"]       = length(csfvec)
    ds.dim["bin"]       = length(csfvec) - 1

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nccsf = defVar(ds,"csf",Float32,("csf",),attrib = Dict(
        "long_name" => "fraction_column_saturation",
        "full_name" => "Column Saturation Fraction",
        "units"     => "kg m^{-2}"
    ))

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","bin"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    nccsf[:] = collect(csfvec); ncbfrq[:] = csfcum;

    close(ds)

    @info "$(Dates.now()) - Binned frequencies for ERA5 Column Saturation Fraction in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates has been saved into $(fnc)."

end
