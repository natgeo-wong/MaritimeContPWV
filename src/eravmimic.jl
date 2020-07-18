using ClimateERA
using ClimateSatellite
using Dates
using GeoRegions
using Logging
using NCDatasets
using StatsBase

include(srcdir("common.jl"))

function eravmimic(
    init::AbstractDict, eroot::AbstractDict, sroot::AbstractString;
    regID::AbstractString="GLB",
    timeID::Union{Integer,Vector}=0,
    pwv::AbstractRange=0:100
)

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    tmod,tpar,ereg,etime = erainitialize(
        init,
        modID="msfc",parID="tcwv",regID=regID,timeID=timeID
    );
    global_logger(ConsoleLogger(stdout,Logging.Info))

    nlon,nlat = ereg["size"]; npwv = length(pwv)
    elon = ereg["lon"]; elat = ereg["lat"]
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    @info "$(Dates.now()) - Preallocating data arrays to compare ERA5 Total Column Water against MIMIC Total Precipitable Water ..."

    evm = zeros(Int32,nlon,nlat,npwv-1,npwv-1);

    for dtii in datevec

        @info "$(Dates.now()) - Extracting ERA5 total column water data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) during $(year(dtii)) $(Dates.monthname(dtii)) ..."

        tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii); tcw = tvar[:]*1; close(tds)

        @info "$(Dates.now()) - Extracting MIMIC total precipitable water data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: 0.25) during $(year(dtii)) $(Dates.monthname(dtii)) ..."

        mds,mvar = clisatrawread("mtpw2m","tpw",dtii,regID,path=sroot);
        tpw = mvar[:]*1; close(mds)

        for ilat = 1 : nlat, ilon = 1 : nlon; mlat = nlat + 1 - ilat

            mtpwii = @view tpw[ilon,mlat,:]
            etcwii = @view tcw[ilon,ilat,:]
            evm[ilon,ilat,:,:] .= fit(Histogram,(mtpwii,etcwii),(pwv,pwv)).weights

        end

        eravmimicsave(evm,pwv,ereg,dtii)

    end

    @info "$(Dates.now()) - Preallocating data arrays to compile the comparison between ERA5 Total Column Water against MIMIC Total Precipitable Water ..."

    evm = zeros(Int32,nlon,nlat,npwv-1,npwv-1);

    for dtii in datevec

        @info "$(Dates.now()) - Extracting binned ERA5 vs MIMIC Total Column Water for $(year(dtii)) $(Dates.monthname(dtii)) ..."

        fol = datadir("eravmimic/$(yr2str(dtii))")
        fnc = joinpath(fol,"eravmimic-$(ereg["fol"])-$(yrmo2str(dtii)).nc")
        ds  = NCDataset(fnc); evm += ds["bin_frq"][:]; close(ds)

    end

    eravmimicsave(evm,pwv,ereg)

    @info "$(Dates.now()) - Preallocating data arrays to find the correlation gridpoint by gridpoint between ERA5 and MIMIC data ..."

    evmcorr = zeros(nlon,nlat);
    nhr  = length(DateTime(etime["Begin"],1):Hour(1):DateTime(etime["End"]+1,1)) - 1
    etmp = zeros(nhr); mtmp = zeros(nhr)

    for ilat = 1 : nlat, ilon = 1 : nlon; ibeg = 1; iend = 0; ; mlat = nlat + 1 - ilat

        for dtii in datevec

            tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii);
            mds,mvar = clisatrawread("mtpw2m","tpw",dtii,regID,path=sroot);

            iend += 24 * daysinmonth(dtii)
            mtmp[ibeg:iend] .= mvar[ilon,mlat,:]
            etmp[ibeg:iend] .= tvar[ilon,ilat,:]

            close(tds); close(mds); ibeg += 24 * daysinmonth(dtii)

        end

        evmcorr[ilon,ilat] .= cor(mtmp,etmp)

    end

    eravmimicsave(evmcorr,ereg)

end

function eravmimicsave(
    evm::Array{<:Real,4}, pwv::AbstractRange,
    ereg::Dict, date::TimeType
)

    @info "$(Dates.now()) - Saving binned frequencies for MIMIC vs ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) ..."

    fol = datadir("eravmimic/$(yr2str(date))"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"eravmimic-$(ereg["fol"])-$(yrmo2str(date)).nc");
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

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","bin","bin"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    ncpwv[:] = collect(pwv); ncbfrq[:] = evm;

    close(ds)

    @info "$(Dates.now()) - Binned frequencies for MIMIC vs ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) has been saved into $(fnc)."

end

function eravmimicsave(
    evm::Array{<:Real,4}, pwv::AbstractRange, ereg::Dict
)

    @info "$(Dates.now()) - Saving binned frequencies for MIMIC vs ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates ..."

    fol = datadir("compiled"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"eravmimic-$(ereg["fol"]).nc");
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

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","bin","bin"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    ncpwv[:] = collect(pwv); ncbfrq[:] = evm;

    close(ds)

    @info "$(Dates.now()) - Binned frequencies for MIMIC vs ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates has been saved into $(fnc)."

end

function eravmimicsave(
    evmcorr::Array{<:Real,4}, ereg::Dict
)

    @info "$(Dates.now()) - Saving correlation for MIMIC vs ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates ..."

    fol = datadir("compiled"); if !isdir(fol); mkpath(fol) end

    fnc = joinpath(fol,"eravmimic-$(ereg["fol"])-corr.nc");
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    ds.dim["longitude"] = ereg["size"][1];
    ds.dim["latitude"]  = ereg["size"][2];

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nccorr = defVar(ds,"rho",Float64,("longitude","latitude"),attrib = Dict(
        "long_name" => "pearson_correlation",
        "full_name" => "Pearson Correlation Coefficient",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]; nccorr[:] = evmcorr;

    close(ds)

    @info "$(Dates.now()) - Correlation for MIMIC vs ERA5 Total Column Water in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for all dates has been saved into $(fnc)."

end
