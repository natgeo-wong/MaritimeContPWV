using ClimateERA
using ClimateSatellite
using Dates
using Logging
using NCDatasets

function pecurve(prcp::AbstractArray,tcwv::AbstractArray,tvec::Vector{<:Real},tsep::Real)


    pvec = zeros(length(tvec)); jj = 0;
    pfrq = zeros(Int32,length(tvec))
    for tii in tvec
        pii = @view prcp[ (tcwv.>(tii-tsep)) .& (tcwv.<(tii+tsep)) ]
        jj = jj + 1; pvec[jj] = mean(pii); pfrq[jj] = length(pii)
    end

    return pvec,pfrq

end

function tcwvVprcp_gpm(
    init::AbstractDict, eroot::AbstractDict, sroot::AbstractString;
    regID::AbstractString="GLB",
    timeID::Union{Integer,Vector}=0
)

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    tmod,tpar,ereg,etime = erainitialize(
        init,
        modID="msfc",parID="tcwv",regID=regID,timeID=timeID
    );
    global_logger(ConsoleLogger(stdout,Logging.Info))

    nlon,nlat = ereg["size"]; lon = ereg["lon"]; lat = ereg["lat"]
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    @info "$(Dates.now()) - Preallocating data arrays to compare precipitation against total column water ..."

    tvec = collect(10:0.5:90); nvec = length(tvec); tstep = (tvec[2]-tvec[1])/2
    pmat = Array{Float32,3}(undef,nlon,nlat,nvec)
    pfrq = Array{Int32,3}(undef,nlon,nlat,nvec)

    @info "$(Dates.now()) - Extracting relevant closest-coordinate points of GPM precipitation for each of the ERA5 total column water grid points ..."

    lon,lat = gpmlonlat(); rlon,rlat,_ = gregiongridvec(reg,lon,lat);
    glon = zeros(nlon); for i = 1 : nlon; glon[i] = argmin(abs.(rlon .- lon[i])) end
    glat = zeros(nlat); for i = 1 : nlat; glat[i] = argmin(abs.(rlat .- lat[i])) end

    for dtii in datevec

        @info "$(Dates.now()) - Extracting ERA5 total column water data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) during $(year(date)) $(Dates.monthname(date)) ..."

        ndy = daysinmonth(dtii)
        tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii); tcwv = tvar[:]*1; close(tds)

        @info "$(Dates.now()) - Extracting GPM Precipitation data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) during $(year(date)) $(Dates.monthname(date)) ..."

        pds,pvar = clisatrawread("gpmimerg","prcp_rate",regID,dtii,sroot);
        prcp  = pvar[:]*1; close(pds)
        itmp1 = Array{Float32,2}(undef,2,24*ndy)
        itmp2 = Array{Float32,1}(undef,24*ndy)

        for ilat = 1 : nlat; ilon = 1 : nlon

            prcpii = @view prcp[rlon[glon[ilon]],rlat[glat[ilat]],:];
            itmp1 .= reshape(prcpii,2,:); itmp2 .= mean(itmp1,dims=1)
            tcwvii = @view tcwv[ilon,ilat,:]
            pmat[ilon,ilat,:],pfrq[ilon,ilat,:] .= pecurve(itmp2,tcwvii,tvec,tstep)

        end

        tcwvVprcpsave(pmat,pfrq,tvec,ereg,dtii,"gpm")

    end

end

function tcwvVprcp_era(
    init::AbstractDict, eroot::AbstractDict;
    regID::AbstractString="GLB",
    timeID::Union{Integer,Vector}=0
)

    global_logger(ConsoleLogger(stdout,Logging.Warn))
    tmod,tpar,ereg,etime = erainitialize(
        init,
        modID="msfc",parID="tcwv",regID=regID,timeID=timeID
    );
    pmod,ppar,____,_____ = erainitialize(init,modID="msfc",parID="prcp_tot");
    global_logger(ConsoleLogger(stdout,Logging.Info))

    nlon,nlat = ereg["size"];
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    @info "$(Dates.now()) - Preallocating data arrays to compare precipitation against total column water ..."

    tvec = collect(10:0.5:90); nvec = length(tvec); tstep = (tvec[2]-tvec[1])/2
    pmat = Array{Float32,3}(undef,nlon,nlat,nvec)
    pfrq = Array{Float32,3}(undef,nlon,nlat,nvec)

    for dtii in datevec

        @info "$(Dates.now()) - Extracting ERA5 total column water and precipitation data for $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) during $(year(date)) $(Dates.monthname(date)) ..."

        tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii); tcwv = tvar[:]*1; close(tds)
        pds,pvar = erarawread(pmod,ppar,ereg,eroot,dtii); prcp = pvar[:]*1; close(pds)

        for ilat = 1 : nlat; ilon = 1 : nlon

            prcpii = @view prcp[ilon,ilat,:]; tcwvii = @view tcwv[ilon,ilat,:]
            pmat[ilon,ilat,:],pfrq[ilon,ilat,:] .= pecurve(prcpii,tcwvii,tvec,tstep)

        end

        tcwvVprcpsave(pmat,pfrq,tvec,ereg,dtii,"era")

    end

end

function tcwvVprcpsave(
    pmat::Array{<:Real,3}, pfrq::Array{<:Real,3}, tvec::Vector{<:Real},
    ereg::Dict, date::TimeType, prefix::AbstractString
)

    @info "$(Dates.now()) - Saving binned averaged precipitation and frequency of bin occurrence in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) ..."

    fnc = joinpath(datadir(),"$(prefix)-$(ereg["fol"])-tcwvVprcp.nc");
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    avgscale,avgoffset = erancoffsetscale(pmat);
    frqscale,frqoffset = erancoffsetscale(pfrq);

    ds.dim["longitude"] = ereg["size"][1];
    ds.dim["latitude"]  = ereg["size"][2];
    ds.dim["tcwv"]      = length(tcwv)

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nctcwv = defVar(ds,"tcwv",Int32,("tcwv",),attrib = Dict(
        "long_name" => "total_column_water_vapour",
        "full_name" => "Total Column Water Vapour",
        "units"     => "kg m^{-2}"
    ))

    ncprcp = defVar(ds,"prcp_avg",Int16,("longitude","latitude","tcwv"),attrib = Dict(
        "long_name" => "averaged_precipitation",
        "full_name" => "Average Precipitation in Bin",
        "units"     => "m",
        "scale_factor"  => avgscale,
        "add_offset"    => avgoffset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
    ))

    ncbfrq = defVar(ds,"bin_frq",Int32,("longitude","latitude","tcwv"),attrib = Dict(
        "long_name" => "bin_frequency",
        "full_name" => "Frequency of Occurrence in Bin",
        "scale_factor"  => frqscale,
        "add_offset"    => frqoffset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    ncprcp[:] = pmat; ncbfrq[:] = pfrq;

    close(ds)

    @info "$(Dates.now()) - Binned averaged precipitation and frequency of bin occurrence in $(gregionfullname(ereg["region"])) (Horizontal Resolution: $(ereg["step"])) for $(year(date)) $(Dates.monthname(date)) has been saved into $(fnc)."

end
