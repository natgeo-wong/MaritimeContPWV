using ClimateERA
using Dates
using GeoRegions
using Logging
using NCDatasets
using StatsBase

include(srcdir("common.jl"))

function sstvcsf(
    init::AbstractDict, eroot::AbstractDict;
    regID::AbstractString="GLB", timeID::Union{Integer,Vector}=0,
    nbins::Integer=100
)

    global_logger(ConsoleLogger(stderr,Logging.Warn))
    tmod,tpar,ereg,etime = erainitialize(
        init,
        modID="dsfc",parID="t_skt",regID=regID,timeID=timeID
    );
    cmod,cpar,____,_____ = erainitialize(init,modID="csfc",parID="csf");
    global_logger(ConsoleLogger(stderr,Logging.Info))

    nlon,nlat = ereg["size"];
    datevec = collect(Date(etime["Begin"],1):Month(1):Date(etime["End"],12));

    @info "$(Dates.now()) - Preallocating data arrays to compare precipitation against column saturation fraction ..."

    cvec = collect(0:nbins); ncbin = length(cvec)
    tvec = collect(285:310); ntbin = length(tvec)
    stvc = zeros(Int32,nlon,nlat,ntbin-1,ncbin-1)

    geo = GeoRegion(ereg["region"])

    for dtii in datevec

        @info "$(Dates.now()) - Extracting ERA5 column saturation fraction and skin temperature data for $(geo.name) (Horizontal Resolution: $(ereg["step"])) during $(year(dtii)) $(Dates.monthname(dtii)) ..."

        cds,cvar = erarawread(cmod,cpar,ereg,eroot,dtii); csf = cvar[:]*1; close(cds)
        tds,tvar = erarawread(tmod,tpar,ereg,eroot,dtii); sst = tvar[:]*1; close(tds)

        @info "$(Dates.now()) - Binning ERA5 column saturation fraction and skin temperature data in $(geo.name) (Horizontal Resolution: $(ereg["step"])) during $(year(dtii)) $(Dates.monthname(dtii)) ..."

        for ilat = 1 : nlat, ilon = 1 : nlon

            csfii = @view csf[ilon,ilat,:]
            sstii = @view sst[ilon,ilat,:]
            stvc[ilon,ilat,:,:] += fit(Histogram,(sstii,csfii),(tvec,cvec)).weights

        end

        csf = []
        sst = []

    end

    sstvcsfsave(stvc,tvec,cvec,ereg)

end

function sstvcsfsave(
    stvc::Array{<:Real,4},
    tvec::Vector{<:Real},
    cvec::Vector{<:Real},
    ereg::Dict
)

    geo = GeoRegion(ereg["region"])
    @info "$(Dates.now()) - Saving binned column saturation fraction and sea surface temperature in $(geo.name) (Horizontal Resolution: $(ereg["step"])) ..."

    fnc = "era5-$(ereg["fol"])-sstvcsf.nc"
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

    ds.dim["longitude"] = ereg["size"][1];
    ds.dim["latitude"]  = ereg["size"][2];
    ds.dim["csf"]       = length(cvec)
    ds.dim["skt"]       = length(tvec)
    ds.dim["csfbin"]    = length(cvec) - 1
    ds.dim["sktbin"]    = length(tvec) - 1

    nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    ncsst = defVar(ds,"skt",Float32,("skt",),attrib = Dict(
        "long_name" => "skin_temperature",
        "full_name" => "Skin Temperature",
        "units"     => "K"
    ))

    nccsf = defVar(ds,"csf",Float32,("csf",),attrib = Dict(
        "long_name" => "fraction_column_saturation",
        "full_name" => "Column Saturation Fraction",
        "units"     => "kg m^{-2}"
    ))

    ncbfrq = defVar(
        ds,"bin_frq",Int32,("longitude","latitude","sktbin","csfbin"),
        attrib = Dict(
            "long_name" => "bin_frequency",
            "full_name" => "Frequency of Occurrence in Bin",
    ))

    nclongitude[:] = ereg["lon"]; nclatitude[:] = ereg["lat"]
    ncsst[:] = tvec; nccsf[:] = cvec; ncbfrq[:] = stvc;

    close(ds)

    @info "$(Dates.now()) - Binned column saturation fraction against sea surface temperature in $(geo.name) (Horizontal Resolution: $(ereg["step"])) has been saved into $(fnc)."

end