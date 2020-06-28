using Dates
using DelimitedFiles
using Logging
using NCDatasets

function ncoffsetscale(data::Array{<:Real})

    dmax = maximum(data); dmin = minimum(data);
    scale = (dmax-dmin) / 65533;
    offset = (dmax+dmin-scale) / 2;

    return scale,offset

end

yr2str(date::TimeType)   = Dates.format(date,dateformat"yyyy")
yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")

function stationvalid(dateint::Integer)

    yr = floor(dateint/10000)
    mo = floor(dateint/100) - yr * 100
    dy = rem(dateint,100)

    if dateint != 99999999
          return Date(yr,mo,dy)
    else; return Date(Dates.now())
    end

end

function retrieveginfo()

    fdir = srcdir("SuGAr.sites")
    @info "$(Dates.now()) - Retrieving SuGAr station info from file $(fdir)"
    info = readdlm(fdir,',',comments=true,comment_char='#')

    info[:,5] .= stationvalid.(info[:,5])
    info[:,6] .= stationvalid.(info[:,6])

    return info

end

function gstationinfo(stn::AbstractString,allinfo::AbstractArray)

    stninfo = Dict{AbstractString,Any}()

    istn = allinfo[:,1] .== stn
    stninfo["name"]      = allinfo[istn,1][1]
    stninfo["longitude"] = allinfo[istn,2][1]
    stninfo["latitude"]  = allinfo[istn,3][1]
    stninfo["height"]    = allinfo[istn,4][1]
    stninfo["start"]     = allinfo[istn,5][1]
    stninfo["end"]       = allinfo[istn,6][1]

    return stninfo

end

function gstationresort(sinfo::AbstractDict)

    dbeg = sinfo["start"]; yrbeg = Year(dbeg); mobeg = Month(dbeg)
    dend = sinfo["end"]; yrend = Year(dend); moend = Month(dend)
    dvec = DateTime(yrbeg,mobeg) : Month(1) : DateTime(yrend,moend);

    gdata = readdlm(datadir("gnss/$(sinfo["name"]).tdpzwd"),Any,comments=true);
    gdate = collect(DateTime(2000,1,1,12,0,0) .+ Second.(@view gdata[:,2]));
    gyear = Dates.year.(gdate); gmonth = Dates.month.(gdate);

    for id in dvec

        iyear = Dates.year(id); imonth = Dates.month(id)
        dtvec = collect(id : Minute(10) : (id + Dates.Month(1))); pop!(dtvec)
        idata = @view gdata[(gyear .== iyear) .& (gmonth .== imonth),3:4]
        idate = @view gdate[(gyear .== iyear) .& (gmonth .== imonth)]

        ndy = daysinmonth(id); ii = 0; ei = 0;
        zdata = zeros(144,ndy); sdata = zeros(144,ndy);

        for idt in dtvec; ii += 1

            jj = findfirst(isequal.(idt,idate));
            if isnothing(jj)
                  zdata[ii] = NaN; sdata[ii] = NaN; ei += 1;
            else; zdata[ii] = idata[jj,1]; sdata[ii] = idata[jj,2]
            end

        end

        if ei != length(dtvec)
            gstationsave(zdata,sdata,dtvec,sinfo)
        else
            @info "$(Dates.now()) - No valid ZWD data is available from the $(sinfo["name"]) station in $(monthname(id)) $(year(id))"
        end

    end

end

function gstationsave(
    zwd::AbstractArray, sig::AbstractArray, dt::AbstractArray,
    sinfo::AbstractDict
)

    fol = datadir("gnss/$(sinfo["name"])/$(yr2str(dt[1]))");
    if !isdir(fol); mkpath(fol); end

    fnc = datadir("$fol/$(sinfo["name"])-$(yrmo2str(dt[1])).nc")
    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end

    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"  => "CF-1.6",
        "Date Created" => "$(Dates.now())"
    ))

    ds.dim["coordinate"] = 1
    ds.dim["step"] = 144
    ds.dim["days"] = size(zwd,2)
    ds.dim["time"] = size(zwd,2) * 144

    nclon = defVar(ds,"longitude",Float64,("coordinate",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float64,("coordinate",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nctime = defVar(ds,"time",Int32,("time",),attrib = Dict(
        "units"     => "minutes since $(Date(dt[1])) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    nczwd = defVar(ds,"zwd",Float32,("step","days"),attrib = Dict(
        "long_name"     => "zenith_wet_delay",
        "full_name"     => "Zenith Wet Delay",
        "units"         => "m",
    ))

    ncsig = defVar(ds,"sigma",Float32,("step","days"),attrib = Dict(
        "long_name"     => "uncertainty_zenith_wet_delay",
        "full_name"     => "Uncertainty in Zenith Wet Delay",
        "units"         => "m",
    ))

    nclon[:]  = sinfo["longitude"]; nclat[:] = sinfo["latitude"]
    nctime[:] = (collect(1:length(dt)).-1) * 10;
    nczwd[:]  = zwd; ncsig[:] = sig;

end
