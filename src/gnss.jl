using Dates
using DelimitedFiles
using Logging
using NCDatasets

function stationvalid(dateint::Integer)

    yr = floor(dateint/10000)
    mo = floor(dateint/100) - yr * 100
    dy = rem(dateint,100)

    if dateint != 99999999
          return Date(yr,mo,dy)
    else; return Date(Dates.now())
    end

end

function retrieveinfo()

    fdir = srcdir("SuGAr.sites")
    @info "$(Dates.now()) - Retrieving SuGAr station info from file $(fdir)"
    info = readdlm(fdir,',',comments=true,comment_char='#')

    info[:,5] .= stationvalid.(info[:,5])
    info[:,6] .= stationvalid.(info[:,6])

    return info

end

function stationinfo(stn::AbstractString,allinfo::AbstractArray)

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

function stationresort(sinfo::AbstractDict)

    dbeg = sinfo["start"]; yrbeg = Year(dbeg); mobeg = Month(dbeg)
    dend = sinfo["start"]; yrend = Year(dend); moend = Month(dend)
    dvec = DateTime(yrbeg,mobeg) : Month(1) : DateTime(yrend,moend);

    for id in dvec

        dtvec = id : Minute(10) : (id + Dates.Month(1))

    end

end
