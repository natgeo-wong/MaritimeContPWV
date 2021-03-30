using Dates
using Statistics

function ncoffsetscale(data::Array{<:Real})

    dataii = @view data[.!isnan.(data)]

    dmax = maximum(dataii); dmin = minimum(dataii);
    scale = (dmax-dmin) / 65533;
    offset = (dmax+dmin-scale) / 2;

    return scale,offset

end

yr2str(date::TimeType)   = Dates.format(date,dateformat"yyyy")
yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")

adderaparams() = eraparameteradd(srcdir("eraparams.txt"));
