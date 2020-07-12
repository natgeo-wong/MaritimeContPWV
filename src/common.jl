using Statistics

function ncoffsetscale(data::Array{<:Real})

    dataii = @view data[.!isnan.(data)]

    dmax = maximum(dataii); dmin = minimum(dataii);
    scale = (dmax-dmin) / 65533;
    offset = (dmax+dmin-scale) / 2;

    return scale,offset

end
