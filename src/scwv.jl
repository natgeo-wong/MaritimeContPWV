using ERA5Reanalysis
using Logging
using NCDatasets
using NumericalIntegration

if !isSingle("scwv",throw=false)
    SingleVariable(
        varID = "scwv",
        lname = "saturated_column_water_vapor",
        vname = "Saturated Column Water Vapor",
        units = "kg m**-2",
        inCDS = false
    )
end

function calcT2e(T::Real)

    if T >= 273.16
    	return 611.21 * exp(17.502 * (T-273.16) / (T-32.19))
    elseif T <= 250.16
    	return 611.21 * exp(22.587 * (T-273.16) / (T+0.7))
    else
        α  = ((T - 250.16) / (273.16 - 250.16))^2
        ei = 611.21 * exp(22.587 * (T-273.16) / (T+0.7))
        ew = 611.21 * exp(17.502 * (T-273.16) / (T-32.19))
        return α * ew + (1-α) * ei
    end

end

function calcT2q(T::Real,p::Real)
    e = calcT2e(T)
    return e * 0.621981 / (p - 0.378019 * e)
end

function calculatescwv(
    e5ds :: ERA5Dataset,
    egeo :: ERA5Region;
    p_bot :: Int = 1000,
    p_top :: Int = 1,
)

    evar = SingleVariable("scwv")

    lsd = getLandSea(e5ds,egeo);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)
    p = era5Pressures(); p = p[p.>=p_top]; p = p[p.<=p_bot]; np = length(p)

    
    qvar_vec = Vector{PressureVariable}(undef,np)
    rvar_vec = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        qvar_vec[ip] = PressureVariable("q",hPa=p[ip])
        rvar_vec[ip] = PressureVariable("r",hPa=p[ip])
    end
    d2mvar   = SingleVariable("d2m")
    t2mvar   = SingleVariable("t2m")
    spvar    = SingleVariable("sp")
    p = convert.(Float32,vcat(0,p*100,0)); np = length(p)

    @info "$(now()) - MaritimeContPWV - Preallocating temporary arrays for calculation ..."
    tint_q = zeros(Int16,nlon,nlat)
    tint_r = zeros(Int16,nlon,nlat)
    tint_t = zeros(Int16,nlon,nlat)
    tint_d = zeros(Int16,nlon,nlat)
    tint_p = zeros(Int16,nlon,nlat)
    tflt_q = zeros(nlon,nlat)
    tflt_r = zeros(nlon,nlat)
    tflt_t = zeros(nlon,nlat)
    tflt_d = zeros(nlon,nlat)
    tflt_p = zeros(nlon,nlat)
    tmp_es = zeros(nlon,nlat,np-1)
    tmp_ev = zeros(np)
    tmp_ip = zeros(Bool,np)
    tmp_s  = zeros(nlon,nlat,31*24)

    qds = Vector{NCDataset}(undef,np)
    rds = Vector{NCDataset}(undef,np)

    for dt in  e5ds.start : Month(1) : e5ds.stop

        @info "$(now()) - MaritimeContPWV - Calculating $(e5ds.lname) $(evar.vname) data in $(egeo.geo.name) (Horizontal Resolution: $(egeo.gres)) for $(Dates.format(dt,dateformat"yyyy-mm")) ..."

        nhr = daysinmonth(dt) * 24

        for ip = 1 : (np-2)
            qds[ip] = read(e5ds,qvar_vec[ip],egeo,dt)
            rds[ip] = read(e5ds,rvar_vec[ip],egeo,dt)
        end
        tds = read(e5ds,t2mvar,egeo,dt)
        dds = read(e5ds,d2mvar,egeo,dt)
        pds = read(e5ds,spvar,egeo,dt)

        nhr = daysinmonth(dt) * 24

        for it = 1 : nhr

            sc = tds["t2m"].attrib["scale_factor"]
            of = tds["t2m"].attrib["add_offset"]
            mv = tds["t2m"].attrib["missing_value"]
            fv = tds["t2m"].attrib["_FillValue"]
            NCDatasets.load!(tds["t2m"].var,tint_t,:,:,it)
            ERA5Reanalysis.int2real!(tflt_t,tint_t,scale=sc,offset=of,mvalue=mv,fvalue=fv)
    
            sc = dds["d2m"].attrib["scale_factor"]
            of = dds["d2m"].attrib["add_offset"]
            mv = dds["d2m"].attrib["missing_value"]
            fv = dds["d2m"].attrib["_FillValue"]
            NCDatasets.load!(dds["d2m"].var,tint_d,:,:,it)
            ERA5Reanalysis.int2real!(tflt_d,tint_d,scale=sc,offset=of,mvalue=mv,fvalue=fv)
    
            sc = pds["sp"].attrib["scale_factor"]
            of = pds["sp"].attrib["add_offset"]
            mv = pds["sp"].attrib["missing_value"]
            fv = pds["sp"].attrib["_FillValue"]
            NCDatasets.load!(pds["sp"].var,tint_p,:,:,it)
            ERA5Reanalysis.int2real!(tflt_p,tint_p,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            for ip = 1 : (np-2)
                sc = qds[ip]["q"].attrib["scale_factor"]
                of = qds[ip]["q"].attrib["add_offset"]
                mv = qds[ip]["q"].attrib["missing_value"]
                fv = qds[ip]["q"].attrib["_FillValue"]
                NCDatasets.load!(qds[ip]["q"].var,tint_q,:,:,it)
                ERA5Reanalysis.int2real!(tflt_q,tint_q,scale=sc,offset=of,mvalue=mv,fvalue=fv)
                sc = rds[ip]["r"].attrib["scale_factor"]
                of = rds[ip]["r"].attrib["add_offset"]
                mv = rds[ip]["r"].attrib["missing_value"]
                fv = rds[ip]["r"].attrib["_FillValue"]
                NCDatasets.load!(rds[ip]["r"].var,tint_q,:,:,it)
                ERA5Reanalysis.int2real!(tflt_r,tint_r,scale=sc,offset=of,mvalue=mv,fvalue=fv)
                @views @. tmp_es[:,:,ip+1] = tflt_q * tflt_r / 100
            end

            for ilat = 1 : nlat, ilon = 1 : nlon
                tii = tflt_t[ilon,ilat]
                dii = tflt_d[ilon,ilat]
                pii = tflt_p[ilon,ilat]
                for ip = 1 : (np-1)
                    tmp_ev[ip] = tmp_es[ilon,ilat,ip]
                    tmp_ip[ip] = p[ip] < pii
                end
                tmp_ev[end] = calcT2q(tii,pii) / calcT2q(dii,pii)
                tmp_ip[end] = true
                tmp_s[ilon,ilat,it] = integrate(view(p,tmp_ip),view(tmp_ev,tmp_ip))
            end

        end

        close(tds); close(dds); close(pds)
        for ip = 1 : (np-2)
            close(qds[ip])
            close(rds[ip])
        end

        ERA5Reanalysis.save(view(tmp_s,:,:,1:nhr),dt,e5ds,evar,egeo,lsd)

    end

end
