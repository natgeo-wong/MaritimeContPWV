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
    p = era5Pressures(); p = p[p.>=p_top]; p = p[p.<=p_bot]

    qvar_vec = PressureVariable.("q",hPa=p)
    rvar_vec = PressureVariable.("r",hPa=p)
    d2mvar   = SingleVariable("d2m")
    t2mvar   = SingleVariable("t2m")
    spvar    = SingleVariable("sp")
    tcwvvar  = SingleVariable("tcwv")
    p = convert.(Float32,vcat(0,p,0)); np = length(p)

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
    tmp_es = zeros(nlon,nlat,np)
    tmp_s  = zeros(nlon,nlat,31*24)

    qds = Vector{NCDataset}(undef,np)
    rds = Vector{NCDataset}(undef,np)

    for dt in  e5ds.start : Month(1) : e5ds.stop

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

            @info "$(now()) - MaritimeContPWV - Calculation of SCWV at $it ..."

            sc = tds["t2m"].attrib["scale_factor"]
            of = tds["t2m"].attrib["add_offset"]
            mv = tds["t2m"].attrib["missing_value"]
            fv = tds["t2m"].attrib["_FillValue"]
            NCDatasets.load!(tds["t2m"].var,tint_t,:,:,it)
            int2real!(tflt_t,tint_t,scale=sc,offset=of,mvalue=mv,fvalue=fv)
    
            sc = dds["d2m"].attrib["scale_factor"]
            of = dds["d2m"].attrib["add_offset"]
            mv = dds["d2m"].attrib["missing_value"]
            fv = dds["d2m"].attrib["_FillValue"]
            NCDatasets.load!(dds["d2m"].var,tint_d,:,:,it)
            int2real!(tflt_d,tint_d,scale=sc,offset=of,mvalue=mv,fvalue=fv)
    
            sc = pds["sp"].attrib["scale_factor"]
            of = pds["sp"].attrib["add_offset"]
            mv = pds["sp"].attrib["missing_value"]
            fv = pds["sp"].attrib["_FillValue"]
            NCDatasets.load!(pds["sp"].var,tint_p,:,:,it)
            int2real!(tflt_p,tint_p,scale=sc,offset=of,mvalue=mv,fvalue=fv)

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
                tii = tflt_tii[ilon,ilat]
                dii = tflt_dii[ilon,ilat]
                pii = tflt_pii[ilon,ilat]
                tmp_es[end] = calcT2q(tii,pii) / calcT2q(dii,pii)
                tmp_s[ilon,ilat,it] = integrate(p,tmp_es)
            end

        end

        close(tds); close(dds); close(wds); close(pds)
        for qdsii in qds
            close(qdsii)
        end
        for rdsii in rds
            close(rdsii)
        end

        save(view(tmp_s,:,:,1:nhr),dt,e5ds,evar,egeo,lsd)

    end

end
