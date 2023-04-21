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

    # if T >= 273.16
    # 	return 611.21 * exp(17.502 * (T-273.16) / (T-32.19))
    # elseif T <= 250.16
    # 	return 611.21 * exp(22.587 * (T-273.16) / (T+0.7))
    # else
    #     α  = ((T - 250.16) / (273.16 - 250.16))^2
    #     ei = 611.21 * exp(22.587 * (T-273.16) / (T+0.7))
    #     ew = 611.21 * exp(17.502 * (T-273.16) / (T-32.19))
    #     return α * ew + (1-α) * ei
    # end

    tb = T - 273.15
    if tb <= 0
    	return exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
    	return exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
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

    
    tvar_vec = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        tvar_vec[ip] = PressureVariable("t",hPa=p[ip])
    end
    t2mvar   = SingleVariable("t2m")
    spvar    = SingleVariable("sp")
    p = convert.(Float32,vcat(0,p*100,0)); np = length(p)

    @info "$(now()) - MaritimeContPWV - Preallocating temporary arrays for calculation ..."
    tint_tair = zeros(Int16,nlon,nlat)
    tint_tsfc = zeros(Int16,nlon,nlat)
    tint_psfc = zeros(Int16,nlon,nlat)
    tflt_tair = zeros(nlon,nlat)
    tflt_tsfc = zeros(nlon,nlat)
    tflt_psfc = zeros(nlon,nlat)
    tmp_es = zeros(nlon,nlat,np-1)
    tmp_ev = zeros(np)
    tmp_ip = zeros(Bool,np)
    tmp_s  = zeros(nlon,nlat,31*24)

    tair_ds = Vector{NCDataset}(undef,np)

    for dt in  e5ds.start : Month(1) : e5ds.stop

        @info "$(now()) - MaritimeContPWV - Calculating $(e5ds.lname) $(evar.vname) data in $(egeo.geo.name) (Horizontal Resolution: $(egeo.gres)) for $(Dates.format(dt,dateformat"yyyy-mm")) ..."

        nhr = daysinmonth(dt) * 24

        for ip = 1 : (np-2)
            tair_ds[ip] = read(e5ds,tvar_vec[ip],egeo,dt)
        end
        tsfc_ds = read(e5ds,t2mvar,egeo,dt)
        psfc_ds = read(e5ds,spvar,egeo,dt)

        nhr = daysinmonth(dt) * 24

        for it = 1 : nhr

            sc = tsfc_ds["t2m"].attrib["scale_factor"]
            of = tsfc_ds["t2m"].attrib["add_offset"]
            mv = tsfc_ds["t2m"].attrib["missing_value"]
            fv = tsfc_ds["t2m"].attrib["_FillValue"]
            NCDatasets.load!(tsfc_ds["t2m"].var,tint_tsfc,:,:,it)
            ERA5Reanalysis.int2real!(tflt_tsfc,tint_tsfc,scale=sc,offset=of,mvalue=mv,fvalue=fv)
    
            sc = psfc_ds["sp"].attrib["scale_factor"]
            of = psfc_ds["sp"].attrib["add_offset"]
            mv = psfc_ds["sp"].attrib["missing_value"]
            fv = psfc_ds["sp"].attrib["_FillValue"]
            NCDatasets.load!(psfc_ds["sp"].var,tint_psfc,:,:,it)
            ERA5Reanalysis.int2real!(tflt_psfc,tint_psfc,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            for ip = 1 : (np-2)
                sc = tair_ds[ip]["t"].attrib["scale_factor"]
                of = tair_ds[ip]["t"].attrib["add_offset"]
                mv = tair_ds[ip]["t"].attrib["missing_value"]
                fv = tair_ds[ip]["t"].attrib["_FillValue"]
                NCDatasets.load!(tair_ds[ip]["t"].var,tint_tair,:,:,it)
                ERA5Reanalysis.int2real!(tflt_tair,tint_tair,scale=sc,offset=of,mvalue=mv,fvalue=fv)
                @views @. tmp_es[:,:,ip+1] = calcT2q(tflt_tair,p[ip+1])
            end

            for ilat = 1 : nlat, ilon = 1 : nlon
                tii = tflt_tsfc[ilon,ilat]
                pii = tflt_psfc[ilon,ilat]
                for ip = 1 : (np-1)
                    tmp_ev[ip] = tmp_es[ilon,ilat,ip]
                    tmp_ip[ip] = p[ip] < pii
                end
                tmp_ev[end] = calcT2q(tii,pii)
                tmp_ip[end] = true
                p[end] = pii
                tmp_s[ilon,ilat,it] = integrate(view(p,tmp_ip),view(tmp_ev,tmp_ip)) / 9.81
            end

        end

        close(tsfc_ds); close(psfc_ds)
        for ip = 1 : (np-2)
            close(tair_ds[ip])
        end

        ERA5Reanalysis.save(view(tmp_s,:,:,1:nhr),dt,e5ds,evar,egeo,lsd)

    end

end
