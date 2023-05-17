using ERA5Reanalysis
using NCDatasets

if !isSingle("csf",throw=false)
    SingleVariable(
        varID = "csf",
        lname = "column_saturation_fraction",
        vname = "Column Saturation Fraction",
        units = "0-1",
        inCDS = false
    )
end

function calculatecsf(
    e5ds :: ERA5Dataset,
    egeo :: ERA5Region
)

    evar = SingleVariable("csf")
    wvar = SingleVariable("tcwv")
    svar = SingleVariable("scwv")
    lsd  = getLandSea(e5ds,egeo)
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    tint_tcwv = zeros(Int16,nlon,nlat,31*24)
    tint_scwv = zeros(Int16,nlon,nlat,31*24)
    tflt_tcwv = zeros(nlon,nlat,31*24)
    tflt_scwv = zeros(nlon,nlat,31*24)
    tflt_csf  = zeros(nlon,nlat,31*24)

    for dt in  e5ds.start : Month(1) : e5ds.stop

        nhr = daysinmonth(dt) * 24
        tcwv_int = view(tint_tcwv,:,:,1:nhr)
        scwv_int = view(tint_scwv,:,:,1:nhr)
        tcwv_flt = view(tflt_tcwv,:,:,1:nhr)
        scwv_flt = view(tflt_scwv,:,:,1:nhr)
        csf_flt  = view(tflt_csf ,:,:,1:nhr)

        wds = read(e5ds,wvar,egeo,dt)
        sds = read(e5ds,svar,egeo,dt)

        sc = wds["tcwv"].attrib["scale_factor"]
        of = wds["tcwv"].attrib["add_offset"]
        mv = wds["tcwv"].attrib["missing_value"]
        fv = wds["tcwv"].attrib["_FillValue"]
        NCDatasets.load!(wds["tcwv"].var,tcwv_int,:,:,:)
        ERA5Reanalysis.int2real!(tcwv_flt,tcwv_int,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["scwv"].attrib["scale_factor"]
        of = sds["scwv"].attrib["add_offset"]
        mv = sds["scwv"].attrib["missing_value"]
        fv = sds["scwv"].attrib["_FillValue"]
        NCDatasets.load!(sds["scwv"].var,scwv_int,:,:,:)
        ERA5Reanalysis.int2real!(scwv_flt,scwv_int,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        @. csf_flt = tcwv_flt ./ scwv_flt

        close(wds)
        close(sds)

        ERA5Reanalysis.save(csf_flt,dt,e5ds,evar,egeo,lsd)

    end

end
