using ERA5Reanalysis
using Logging
using NCDatasets
using NumericalIntegration

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

    for dt in  e5ds.start : Month(1) : e5ds.stop

        wds = read(e5ds,wvar,egeo,dt)
        sds = read(e5ds,svar,egeo,dt)

        tcwv = nomissing(wds["tcwv"][:],NaN)
        scwv = nomissing(sds["scwv"][:],NaN)
        csf  = tcwv ./ scwv

        close(wds)
        close(sds)

        ERA5Reanalysis.save(csf,dt,e5ds,evar,egeo,lsd)

    end

end
