using DrWatson
@quickactivate "MaritimeContPWV"
using ERA5Reanalysis

addGeoRegions(srcdir("addRectRegion_TRP.txt"))

e5ds = ERA5Hourly(dtbeg=Date(1979),dtend=Date(2021,12),eroot=datadir())
sgeo = GeoRegion("DTP")
egeo = ERA5Region("TRP")
evar = SingleVariable("csf")

extract(sgeo,e5ds,evar,egeo)