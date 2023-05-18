using DrWatson
@quickactivate "MaritimeContPWV"
using ERA5Reanalysis

addGeoRegions(srcdir("addRectRegion_TRP.txt"))

e5ds = ERA5Hourly(start=Date(1979),stop=Date(2021,12),path=datadir())
sgeo = GeoRegion("DTP")
egeo = ERA5Region("TRP")
evar = SingleVariable("csf")

extract(sgeo,e5ds,evar,egeo)