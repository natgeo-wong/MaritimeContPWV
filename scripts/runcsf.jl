using DrWatson
@quickactivate "MaritimeContPWV"

using ERA5Reanalysis

include(srcdir("scwv.jl"))

e5ds = ERA5Hourly(start=Date(2021),stop=Date(2021),path=datadir())
egeo = ERA5Region("TRP",gres=0.25)

calculatecsf(e5ds,egeo)