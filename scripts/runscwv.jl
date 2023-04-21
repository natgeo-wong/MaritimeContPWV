using DrWatson
@quickactivate "MaritimeContPWV"

using ERA5Reanalysis

include(srcdir("scwv.jl"))

e5ds = ERA5Hourly(start=Date(1979),stop=Date(2021,12),path=datadir())
egeo = ERA5Region("TRP",gres=0.25)

calculatescwv(e5ds,egeo,p_top=10)