using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using GeoRegions

include(srcdir("eravssugar.jl"))

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/")
gregioninfoadd(srcdir("gregionsadd.txt"))
eravsugar(init,eroot)
