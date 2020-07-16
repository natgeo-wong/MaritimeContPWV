using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using GeoRegions

include(srcdir("gnss.jl"))
include(srcdir("eravgnss.jl"))

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/")
mkpath(datadir("compiled/era5"));
mkpath(datadir("compiled/erai"));
mkpath(datadir("compiled/era6"));
gregioninfoadd(srcdir("gregionsadd.txt"))

gstns = retrieveginfo()[:,1]
for gstn in gstns
    init,croot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/",welcome=false)
    init,proot = erastartup(
        aID=2,dID=1,
        path="/n/kuangdss01/users/nwong/PiPWV",
        welcome=false
    )
    eravgnss(gstn,init,croot,proot)
end
