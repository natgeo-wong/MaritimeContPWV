using DrWatson
@quickactivate "MaritimeContPWV"

using ClimateERA
using GeoRegions

include(srcdir("eravmimic.jl")); droot = "/n/kuangdss01/lab/"

init,eroot = erastartup(aID=2,dID=1,path=droot);
eravmimic(init,eroot,droot,regID="SEA",timeID=[2017,2019])
eramimicrho(init,eroot,droot,regID="SEA",timeID=[2017,2019])
