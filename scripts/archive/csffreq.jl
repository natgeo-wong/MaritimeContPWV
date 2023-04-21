using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateERA

include(srcdir("csffreq.jl"))

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/")
csffreq(init,eroot,regID="SEA",timeID=[1980,2019])
