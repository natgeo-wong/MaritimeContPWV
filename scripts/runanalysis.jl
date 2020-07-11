using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateERA

include(srcdir("tcwvVprcp.jl"));
droot = "/n/kuangdss01/lab/"

init,eroot = erastartup(aID=2,dID=1,path=droot);
tcwvVprcp_gpm(init,eroot,droot,regID="SEA",timeID=[2001,2019])
tcwvVprcp_era(init,eroot,regID="SEA",timeID=[1980,2019])

init,eroot = erastartup(aID=2,dID=2,path=droot);
tcwvVprcp_era(init,eroot,regID="SEA",timeID=[1980,2019])
