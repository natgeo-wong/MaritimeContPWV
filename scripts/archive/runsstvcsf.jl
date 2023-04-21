using DrWatson
@quickactivate "MaritimeContPWV"

include(srcdir("sstvcsf.jl")); adderaparams()

init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/");
sstvcsf(init,eroot,regID="TRP",timeID=1980)
