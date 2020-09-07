using DrWatson
@quickactivate "MaritimeContPWV"

include(srcdir("swpcsf.jl"))

init,eroot = erastartup(aID=2,dID=1,path="/n/holyscratch01/kuang_lab/");
adderaparams();
swp(init,eroot,eroot,regID="SEA");
csf(init,eroot,eroot,regID="SEA")
