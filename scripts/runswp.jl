using DrWatson
@quickactivate "PiPWV"

include(srcdir("swp.jl"))

init,eroot = erastartup(aID=2,dID=1,path="/n/holyscratch01/kuang_lab/");
adderaparams(); swp(init,eroot,proot,regID="SEA")
