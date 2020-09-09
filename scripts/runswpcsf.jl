using DrWatson
@quickactivate "MaritimeContPWV"

include(srcdir("swpcsf.jl")); adderaparams()

for yr = 1979 : 2019
    init,eroot = erastartup(aID=2,dID=1,path="/n/holyscratch01/kuang_lab/",welcome=false)
    swp(init,eroot,eroot,regID="SEA",timeID=yr)
    csf(init,eroot,eroot,regID="SEA",timeID=yr)
end
