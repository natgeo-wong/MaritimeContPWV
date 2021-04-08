using DrWatson
@quickactivate "MaritimeContPWV"

include(srcdir("calcswp.jl")); adderaparams()

for yr = 1979 : 2019
    init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/",welcome=false)
    swp(init,eroot,eroot,regID="SEA",timeID=yr)
    eraanalysis(init,eroot,modID="csfc",parID="swp",regID="SEA",plvls="sfc",timeID=yr)
end

eracompile(init,eroot,modID="csfc",parID="swp",regID="SEA",plvls="sfc")