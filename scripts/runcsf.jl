using DrWatson
@quickactivate "MaritimeContPWV"

include(srcdir("calccsf.jl")); adderaparams()

for yr = 1979 : 2019
    init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/",welcome=false)
    csf(init,eroot,eroot,regID="SEA",timeID=yr)
    eraanalysis(init,eroot,modID="csfc",parID="csf",regID="SEA")
end
