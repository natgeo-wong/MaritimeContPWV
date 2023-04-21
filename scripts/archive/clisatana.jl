using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateSatellite

droot = "/n/kuangdss01/lab/"

rvec = ["TRP","SEA","SMT"]; ddir = "/n/kuangdss01/lab/"

for yr = 2001 : 2018, reg in rvec
    clisatanalysis("gpmimerg",yr,varname="prcp_rate",region=reg,path=ddir);
end

for yr = 2001 : 2019, reg in rvec
    clisatanalysis("gpmlate",yr,varname="prcp_rate",region=reg,path=ddir);
end

for yr = 2017 : 2019, reg in rvec
    clisatanalysis("mtpw2m",yr,varname="tpw",region=reg,path=ddir);
end
