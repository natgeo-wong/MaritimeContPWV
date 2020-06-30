using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateERA
# using ClimateSatellite

init,eroot = erastartup(aID=1,dID=1,path="/n/kuangdss01/lab/ecmwf/");
# eradownload(init,eroot,modID="msfc",parID="tcw",regID="TRP");
eradownload(init,eroot,modID="msfc",parID="tcwv",regID="TRP");

init,eroot = erastartup(aID=1,dID=2,path="/n/kuangdss01/lab/ecmwf/");
eradownload(init,eroot,modID="msfc",parID="tcw",regID="SEA");
eradownload(init,eroot,modID="msfc",parID="tcwv",regID="SEA");

# for yr = 2017 : 2019, mo = 1 : 12
    # clisatdownload("mtpw2m",Date(yr,mo),regions=["TRP"],path="/n/kuangdss01/lab/clisat/");
# end
