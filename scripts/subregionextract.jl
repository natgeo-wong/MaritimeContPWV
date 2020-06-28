using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateERA
using GeoRegions

gregioninfoadd(srcdir("gregionsadd.txt"))

init,eroot = erastartup(aID=1,dID=1,path="/n/kuangdss01/lab/ecmwf/");
erasubregion(init,eroot,modID="msfc",parID="tcw",regID="SEA",pregID="TRP");
erasubregion(init,eroot,modID="msfc",parID="tcw",regID="SMT",pregID="TRP");
erasubregion(init,eroot,modID="msfc",parID="tcwv",regID="SEA",pregID="TRP");
erasubregion(init,eroot,modID="msfc",parID="tcwv",regID="SMT",pregID="TRP");

init,eroot = erastartup(aID=1,dID=2,path="/n/kuangdss01/lab/ecmwf/");
erasubregion(init,eroot,modID="msfc",parID="tcw",regID="SMT",pregID="SEA");
erasubregion(init,eroot,modID="msfc",parID="tcwv",regID="SMT",pregID="SEA");

# for yr = 2017 : 2019, mo = 1 : 12
    # clisatdownload("mtpw2m",Date(yr,mo),regions=["TRP"],path="/n/kuangdss01/lab/clisat/");
# end
