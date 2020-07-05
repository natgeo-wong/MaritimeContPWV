using DrWatson
@quickactivate "MaritimeContPWV"
using ClimateERA
# using ClimateSatellite

init,eroot = erastartup(aID=1,dID=1,path="/n/kuangdss01/lab/");
# eradownload(init,eroot,modID="msfc",parID="tcw",regID="TRP");
eradownload(init,eroot,modID="msfc",parID="tcwv",regID="TRP");

init,eroot = erastartup(aID=1,dID=2,path="/n/kuangdss01/lab/");
eradownload(init,eroot,modID="msfc",parID="tcw",regID="SEA");
eradownload(init,eroot,modID="msfc",parID="tcwv",regID="SEA");
