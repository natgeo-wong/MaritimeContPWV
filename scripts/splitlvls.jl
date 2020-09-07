using ClimateERA
init,eroot = erastartup(aID=2,dID=1,path="/n/kuangdss01/lab/");
eralvlsplit(init,eroot,modID="dpre",parID="t_air",regID="SEA")
