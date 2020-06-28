using DrWatson
@quickactivate "MaritimeContPWV"

include(srcdir("gnss.jl"))

function gnssresortall()

    allinfo = retrieveginfo(); stns = allinfo[:,1]
    for stn in stns; gstationresort(gstationinfo(stn,allinfo)) end

end

# gnssresortall()
for stn in ["NTUS","SING"]; gstationresort(gstationinfo(stn,retrieveginfo())) end
