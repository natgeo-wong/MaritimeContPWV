using DrWatson
@quickactivate "MaritimeContPWV"
using NCDatasets

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plotfreqhist(
    thr::Real=0.5
)

    ds = NCDataset(datadir("compiled/tcwfreq-SEAx0.25.nc"))
    freqall = ds["bin_frq"][:]; pwv = ds["pwv"][:]; close(ds)
    pwv = pwv[1:(end-1)] .+ 0.5; nbins = length(pwv)

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
    lsmsea = lsm .< thr; lsmlnd = lsm .> thr;

    freqsea = dropdims(mean(freqall .* lsmsea,dims=(1,2)),dims=(1,2));
    freqlnd = dropdims(mean(freqall .* lsmlnd,dims=(1,2)),dims=(1,2));
    freqall = dropdims(mean(freqall,dims=(1,2)),dims=(1,2));

    pplt.close(); f,axs = pplt.subplots(nrows=1,aspect=2,axwidth=4);

    axs[1].plot(pwv,freqall / sum(freqall) * nbins,lw=1,label="All",legend="ul",c="k")
    axs[1].plot(pwv,freqsea / sum(freqsea) * nbins,lw=1,label="Sea",legend="ul",c="b")
    axs[1].plot(pwv,freqlnd / sum(freqlnd) * nbins,lw=1,label="Land",legend="ul",c="g")

    # axs[2].plot(pwv,freqall,lw=1,label="All",legend="ul",c="k")
    # axs[2].plot(pwv,freqsea,lw=1,label="Sea",legend="ul",c="b")
    # axs[2].plot(pwv,freqlnd,lw=1,label="Land",legend="ul",c="g")

    axs[1].format(
        xlim=(0,80),xlabel="Precipitable Water / mm",
        ylim=(0,6),ylabel="Density",
    )

    # axs[2].format(
    #     xlim=(0,80),xlabel="Precipitable Water / mm",
    #     ylabel="Frequency",
    # )

    f.savefig(plotsdir("tcwfreq-SEAx0.25.png"),transparent=false,dpi=200)

end

plotfreqhist()
