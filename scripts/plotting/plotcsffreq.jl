using DrWatson
@quickactivate "MaritimeContPWV"
using NCDatasets

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plotfreqhist(
    thr::Real=0.5
)

    ds = NCDataset(datadir("compiled/csffreqsave-SEAx0.25.nc"))
    freqall = ds["bin_frq"][:]; csf = ds["csf"][:]; close(ds)
    csf = csf[1:(end-1)] .+ (csf[2]-csf[1])/2; nbins = length(csf)

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
    lsmsea = lsm .< thr; lsmlnd = lsm .> thr;

    freqsea = dropdims(sum(freqall .* lsmsea,dims=(1,2)),dims=(1,2));
    freqlnd = dropdims(sum(freqall .* lsmlnd,dims=(1,2)),dims=(1,2));
    freqall = dropdims(sum(freqall,dims=(1,2)),dims=(1,2));

    pplt.close(); f,axs = pplt.subplots(nrows=1,aspect=2,axwidth=4);

    axs[1].plot(csf,freqall / sum(freqall) * nbins,lw=1,label="All",legend="lr",c="k")
    axs[1].plot(csf,freqsea / sum(freqsea) * nbins,lw=1,label="Sea",legend="lr",c="b")
    axs[1].plot(csf,freqlnd / sum(freqlnd) * nbins,lw=1,label="Land",legend="lr",c="g")

    # axs[2].plot(csf,freqall,lw=1,label="All",legend="lr",c="k")
    # axs[2].plot(csf,freqsea,lw=1,label="Sea",legend="lr",c="b")
    # axs[2].plot(csf,freqlnd,lw=1,label="Land",legend="lr",c="g")

    axs[1].format(
        xlim=(0,1),xlabel="Column Saturation Fraction",
        ylim=(0,5),ylabel="Density",
    )

    # axs[2].format(
    #     xlim=(0,1),xlabel="Column Saturation Fraction",
    #     ylabel="Frequency",
    # )

    f.savefig(plotsdir("csffreq-SEAx0.25.png"),transparent=false,dpi=200)

end

plotfreqhist()
