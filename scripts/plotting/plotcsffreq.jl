using DrWatson
@quickactivate "MaritimeContPWV"
using Dates
using NCDatasets

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

include(srcdir("common.jl"))

function plotfreqhist(
    coast::Real=0.5, thr::Real=0.05
)

    ds = NCDataset(datadir("compiled/csffreqsave-SEAx0.25.nc"))
    freqall = ds["bin_frq"][:]; csf = ds["csf"][:]; close(ds)
    csf = csf[1:(end-1)] .+ (csf[2]-csf[1])/2; nbins = length(csf)

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
    lsmsea = lsm .< coast; lsmlnd = lsm .> coast;

    freqsea = dropdims(sum(freqall .* lsmsea,dims=(1,2)),dims=(1,2));
    freqlnd = dropdims(sum(freqall .* lsmlnd,dims=(1,2)),dims=(1,2));
    freqall = dropdims(sum(freqall,dims=(1,2)),dims=(1,2));

    pplt.close(); f,axs = pplt.subplots(nrows=1,aspect=2,axwidth=4);

    axs[1].plot(csf,freqall / sum(freqall) * nbins,lw=1,label="All",legend="ul",c="k")
    axs[1].plot(csf,freqsea / sum(freqsea) * nbins,lw=1,label="Sea",legend="ul",c="b")
    axs[1].plot(csf,freqlnd / sum(freqlnd) * nbins,lw=1,label="Land",legend="ul",c="g")
    axs[1].plot([0,nbins],[1,1]*thr,lw=0.5,c="k")

    # axs[1].plot(csf,freqall,lw=1,label="All",legend="lr",c="k")
    # axs[1].plot(csf,freqsea,lw=1,label="Sea",legend="lr",c="b")
    # axs[1].plot(csf,freqlnd,lw=1,label="Land",legend="lr",c="g")

    axs[1].format(
        xlim=(0,nbins),xlabel="Column Saturation Fraction",
        ylim=(0,3),ylabel="Density"
    )

    # axs[2].format(
    #     xlim=(0,1),xlabel="Column Saturation Fraction",
    #     ylabel="Frequency",
    # )

    f.savefig(plotsdir("csffreq-SEAx0.25.png"),transparent=false,dpi=200)

end

function plotfreqhistmonthly(
    coast::Real=0.5, thr::Real=0.05
)

    ds = NCDataset(datadir("compiled/csffreqsave-SEAx0.25.nc"))
    nlon = ds.dim["longitude"]
    nlat = ds.dim["latitude"]
    csf  = ds["csf"][:]; nbins = length(csf) - 1
    close(ds)

    frq = zeros(nlon,nlat,nbins,12)
    csf = csf[1:(end-1)] .+ (csf[2]-csf[1])/2

    for yr = 1980 : 2019, mo = 1 : 12

        @info "$(Dates.now()) - Extracting binned frequencies for $yr $(monthname(mo)) ..."
        dt = Date(yr,mo)
        fol = datadir("compiled/csffreq/$yr")
        ds = NCDataset("$fol/csffreqsave-SEAx0.25-$(yrmo2str(dt)).nc")
        frq[:,:,:,mo] += ds["bin_frq"][:]*1
        close(ds)

    end

    frq = frq / 40

    mds = Dataset(datadir("era5-SEAx0.25-lsm.nc")); lsm = mds["lsm"][:]*1; close(mds)
    lsmsea = lsm .< coast; lsmlnd = lsm .> coast;

    for mo = 1 : 12

        freqsea = dropdims(sum(frq[:,:,:,mo] .* lsmsea,dims=(1,2)),dims=(1,2));
        freqlnd = dropdims(sum(frq[:,:,:,mo] .* lsmlnd,dims=(1,2)),dims=(1,2));
        freqall = dropdims(sum(frq[:,:,:,mo],dims=(1,2)),dims=(1,2));

        pplt.close(); f,axs = pplt.subplots(nrows=1,aspect=2,axwidth=4);

        axs[1].plot(csf,freqall / sum(freqall) * nbins,lw=1,label="All",legend="ul",c="k")
        axs[1].plot(csf,freqsea / sum(freqsea) * nbins,lw=1,label="Sea",legend="ul",c="b")
        axs[1].plot(csf,freqlnd / sum(freqlnd) * nbins,lw=1,label="Land",legend="ul",c="g")
        axs[1].plot([0,nbins],[1,1]*thr,lw=0.5,c="k")

        # axs[1].plot(csf,freqall,lw=1,label="All",legend="lr",c="k")
        # axs[1].plot(csf,freqsea,lw=1,label="Sea",legend="lr",c="b")
        # axs[1].plot(csf,freqlnd,lw=1,label="Land",legend="lr",c="g")

        axs[1].format(
            xlim=(0,nbins),xlabel="Column Saturation Fraction",
            ylim=(0,3),ylabel="Density"
        )

        # axs[2].format(
        #     xlim=(0,1),xlabel="Column Saturation Fraction",
        #     ylabel="Frequency",
        # )

        f.savefig(plotsdir("csffreq-SEAx0.25-month$mo.png"),transparent=false,dpi=200)

    end

end

plotfreqhist()
plotfreqhistmonthly()
