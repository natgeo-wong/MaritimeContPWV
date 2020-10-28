# **<div align="center">MaritimeContPWV: A comparison of Satellite and Reanalysis Precipitable Water products over the Maritime Continent</div>**

This repository contains the analysis scripts and output for the **MaritimeContPWV** project, written in Julia.  We aim to investigate if ERA5 total column water data can be used as a good approximation to satellite and GNSS observations of precipitable water over the Maritime Continent.

**Created/Mantained By:** Nathanael Wong (nathanaelwong@fas.harvard.edu)\
**Other Contributors:** Lujia Feng (lfeng@ntu.edu.sg)

> In this project I compare the MIMIC-TPW2m satellite and ERA5 reanalysis precipitable water prod- ucts using GNSS-derived precipitable water as the ground truth. Despite the fact that ERA5 is the product of reanalysis model output, it is closer to GNSS-derived precipitable water values than MIMIC-TPW2m. We postulate that this due to the fact that MIMIC-TPW2m algorithms do not work well over land as compared to the ocean, and therefore reanalysis datasets are still valuable in helping us characterize the climatology of precipitable water vapour over the Maritime Continent. Future studies using precipitable water data should keep this in mind.

## Progress

* [x] Downloaded ERA5 and ERA-Interim Raw Data from CDS
* [x] Downloaded GPM and TRMM Raw Data from PMM website
* [x] Retrieved GNSS Raw Data from EOS

* [x] Compare between different Precipitable Water datasets
  - [x] Compare ERA5 and ERA-Interim TCWV with TCW
  - [x] Over land, compare GNSS PWV with ERA5 and ERA-Interim TCWV
  - [x] Over ocean, compare ERA5 and ERA-Interim TCWV with MIMIC-TPW2m
  
* [ ] Produce TPW-Precipitation curve over land, ocean and coastline
  - [x] General overview curve
  - [ ] Compare curve in different times of year / season
  - [ ] Compare against MJO-Phase/Amplitude data

## Installation

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> ] activate .
    Activating environment at `~/Projects/MaritimeContPWV/Project.toml`

   (MaritimeContPWV) pkg> instantiate
   (MaritimeContPWV) pkg> add GeoRegions#master
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

*(Note: You need to install the #master versions of GeoRegions.jl as of now.)*

## **Other Acknowledgements**
> Project Repository Template generated using [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) created by George Datseris.
> Work from this project was funded by the [Earth Observatory of Singapore](https://earthobservatory.sg/).
