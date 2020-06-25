# **<div align="center">MaritimeContPWV: A comparison of Satellite and Reanalysis Precipitable Water products over the Maritime Continent</div>**

This repository contains the basic setup for the **MaritimeContPWV** project, used to xxx

**Created/Mantained By:** Nathanael Wong (nathanaelwong@fas.harvard.edu)\
**Other Contributors:** Lujia Feng (lfeng@ntu.edu.sg)

> Introductory Text Here.

## Current Status

**Data Downloads**
* [x] Downloaded ERA5 and ERA-Interim Raw Data from CDS
* [x] Downloaded GPM and TRMM Raw Data from PMM website
* [x] Retrieved GNSS Raw Data from EOS

**Analysis/Comparison between Datasets**

* [ ] Compare between different Precipitable Water datasets
  - [ ] Compare ERA5 and ERA-Interim TCWV with TCW
  - [ ] Over land, compare GNSS PWV with ERA5 and ERA-Interim TCWV
  - [ ] Over ocean, compare ERA5 and ERA-Interim TCWV with MIMIC-TPW2m
* [ ] Produce TPW-Precipitation curve over land, ocean and coastline
  - [ ] General overview curve
  - [ ] Compare curve in different times of year / season
  - [ ] Compare against MJO-Phase/Amplitude data

## Project Setup

In order to obtain the relevant data, please follow the instructions listed below:

### 1) Required Julia Dependencies

The entire codebase is written in Julia.  If the data files are downloaded, you should be able to produce my results in their entirety.  The following are the most important Julia packages that were used in this project:
* Satellite Data Handling: `ClimateSatellite.jl`
* ECWMF Data Handling: `ClimateERA.jl`
* GNSS ZWD Data Handling: `ClimateGNSS.jl`
* NetCDF Data Handling: `NCDatasets.jl`

In order to reproduce the results, first you have to clone the repository, and instantiate the project environment in the Julia REPL in order to install the required packages:
```
git clone https://github.com/natgeo-wong/PiPWV.git
] activate .
instantiate
```

## **Other Acknowledgements**
> Project Repository Template generated using [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) created by George Datseris.
> Work from this project was funded by the [Earth Observatory of Singapore](https://earthobservatory.sg/).
