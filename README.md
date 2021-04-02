# A simple way to evaluate C3S seasonal forecasts monthly data

This repository contains a set of Python functions, contained in the folder `c3sdl`, to:

1. Retrieve seasonal monthly data from the Copernicus Climate Data Store (CDS)
2. Retrieve ERA5 monthly data
3. Regrid ERA5 to the forecast's grid
4. Calculate deterministic and probabilistic skills
5. Save the validation data in a NetCDF
6. Quickly create plots with the computed metrics

## Limitations
This code has been created to implement a quick workflow for the users that want to compute the skills in a quick way. Currently, it lacks the following features:

- Accessing seasonal and reanalysis data from local files
- Some well-known metrics (e.g. CRPSS)

However, both the features can be easily implemented and the classes have been designed considering this extension. 

Feel free to extend this code and contribute. 

## Requirements
This code doesn't run on Windows (because `xesmf` supports only Unix-based systems).
The packages needed to run this code are contained in the `environment.yml` file, if you are an Anaconda user you can type: 
```
conda env create -f environment.yml
```

The user also needs an account to the [Copernicus Climate Change CDS](https://cds.climate.copernicus.eu/) and the CDS API installed ([see here](https://cds.climate.copernicus.eu/api-how-to)) 

If today we can implement this workflow so easily, we have to thanks the following packages and their developers:
- [cfgrib](https://github.com/ecmwf/cfgrib): the Python interface to GRIB files provided by ECMWF
- [xesmf](https://xesmf.readthedocs.io/en/latest/): to regrid xarray data 
- [xskillscore](https://xskillscore.readthedocs.io/en/stable/): for the skill score calculation
- xarray, pandas, numpy: the open projects that we researchers should try [to support](https://numfocus.org/donate) as much as possible

## Examples

Let's see how to calculate the skill of ECMWF C3S seasonal forecast for JAS with a lead-time of one month. The first step is downloading the data:
```
import c3sdl
d = c3sdl.C3S_data_loader('ecmwf', '2m_temperature', 6, [2,3,4])
```
Unless we specify another folder in the `config.yml` file, a NetCDF named `ecmwf-2m_temperature-S6-L2-4.nc` will appear in the local directory with this structure:

```
<xarray.Dataset>
Dimensions:  (lat: 181, lon: 360, number: 25, year: 23)
Coordinates:
  * lon      (lon) float64 0.0 1.0 2.0 3.0 4.0 ... 355.0 356.0 357.0 358.0 359.0
  * lat      (lat) float64 90.0 89.0 88.0 87.0 86.0 ... -87.0 -88.0 -89.0 -90.0
  * number   (number) int32 0 1 2 3 4 5 6 7 8 9 ... 16 17 18 19 20 21 22 23 24
  * year     (year) int32 1993 1994 1995 1996 1997 ... 2011 2012 2013 2014 2015
    surface  int32 ...
    step     timedelta64[ns] ...
Data variables:
    obs_var  (year, lat, lon) float64 ...
    fct_var  (year, number, lat, lon) float32 ...
```
The two data variables contain respectively the JAS average of reanalysis data and the forecast ensemble (25 members).
We can evaluate the deterministic and probabilistic skills:
```
e = c3sdl.evaluator(d)
e.compute_deterministic()
e.compute_probabilistic()
```
It takes a couple of minutes and then, as for the retrieval, a NetCDF with the computer metrics will appear in the local directory. The NetCDF and the object returned by `e.get_skill_data()` contain all the computed metrics and the needed information for maps and charts.
However, a quick way to generate the maps is the following:
```
e.save_maps()
```
Then some PNGs will appear showing the global metrics and skill scores.

By default, the functions will calculate the metrics on the entire domain for all the hindcast years. However, using `e.get_data()` and `e.set_data()`, it is possible to get the forecast/reanalysis data, subset it (in space and/or time) and then go ahead with the computation. 

