
GetERA5-20161021-20161023-ml.py --- for downloading model level data
GetERA5-20161021-20161023-sfc.py --- for downloading sfc data

You need .cdsapirc in your home dir, and cdsapi package installed (through pip)

The get-era5.py code automates the downloading a little bit. --- 12/06/2023

To convert the .grb file to .nc, simply do `grib_to_netcdf -o .nc .grb`