## forcing_extraction.ipynb

- driver for generating forcings (ERA5 and MERRA2)

## extract_forcing.py

- only read in 3 days of ERA5/MERRA2 data
- contains the code to actually do the job 
- uses ARM sounding as initial condition
- can choose the initial time (`hrs_to_shift`, within a case day) and the intended duration of the simulation (`hrs_to_simulate`)
- can choose an averaging box (`dx`, half width of the box)
- the model top (ztop) and `dz` are hardwired as of now