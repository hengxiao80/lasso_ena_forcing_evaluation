import pandas as pd

# 2016-10-22: Cu-under-St, CF > 0.5
# 2018-10-28: Overcast (except near the bigger islands), CF ~ 1
# 2018-11-21: Open cell Cu

periods = [
    ['20161021', '20161023'],
    ['20181027', '20181029'],
    ['20181120', '20181122'],
]
for d1, d2 in periods:
    with open(f'dl_{d1}_{d2}.txt', 'w') as f:
        dates = pd.date_range(start=d1, end=d2)
        for adate in dates:
            date_string = adate.strftime('%Y%m%d') 
            dir_string = adate.strftime('%Y/%m') 
            f.write(f'https://data.gesdisc.earthdata.nasa.gov/data/MERRA2/M2I3NVASM.5.12.4/{dir_string}/MERRA2_400.inst3_3d_asm_Nv.{date_string}.nc4\n')
            f.write(f'https://data.gesdisc.earthdata.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{dir_string}/MERRA2_400.tavg1_2d_flx_Nx.{date_string}.nc4\n')
            f.write(f'https://data.gesdisc.earthdata.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4/{dir_string}/MERRA2_400.tavg1_2d_slv_Nx.{date_string}.nc4\n')