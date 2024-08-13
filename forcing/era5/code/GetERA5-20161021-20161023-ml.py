# Template script for downloading ERA5 model-level data for WRF/WPS.
# From https://dreambooker.site/2018/04/20/Initializing-the-WRF-model-with-ERA5 on 6-Mar-2020.

import cdsapi
c = cdsapi.Client()
c.retrieve('reanalysis-era5-complete',{
    'class':'ea',
    'date':'20161021/to/20161023',
    'area':'46.5/-35.5/31.5/-20.5',
    'expver':'1',
    'levelist': '1/to/137',
    'levtype':'ml',
    'param':'129/130/131/132/133/152/75/76/135/246/247/248',
    'stream':'oper',
    'time':'00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
    'type':'an',
    'grid':"0.25/0.25",
},'ERA5-20161021-20161023-ml.grb')

