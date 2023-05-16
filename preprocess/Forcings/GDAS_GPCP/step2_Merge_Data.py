# -*- coding: utf-8 -*-
"""
A libray with Python functions for merging GLDAS data / need to be updated
****
"""
__author__ = "Zhongwang Wei / zhongwang007@gmail.com"
__version__ = "0.1"
__release__ = "0.1"
__date__ = "Oct 2020"


import xarray as xr
from pylab import *


combined=xr.open_mfdataset('GLDAS_*.nc4', combine='nested',concat_dim="time",chunks={'time': 365})#.chunk(chunks='auto', token=None, lock=False)
#combined.to_netcdf('test.nc',engine='netcdf4',encoding={'time':{'units':'days since 2000-01-01 00:00:00'}})
comp = dict(zlib=True, complevel=6, _FillValue= -9999) #,chunksizes=(1,150,360)
encoding = {var: comp for var in combined.data_vars}
#combined.to_netcdf('test.nc',engine='netcdf4',encoding=encoding)
combined.to_netcdf('test.nc',engine='netcdf4',encoding=encoding) #scipy
#combined.to_netcdf('test.nc',engine='netcdf4')