# -*- coding: utf-8 -*-
import os
import pandas as pd
import xarray as xr
import numpy as np
import sys
import shutil 
from joblib import Parallel, delayed

def remove_unused_variables(modname, namelist):
    del namelist[modname+'LIST']['Obs_Dir']
    del namelist[modname+'LIST']['figplot']
    del namelist[modname+'LIST']['Obs_GeoRes']
    del namelist[modname+'LIST']['Sim_Dir']
    del namelist[modname+'LIST']['Sim_TimRes']
    del namelist[modname+'LIST']['Sim_DataGroupby']
    del namelist[modname+'LIST']['Sim_Suffix']
    del namelist[modname+'LIST']['Sim_Prefix']
    del namelist[modname+'LIST']['Sim_GeoRes']
    del namelist[modname+'LIST']['Obs_source']
    del namelist[modname+'LIST']['Sim_Syear']
    del namelist[modname+'LIST']['Sim_Eyear']
    del namelist[modname+'LIST']['Obs_Syear']
    del namelist[modname+'LIST']['Obs_Eyear']
    del namelist[modname+'LIST']['Obs_DataGroupby']
    del namelist[modname+'LIST']['Obs_Suffix']
    del namelist[modname+'LIST']['Obs_Prefix']
    del namelist[modname+'LIST']['Obs_TimRes']
    variables = namelist[modname+'LIST']
    variables = {k: v for k, v in variables.items() if isinstance(v, str)}
    if len(variables) == 0:
        print(f"Error: {modname}LIST is empty!")
        sys.exit()
    return variables

def select_metrics(namelist):
    select_metrics = {k: v for k, v in namelist['metrics'].items() if v}
    return select_metrics
#                    
def split_year(info,TimRes,vars,Dir,Suffix,Prefix):
    os.makedirs(info.casedir+'/scratch', exist_ok=True)
    print(info.casedir+'/scratch')
    # Open the netCDF file
    VarFile = os.path.join(Dir, f'{Suffix}{Prefix}.nc')
    ds = xr.open_dataset(VarFile)
    num=len(ds['time'])
                
    if (TimRes=="Hour"):
        if any(ds['time'].dt.dayofyear) == 366:
            ds['time'] = pd.date_range(f"{info.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
        else:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

    elif (TimRes=="Day"):
        if any(ds['time'].dt.dayofyear) == 366:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
        else:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
    elif (TimRes=="Month"):
        if any(ds['time'].dt.dayofyear) == 366:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
            print('leap')
        else:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
            print('noleap')
    elif (TimRes=="Year"):
        if any(ds['time'].dt.dayofyear) == 366:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="Y", periods=num, calendar="standard") 
            print('leap')
        else:
            ds['time'] = xr.cftime_range(start=f"{info.Sim_Syear}-01-01", freq="Y", periods=num, calendar="noleap") 
    else:
        sys.exit(1)
           
           
    # Split the data into yearly files
    for year in range(info.use_Syear, info.use_Eyear+1):
        ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[vars]
        # Remove all attributes
        ds_year.attrs = {}
        ds_year.to_netcdf(os.path.join(info.casedir,'scratch',f'{Suffix}{year}{Prefix}.nc'))
    ds.close()

class get_general_info:
    def __init__(self,vanme,namelist):
        self.name = self.__class__.__name__  
        self.Minimum_lenghth         =  (namelist['General']['Min_year'])
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']
        self.Sim_Dir                 =  namelist[vanme+'LIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist[vanme+'LIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist[vanme+'LIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist[vanme+'LIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist[vanme+'LIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist[vanme+'LIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist[vanme+'LIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist[vanme+'LIST']['Sim_Eyear']

        self.obs_source              =  namelist[vanme+'LIST']['Obs_source']
        self.Obs_TimRes              =  namelist[vanme+'LIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist[vanme+'LIST']['Obs_GeoRes']
        self.Obs_Dir                  =  namelist[vanme+'LIST']['Obs_Dir'] 
        self.Obs_Syear               =  namelist[vanme+'LIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist[vanme+'LIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist[vanme+'LIST']['Obs_DataGroupby']
        self.figplot                 =  namelist[vanme+'LIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/{vanme}/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']

        # remove_unused_variables
        self.variables=remove_unused_variables(vanme,namelist)
        self.metrics = select_metrics(namelist)

        if self.Sim_DataGroupby=='single':
            k1=list(self.variables.values())
            split_year(self,self.Sim_TimRes, k1, self.Sim_Dir,self.Sim_Suffix,self.Sim_Prefix)
            self.Sim_Dir=self.casedir+'/scratch'
            
        if self.Obs_DataGroupby=='single':
            k1=list(self.variables.keys())
            split_year(self,self.Obs_TimRes, k1, self.Obs_Dir,self.Obs_Suffix,self.Obs_Prefix)
            self.Obs_Dir=self.casedir+'/scratch'

