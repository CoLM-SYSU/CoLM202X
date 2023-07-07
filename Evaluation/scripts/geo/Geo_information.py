# -*- coding: utf-8 -*-
import os
import pandas as pd
import xarray as xr
import numpy as np
import sys
import shutil 
from joblib import Parallel, delayed

class Inundation:
    def __init__(self,namelist):
        self.name = 'Inundation_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['InundationLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['InundationLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['InundationLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['InundationLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['InundationLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['InundationLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['InundationLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['InundationLIST']['Sim_Eyear']

        self.obs_source              =  namelist['InundationLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['InundationLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['InundationLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['InundationLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['InundationLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['InundationLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['InundationLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['InundationLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Inundation/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='GIEMS_v2_2020'):
            if ((namelist['General']['compare_Tres']!="Month")):
                print('The compare_Tres should be "Month" for GIEMS')
                print('compare_Tres is set to be "Month"')
                namelist['General']['compare_Tres']="Month"
            if ((namelist['General']['compare_Gres']!="15min")):
                print('The compare_Gres should be "15min" for GIEMS')
                print('compare_Gres is set to be "15min"')
                namelist['General']['compare_Gres']="15min"

        elif (self.obs_source=='GSWO_2023'): 
            if (namelist['General']['compare_Tres']=="Hour")|(namelist['General']['compare_Tres']=="Day"):
                print(' the compare_Tres should be "Month" or longer ')
                sys.exit(1) 
        elif (self.obs_source=='GSWO_Pekel_2016'):
            if (namelist['General']['compare_Tres']=="Hour")|(namelist['General']['compare_Tres']=="Day"):
                print(' the compare_Tres should be "Month" or longer ')
                sys.exit(1) 
        elif (self.obs_source=='Pekel'):
            if (namelist['General']['compare_Tres']=="Hour")|(namelist['General']['compare_Tres']=="Day"):
                print(' the compare_Tres should be "Month" or longer ')
                sys.exit(1) 
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['InundationLIST']['OBSDIR']      , namelist['InundationLIST']['figplot'], namelist['InundationLIST']['Obs_GeoRes']
        del  namelist['InundationLIST']['Sim_Dir']     , namelist['InundationLIST']['Sim_TimRes']  , namelist['InundationLIST']['Sim_DataGroupby']
        del  namelist['InundationLIST']['Sim_Suffix']  , namelist['InundationLIST']['Sim_Prefix']  , namelist['InundationLIST']['Sim_GeoRes'], namelist['InundationLIST']['Obs_source']
        del  namelist['InundationLIST']['Sim_Syear'], namelist['InundationLIST']['Sim_Eyear']  , namelist['InundationLIST']['Obs_Syear'], namelist['InundationLIST']['Obs_Eyear']
        del  namelist['InundationLIST']['Obs_DataGroupby'], namelist['InundationLIST']['Obs_Suffix'],namelist['InundationLIST']['Obs_Prefix']
        del  namelist['InundationLIST']['Obs_TimRes']
        self.variables               =  namelist['InundationLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: InundationLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

class Evapotranspiration:
    def __init__(self,namelist):
        self.name = 'Evapotranspiration_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['EvapotranspirationLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['EvapotranspirationLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['EvapotranspirationLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['EvapotranspirationLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['EvapotranspirationLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['EvapotranspirationLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['EvapotranspirationLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['EvapotranspirationLIST']['Sim_Eyear']

        self.obs_source              =  namelist['EvapotranspirationLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['EvapotranspirationLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['EvapotranspirationLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['EvapotranspirationLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['EvapotranspirationLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['EvapotranspirationLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['EvapotranspirationLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['EvapotranspirationLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Evapotranspiration/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='GLEAM'):
            print('GLEAM .')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['EvapotranspirationLIST']['OBSDIR']      , namelist['EvapotranspirationLIST']['figplot'], namelist['EvapotranspirationLIST']['Obs_GeoRes']
        del  namelist['EvapotranspirationLIST']['Sim_Dir']     , namelist['EvapotranspirationLIST']['Sim_TimRes']  , namelist['EvapotranspirationLIST']['Sim_DataGroupby']
        del  namelist['EvapotranspirationLIST']['Sim_Suffix']  , namelist['EvapotranspirationLIST']['Sim_Prefix']  , namelist['EvapotranspirationLIST']['Sim_GeoRes'], namelist['EvapotranspirationLIST']['Obs_source']
        del  namelist['EvapotranspirationLIST']['Sim_Syear'], namelist['EvapotranspirationLIST']['Sim_Eyear']  , namelist['EvapotranspirationLIST']['Obs_Syear'], namelist['EvapotranspirationLIST']['Obs_Eyear']
        del  namelist['EvapotranspirationLIST']['Obs_DataGroupby'], namelist['EvapotranspirationLIST']['Obs_Suffix'],namelist['EvapotranspirationLIST']['Obs_Prefix']
        del  namelist['EvapotranspirationLIST']['Obs_TimRes']
     

        self.variables               =  namelist['EvapotranspirationLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: EvapotranspirationLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

        if self.Sim_DataGroupby=='single':
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            print(self.casedir+'/scratch')
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
        
            num=len(ds['time'])
                
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="Y", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="Y", periods=num, calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'

class Transpiration:
    def __init__(self,namelist):
        self.name = 'Transpiration_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['TranspirationLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['TranspirationLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['TranspirationLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['TranspirationLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['TranspirationLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['TranspirationLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['TranspirationLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['TranspirationLIST']['Sim_Eyear']

        self.obs_source              =  namelist['TranspirationLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['TranspirationLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['TranspirationLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['TranspirationLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['TranspirationLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['TranspirationLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['TranspirationLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['TranspirationLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Transpiration/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='GLEAM'):
            print('GLEAM .')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['TranspirationLIST']['OBSDIR']      , namelist['TranspirationLIST']['figplot'], namelist['TranspirationLIST']['Obs_GeoRes']
        del  namelist['TranspirationLIST']['Sim_Dir']     , namelist['TranspirationLIST']['Sim_TimRes']  , namelist['TranspirationLIST']['Sim_DataGroupby']
        del  namelist['TranspirationLIST']['Sim_Suffix']  , namelist['TranspirationLIST']['Sim_Prefix']  , namelist['TranspirationLIST']['Sim_GeoRes'], namelist['TranspirationLIST']['Obs_source']
        del  namelist['TranspirationLIST']['Sim_Syear'], namelist['TranspirationLIST']['Sim_Eyear']  , namelist['TranspirationLIST']['Obs_Syear'], namelist['TranspirationLIST']['Obs_Eyear']
        del  namelist['TranspirationLIST']['Obs_DataGroupby'], namelist['TranspirationLIST']['Obs_Suffix'],namelist['TranspirationLIST']['Obs_Prefix']
        del  namelist['TranspirationLIST']['Obs_TimRes']
                    
        self.variables               =  namelist['TranspirationLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: TranspirationLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

        if self.Sim_DataGroupby=='single':
            shutil.rmtree(self.casedir+'/scratch',ignore_errors=True)
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
            num=len(ds['time'])
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="Y", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="Y", periods=num, calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'


class Interception:
    def __init__(self,namelist):
        self.name = 'Interception_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['InterceptionLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['InterceptionLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['InterceptionLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['InterceptionLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['InterceptionLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['InterceptionLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['InterceptionLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['InterceptionLIST']['Sim_Eyear']

        self.obs_source              =  namelist['InterceptionLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['InterceptionLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['InterceptionLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['InterceptionLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['InterceptionLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['InterceptionLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['InterceptionLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['InterceptionLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Interception/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='GLEAM'):
            print('GLEAM .')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['InterceptionLIST']['OBSDIR']      , namelist['InterceptionLIST']['figplot'], namelist['InterceptionLIST']['Obs_GeoRes']
        del  namelist['InterceptionLIST']['Sim_Dir']     , namelist['InterceptionLIST']['Sim_TimRes']  , namelist['InterceptionLIST']['Sim_DataGroupby']
        del  namelist['InterceptionLIST']['Sim_Suffix']  , namelist['InterceptionLIST']['Sim_Prefix']  , namelist['InterceptionLIST']['Sim_GeoRes'], namelist['InterceptionLIST']['Obs_source']
        del  namelist['InterceptionLIST']['Sim_Syear'], namelist['InterceptionLIST']['Sim_Eyear']  , namelist['InterceptionLIST']['Obs_Syear'], namelist['InterceptionLIST']['Obs_Eyear']
        del  namelist['InterceptionLIST']['Obs_DataGroupby'], namelist['InterceptionLIST']['Obs_Suffix'],namelist['InterceptionLIST']['Obs_Prefix']
        del  namelist['InterceptionLIST']['Obs_TimRes']
                    
        self.variables               =  namelist['InterceptionLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: InterceptionLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

        if self.Sim_DataGroupby=='single':
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
        
            num=len(ds['time'])
                
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar='standard')
                else:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'

class SoilEvaporation:
    def __init__(self,namelist):
        self.name = 'SoilEvaporation_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['SoilEvaporationLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['SoilEvaporationLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['SoilEvaporationLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['SoilEvaporationLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['SoilEvaporationLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['SoilEvaporationLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['SoilEvaporationLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['SoilEvaporationLIST']['Sim_Eyear']

        self.obs_source              =  namelist['SoilEvaporationLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['SoilEvaporationLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['SoilEvaporationLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['SoilEvaporationLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['SoilEvaporationLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['SoilEvaporationLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['SoilEvaporationLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['SoilEvaporationLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/SoilEvaporation/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='GLEAM'):
            print('GLEAM.')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['SoilEvaporationLIST']['OBSDIR']      , namelist['SoilEvaporationLIST']['figplot'], namelist['SoilEvaporationLIST']['Obs_GeoRes']
        del  namelist['SoilEvaporationLIST']['Sim_Dir']     , namelist['SoilEvaporationLIST']['Sim_TimRes']  , namelist['SoilEvaporationLIST']['Sim_DataGroupby']
        del  namelist['SoilEvaporationLIST']['Sim_Suffix']  , namelist['SoilEvaporationLIST']['Sim_Prefix']  , namelist['SoilEvaporationLIST']['Sim_GeoRes'], namelist['SoilEvaporationLIST']['Obs_source']
        del  namelist['SoilEvaporationLIST']['Sim_Syear'], namelist['SoilEvaporationLIST']['Sim_Eyear']  , namelist['SoilEvaporationLIST']['Obs_Syear'], namelist['SoilEvaporationLIST']['Obs_Eyear']
        del  namelist['SoilEvaporationLIST']['Obs_DataGroupby'], namelist['SoilEvaporationLIST']['Obs_Suffix'],namelist['SoilEvaporationLIST']['Obs_Prefix']
        del  namelist['SoilEvaporationLIST']['Obs_TimRes']
                    
        self.variables               =  namelist['SoilEvaporationLIST']
        print(self.variables)


        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: SoilEvaporationLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

        if self.Sim_DataGroupby=='single':
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
        
            num=len(ds['time'])
                
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar='standard')
                else:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'



class SoilMoisture:
    def __init__(self,namelist):
        self.name = 'SoilMoisture_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['SoilMoistureLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['SoilMoistureLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['SoilMoistureLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['SoilMoistureLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['SoilMoistureLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['SoilMoistureLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['SoilMoistureLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['SoilMoistureLIST']['Sim_Eyear']

        self.obs_source              =  namelist['SoilMoistureLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['SoilMoistureLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['SoilMoistureLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['SoilMoistureLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['SoilMoistureLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['SoilMoistureLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['SoilMoistureLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['SoilMoistureLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/SoilMoisture/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='GLEAM'):
            print('GLEAM .')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['SoilMoistureLIST']['OBSDIR']      , namelist['SoilMoistureLIST']['figplot'], namelist['SoilMoistureLIST']['Obs_GeoRes']
        del  namelist['SoilMoistureLIST']['Sim_Dir']     , namelist['SoilMoistureLIST']['Sim_TimRes']  , namelist['SoilMoistureLIST']['Sim_DataGroupby']
        del  namelist['SoilMoistureLIST']['Sim_Suffix']  , namelist['SoilMoistureLIST']['Sim_Prefix']  , namelist['SoilMoistureLIST']['Sim_GeoRes'], namelist['SoilMoistureLIST']['Obs_source']
        del  namelist['SoilMoistureLIST']['Sim_Syear'], namelist['SoilMoistureLIST']['Sim_Eyear']  , namelist['SoilMoistureLIST']['Obs_Syear'], namelist['SoilMoistureLIST']['Obs_Eyear']
        del  namelist['SoilMoistureLIST']['Obs_DataGroupby'], namelist['SoilMoistureLIST']['Obs_Suffix'],namelist['SoilMoistureLIST']['Obs_Prefix']
        del  namelist['SoilMoistureLIST']['Obs_TimRes']
                    
        self.variables               =  namelist['SoilMoistureLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: SoilMoistureLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

        if self.Sim_DataGroupby=='single':
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
        
            num=len(ds['time'])
                
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar='standard')
                else:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'

class Runoff:
    def __init__(self,namelist):
        self.name = 'Runoff_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['RunoffLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['RunoffLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['RunoffLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['RunoffLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['RunoffLIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['RunoffLIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['RunoffLIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['RunoffLIST']['Sim_Eyear']

        self.obs_source              =  namelist['RunoffLIST']['Obs_source']
        self.Obs_TimRes              =  namelist['RunoffLIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['RunoffLIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['RunoffLIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['RunoffLIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['RunoffLIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['RunoffLIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['RunoffLIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Runoff/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='ERA5-Land'):
            print('ERA5-Land .')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['RunoffLIST']['OBSDIR']      , namelist['RunoffLIST']['figplot'], namelist['RunoffLIST']['Obs_GeoRes']
        del  namelist['RunoffLIST']['Sim_Dir']     , namelist['RunoffLIST']['Sim_TimRes']  , namelist['RunoffLIST']['Sim_DataGroupby']
        del  namelist['RunoffLIST']['Sim_Suffix']  , namelist['RunoffLIST']['Sim_Prefix']  , namelist['RunoffLIST']['Sim_GeoRes'], namelist['RunoffLIST']['Obs_source']
        del  namelist['RunoffLIST']['Sim_Syear'], namelist['RunoffLIST']['Sim_Eyear']  , namelist['RunoffLIST']['Obs_Syear'], namelist['RunoffLIST']['Obs_Eyear']
        del  namelist['RunoffLIST']['Obs_DataGroupby'], namelist['RunoffLIST']['Obs_Suffix'],namelist['RunoffLIST']['Obs_Prefix']
        del  namelist['RunoffLIST']['Obs_TimRes']
                    

        self.variables               =  namelist['RunoffLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: RunoffLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)


        if self.Sim_DataGroupby=='single':
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
        
            num=len(ds['time'])
                
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar='standard')
                else:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'

class LAI:
    def __init__(self,namelist):
        self.name = 'LAI_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_Tres            =  namelist['General']['compare_Tres']
        self.compare_Gres            =  namelist['General']['compare_Gres']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['LAILIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['LAILIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['LAILIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['LAILIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['LAILIST']['Sim_Prefix']
        self.Sim_GeoRes              =  namelist['LAILIST']['Sim_GeoRes']
        self.Sim_Syear               =  namelist['LAILIST']['Sim_Syear']
        self.Sim_Eyear               =  namelist['LAILIST']['Sim_Eyear']

        self.obs_source              =  namelist['LAILIST']['Obs_source']
        self.Obs_TimRes              =  namelist['LAILIST']['Obs_TimRes']
        self.Obs_GeoRes              =  namelist['LAILIST']['Obs_GeoRes']
        self.OBSDIR                  =  namelist['LAILIST']['OBSDIR'] 
        self.Obs_Syear               =  namelist['LAILIST']['Obs_Syear'] 
        self.Obs_Eyear               =  namelist['LAILIST']['Obs_Eyear'] 
        self.Obs_DataGroupby         =  namelist['LAILIST']['Obs_DataGroupby']
        self.figplot                 =  namelist['LAILIST']['figplot']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/LAI/{self.obs_source}/'
        self.use_Syear=max(self.Obs_Syear,self.Sim_Syear,self.Syear)
        self.use_Eyear=min(self.Obs_Eyear,self.Sim_Eyear,self.Eyear)
        self.num_cores               =  namelist['General']['num_cores']
        if (self.obs_source=='Yuan_etal'):
            print('Yuan_etal')
        else:
            print('The source is not exist.')
            sys.exit(1)

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['LAILIST']['OBSDIR']      , namelist['LAILIST']['figplot'], namelist['LAILIST']['Obs_GeoRes']
        del  namelist['LAILIST']['Sim_Dir']     , namelist['LAILIST']['Sim_TimRes']  , namelist['LAILIST']['Sim_DataGroupby']
        del  namelist['LAILIST']['Sim_Suffix']  , namelist['LAILIST']['Sim_Prefix']  , namelist['LAILIST']['Sim_GeoRes'], namelist['LAILIST']['Obs_source']
        del  namelist['LAILIST']['Sim_Syear'], namelist['LAILIST']['Sim_Eyear']  , namelist['LAILIST']['Obs_Syear'], namelist['LAILIST']['Obs_Eyear']
        del  namelist['LAILIST']['Obs_DataGroupby'], namelist['LAILIST']['Obs_Suffix'],namelist['LAILIST']['Obs_Prefix']
        del  namelist['LAILIST']['Obs_TimRes']
                     
        self.variables               =  namelist['LAILIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: LAILIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)
        if self.Sim_DataGroupby=='single':
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            # Open the netCDF file
            VarFile = os.path.join(self.Sim_Dir, f'{self.Sim_Suffix}{self.Sim_Prefix}.nc')
            ds = xr.open_dataset(VarFile)
        
            num=len(ds['time'])
                
            if (self.Sim_TimRes=="Hour"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="H", periods=num,calendar='standard')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="H", periods=num, calendar="noleap") 

            elif (self.Sim_TimRes=="Day"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="standard") 
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="D", periods=num, calendar="noleap") 
            elif (self.Sim_TimRes=="Month"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="standard") 
                    print('leap')
                else:
                    ds['time'] = xr.cftime_range(start=f"{self.Sim_Syear}-01-01", freq="M", periods=num, calendar="noleap") 
                    print('noleap')
            elif (self.Sim_TimRes=="Year"):
                if any(ds['time'].dt.dayofyear) == 366:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar='standard')
                else:
                    ds['time'] = pd.date_range(f"{self.Sim_Syear}-01-01", freq="Y", periods=num,calendar="noleap") 
            else:
                sys.exit(1)
           
           
            # Split the data into yearly files
            for year in range(self.use_Syear, self.use_Eyear+1):
                ds_year = ds.sel(time=slice(f'{year}-01-01T00:00:00',f'{year}-12-31T23:59:59'))[list(self.variables.values())]
                # Remove all attributes
                ds_year.attrs = {}
                ds_year.to_netcdf(os.path.join(self.casedir,'scratch',f'{self.Sim_Suffix}{year}{self.Sim_Prefix}.nc'))
            ds.close()
        self.Sim_Dir=self.casedir+'/scratch'