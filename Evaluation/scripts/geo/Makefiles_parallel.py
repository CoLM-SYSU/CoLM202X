# -*- coding: utf-8 -*-
import os
import pandas as pd
import xarray as xr
import numpy as np
import sys
import shutil 
from dask.diagnostics import ProgressBar
from joblib import Parallel, delayed
os.environ['PYTHONWARNINGS']='ignore::FutureWarning'

# Check the platform
if not sys.platform.startswith('win'):
    import xesmf as xe

class Makefiles_parallel:
    def __init__(self, info):
        self.name = 'Makefile'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  info.Minimum_lenghth
        self.Max_lat                 =  info.Max_lat
        self.Min_lat                 =  info.Min_lat
        self.Max_lon                 =  info.Max_lon
        self.Min_lon                 =  info.Min_lon
        self.Syear                   =  info.Syear
        self.Eyear                   =  info.Eyear
        self.compare_Tres            =  info.compare_Tres
        self.compare_Gres            =  info.compare_Gres

        self.casename                =  info.casename

        self.Sim_Dir                 =  info.Sim_Dir
        self.Sim_TimRes              =  info.Sim_TimRes
        self.Sim_DataGroupby         =  info.Sim_DataGroupby
        self.Sim_Suffix              =  info.Sim_Suffix
        self.Sim_Prefix              =  info.Sim_Prefix
        self.Sim_GeoRes              =  info.Sim_GeoRes
        self.Sim_Syear               =  info.Sim_Syear
        self.Sim_Eyear               =  info.Sim_Eyear

        self.obs_source              =  info.obs_source
        self.Obs_GeoRes              =  info.Obs_GeoRes
        self.Obs_TimRes              =  info.Obs_TimRes
        self.Obs_Dir                  =  info.Obs_Dir
        self.Obs_Syear               =  info.Obs_Syear
        self.Obs_Eyear               =  info.Obs_Eyear
        self.Obs_DataGroupby         =  info.Obs_DataGroupby
        self.Obs_Suffix              =  info.Obs_Suffix
        self.Obs_Prefix              =  info.Obs_Prefix
        self.figplot                 =  info.figplot
        self.casedir                 =  info.casedir
        self.variables               =  info.variables
        self.metrics                 =  info.metrics
        self.use_Syear               =  info.use_Syear
        self.use_Eyear               =  info.use_Eyear
        self.num_cores               =  info.num_cores
        self.res                     = int(self.compare_Gres[:-3])/60. #degree

       # shutil.rmtree(self.casedir+'/scratch',ignore_errors=True)
       # os.makedirs(self.casedir+'/scratch', exist_ok=True)

    def make_obs_combine_parallel(self,ii):
        for key in self.variables.keys():
            if self.Obs_DataGroupby == 'Year':
                if (self.obs_source=='GIEMS_v2_2020') :
                    VarFiles = os.path.join(self.Obs_Dir, f'GIEMS_v2_2020_{key}_15min_{ii}.nc')
                elif (self.obs_source=='GLEAM'):
                    VarFiles = os.path.join(self.Obs_Dir, 'daily', f'{key}_{ii}_GLEAM_v3.7a.nc')
                elif (self.obs_source=='GLDAS'):
                    print('not ready yet')
                    sys.exit(1)
                elif ((self.obs_source=='ERA5-Land')and(key=='ro')):
                    VarFiles = os.path.join(self.Obs_Dir, f'ERA5LAND_runoff_{ii}.nc4') ###check here_2018_GLEAM_v3.7a.nc
                elif ((self.obs_source=='Yuan_etal')and(key=='lai')):
                    VarFiles = os.path.join(self.Obs_Dir, f'lai_8-day_30s_{ii}.nc4') ###check here_2018_GLEAM_v3.7a.nc
                else:
                    VarFiles = os.path.join(self.Obs_Dir, f'{self.Obs_Suffix}{ii}{self.Obs_Prefix}.nc')
                obsx0 = xr.open_dataset(VarFiles)
            else:
                print ('The Obs_DataGroupby is not Year-->combine it to Year')
                VarFiles = os.path.join(self.Obs_Dir, f'{self.Obs_Suffix}{ii}*{self.Obs_Prefix}.nc')
                obsx0= xr.open_mfdataset(VarFiles, combine='nested',concat_dim="time",decode_times=False,chunks={'time': 30},
                                   preprocess=lambda obsx: obsx[f'{key}'].astype('float32'))
            obsx0['lon'] = xr.where(obsx0.lon > 180, obsx0.lon - 360, obsx0.lon)

            lon_new = xr.DataArray(
                data=np.arange(self.Min_lon+self.res/2, self.Max_lon, self.res),
                dims=('lon',),
                coords={'lon': np.arange(self.Min_lon+self.res/2, self.Max_lon, self.res)},
                attrs={'units': 'degrees_east', 'long_name': 'longitude'}
                )
            lat_new = xr.DataArray(
                data=np.arange(self.Min_lat+self.res/2, self.Max_lat, self.res),
                dims=('lat',),
                coords={'lat': np.arange(self.Min_lat+self.res/2, self.Max_lat, self.res)},
                attrs={'units': 'degrees_north', 'long_name': 'latitude'}
                )
            new_grid = xr.Dataset({'lon': lon_new, 'lat': lat_new})

            if not sys.platform.startswith('win'):
                # Get the data variable
                #obsx = obsx[f'{key}'] #['fldfrc']
                # Create the regridder
                # Define the path to the weights file
                if (self.Min_lon==-180 and self.Max_lon==180 and self.Min_lat==-90 and self.Max_lat==90):
                    #if ii==self.use_Syear:
                    #    regridder = xe.Regridder(obsx0, new_grid, 'bilinear', periodic=True)
                    #else:
                    #    regridder = xe.Regridder(obsx0, new_grid, 'bilinear', periodic=True, reuse_weights=True)
                    regridder = xe.Regridder(obsx0, new_grid, 'bilinear', periodic=True)
                else:
                    #if ii==self.use_Syear:
                    #    regridder = xe.Regridder(obsx0, new_grid, 'bilinear', periodic=False)
                    #else:
                    #    regridder = xe.Regridder(obsx0, new_grid, 'bilinear', periodic=False, reuse_weights=True)
                    regridder = xe.Regridder(obsx0, new_grid, 'bilinear', periodic=False)

                # Perform the remapping
                obsx = regridder(obsx0)
            else:
                obsx = obsx[f'{key}'] #['fldfrc']
                obsx = obsx.interp(coords=new_grid.coords) 
            
            #need to improve here
            num=len(obsx['time'])
            if (self.Obs_TimRes=="Hour"):
                obsx['time'] = pd.date_range(f"{ii}-01-01", freq="H", periods=num)
            elif (self.Obs_TimRes=="Day"):
                obsx['time'] = pd.date_range(f"{ii}-01-01", freq="D", periods=num)
            elif (self.Obs_TimRes=="Month"):
                obsx['time'] = pd.date_range(f"{ii}-01-01", freq="M", periods=num)
            elif (self.Obs_TimRes=="Year"):
                obsx['time'] = pd.date_range(f"{ii}-01-01", freq="Y", periods=num)
            else:
                sys.exit(1)
            if (self.compare_Tres =="Month"):
                obsx=obsx.resample(time='1M').mean() 
            elif (self.compare_Tres =="Day"):
                obsx=obsx.resample(time='1D').mean() 
            elif (self.compare_Tres =="Hour"):
                obsx=obsx.resample(time='1H').mean() 
            else:
                sys.exit(1)   

            obsx=obsx.sel(time=slice(f'{ii}-01-01T00:00:00',f'{ii}-12-31T23:59:59'))

            # Save the output dataset to a netcdf file
            out_file = f'{self.casedir}/tmp/obs/'+f'obs_{key}_remap_{ii}.nc'
            obsx.to_netcdf(out_file)
            print(f"Done with Year {ii}")

    def make_sim_combine_parallel(self,ii):
        if self.Sim_DataGroupby == 'Year':
            VarFiles=(f'{self.Sim_Dir}/{self.Sim_Suffix}{ii}{self.Sim_Prefix}.nc')
            simx0= xr.open_dataset(VarFiles)

        else:
            VarFiles=(f'{self.Sim_Dir}/{self.Sim_Suffix}{ii}*{self.Sim_Prefix}.nc')
            print(VarFiles)
            simx0=xr.open_mfdataset(VarFiles, combine='nested',concat_dim="time",decode_times=False,chunks={'time': 30},
                                   preprocess=lambda simx0: simx0[list(self.variables.values())].astype('float64'))
            
        simx0['lon'] = (simx0['lon'] + 180) % 360 - 180.0
        lon_new = xr.DataArray(
                data=np.arange(self.Min_lon+self.res/2, self.Max_lon, self.res),
                dims=('lon',),
                coords={'lon': np.arange(self.Min_lon+self.res/2, self.Max_lon, self.res)},
                attrs={'units': 'degrees_east', 'long_name': 'longitude'}
                )
        lat_new = xr.DataArray(
                data=np.arange(self.Min_lat+self.res/2, self.Max_lat, self.res),
                dims=('lat',),
                coords={'lat': np.arange(self.Min_lat+self.res/2, self.Max_lat, self.res)},
                attrs={'units': 'degrees_north', 'long_name': 'latitude'}
                )
        new_grid = xr.Dataset({'lon': lon_new, 'lat': lat_new})
            
        if not sys.platform.startswith('win'):
            # Get the data variable
            #obsx = obsx[f'{key}'] #['fldfrc']
            # Create the regridder
            # Define the path to the weights file
            if (self.Min_lon==-180 and self.Max_lon==180 and self.Min_lat==-90 and self.Max_lat==90):
                regridder = xe.Regridder(simx0, new_grid, 'bilinear', periodic=True)
            else:
                regridder = xe.Regridder(simx0, new_grid, 'bilinear', periodic=False)
                # Perform the remapping
            simx = regridder(simx0)
        else:
            simx = simx0.interp(coords=new_grid.coords) 


        num=len(simx['time'])
        if (self.Sim_TimRes=="Hour"):
            simx['time'] = pd.date_range(f"{ii}-01-01", freq="H", periods=num)
        elif (self.Sim_TimRes=="Day"):
            simx['time'] = pd.date_range(f"{ii}-01-01", freq="D", periods=num)
        elif (self.Sim_TimRes=="Month"):
            simx['time'] = pd.date_range(f"{ii}-01-01", freq="M", periods=num)

        if (self.compare_Tres =="Month"):
            simx=simx.resample(time='1M').mean() 
        elif (self.compare_Tres =="Day"):
            simx=simx.resample(time='1D').mean() 
        elif (self.compare_Tres =="Hour"):
            simx=simx.resample(time='1H').mean() 
        else:
            sys.exit(1)     
        simx=simx.sel(time=slice(f'{ii}-01-01T00:00:00',f'{ii}-12-31T23:59:59'))
        for key in self.variables.keys():
            variable_value = self.variables[key]
            simxvar=simx[variable_value]
            # Save the output dataset to a netcdf file
            out_file = f'{self.casedir}/tmp/sim/'+f'sim_{variable_value}_remap_{ii}.nc'
            simxvar.to_netcdf(out_file)
        print(f"Done with Year {ii}")
    
    def Makefiles_parallel(self):
        print("=======================================")
        print("Create directory!")
        print(" ")
        print(" ")
        timeout=9999999
        #remove tmp directory if exist
        shutil.rmtree(f'{self.casedir}/tmp/sim',ignore_errors=True)
        shutil.rmtree(f'{self.casedir}/tmp/obs',ignore_errors=True)
        #creat tmp directory
        os.makedirs(f'{self.casedir}/tmp/sim', exist_ok=True)
        os.makedirs(f'{self.casedir}/tmp/obs', exist_ok=True)
        if self.figplot==True:
            shutil.rmtree(f'{self.casedir}/tmp/plt',ignore_errors=True)
            os.makedirs(f'{self.casedir}/tmp/plt', exist_ok=True)

        print(f"tmp directory: {self.casedir}/tmp has been created!")
        print("=======================================")
        print(" ")
        print(" ")

        minyear=self.use_Syear
        maxyear=self.use_Eyear
        #if not sys.platform.startswith('win'):
        #    num_cores = 1
        #else:
        num_cores = self.num_cores #os.cpu_count()  ##用来计算现在可以获得多少cpu核心。 也可以用multipocessing.cpu_count(),或者随意设定<=cpu核心数的数值
        #deal with observation data
        print("=======================================")
        print("deal with observation data")
        print(" ")
        print(" ")

        Parallel(n_jobs=num_cores, timeout=timeout)(delayed(self.make_obs_combine_parallel)(i) for i in range((minyear),(maxyear)+1))
        for key in self.variables.keys():
            VarFiles=(f'{self.casedir}/tmp/obs/'+f'obs_{key}_remap_*.nc')
            with xr.open_mfdataset(VarFiles, combine='nested',concat_dim="time",chunks={'time': 30}) as ds1: #,parallel=True,autoclose=True
                delayed_obj=ds1.to_netcdf(f'{self.casedir}/tmp/obs/obs_{key}.nc', compute=False)
                with ProgressBar():
                    delayed_obj.compute()
            del ds1   
        #    for ii in range((minyear),(maxyear)+1):
        #        os.remove(f'{self.casedir}/tmp/obs/'+f'obs_{key}_remap_{ii}.nc')

        print ('observation data prepared!')
        print("=======================================")
        print(" ")
        print(" ")   


        #deal with simulation data   
        print("=======================================")
        print("deal with simulation data")
        print(" ")
        print(" ") 
        #if self.Sim_GeoRes=='01min':
        #    num_cores=1


        Parallel(n_jobs=num_cores, timeout=timeout)(delayed(self.make_sim_combine_parallel)(ii) for ii in range((minyear),(maxyear)+1))
        for variable_value in self.variables.values():
            VarFiles=f'{self.casedir}/tmp/sim/'+f'sim_{variable_value}_remap_*.nc'
            with xr.open_mfdataset(VarFiles, combine='nested',concat_dim="time",chunks={'time': 30}) as ds1: #,parallel=True,autoclose=True
                delayed_obj=ds1[f'{variable_value}'].to_netcdf(f'{self.casedir}/tmp/sim/sim_{variable_value}.nc', compute=False)
                with ProgressBar():
                    delayed_obj.compute()
                del delayed_obj
            del ds1
            #for ii in range((minyear),(maxyear)+1):
            #    os.remove(f'{self.casedir}/tmp/sim/'+f'sim_{variable_value}_remap_{ii}.nc')

        print ('simulation data prepared!')
        print("=======================================")
        print(" ")
        print(" ")   
        return
    

